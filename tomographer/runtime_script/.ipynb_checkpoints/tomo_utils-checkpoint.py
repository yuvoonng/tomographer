import numpy as np
from numpy.linalg import inv, det
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, NullFormatter
from pylab import cm
import warnings
from pathlib import Path
import logging
import os
import tomllib
import tomli_w
import difflib

import hashlib
import aiohttp
from tqdm.auto import tqdm
import asyncio
import nest_asyncio

cmap = cm.RdYlBu_r
cmap.set_under('w') 
cmap.set_bad('gray') 

# ==== Check Configuration File ====
class ConfigChecker():
    """
    Load a configuration file and validate its input parameters.
    """
    def __init__(self, config):
        self.config = config 
        self.run_type = self.config.get('Reference Sample', 'run_type').lower()
        self.ref_sample = 'sdss_lss_plus'
        self.test_type = self.config.get('Test Sample', 'test_type').lower()
        self.test_data_file = self.config.get('Test Sample', 'test_data_file')
        self.test_random_file = self.config.get('Test Sample', 'test_random_file')
        self.beam_fwhm_arcmin = 1.24 # Default
        self.auto_footprint_detection = False # If True, the code will run auto_footprint() for the map
        self.footprint_definition = self.config.get('Test Sample', 'footprint_definition').lower()
        self.weight_colname = self.config.get('Test Sample', 'per_source_weight_colname')
        self.weight_path = []
        self.weight_ord = []
        self.template_cleaning = self.config.get('Test Sample Processing', 'template_cleaning').upper()
        self.filter_fwhm_deg = self.config.get('Test Sample Processing', 'filter_fwhm_deg', fallback='').lower()
        self.N_BS = config.getint('Error Estimation', 'N_BS')
        self.N_prl = 4
        self.output_filename = self.config.get('Output', 'filename')

    def check_all(self):
        self.check_reference_sample()
        self.check_test_sample()
        self.check_test_sample_processing()
        self.check_error_estimation()
        self.check_optional_input()
        return self
        
    def check_reference_sample(self):
        options = ['quick', 'fiducial']
        if self.run_type not in options:
            valid_str = " or ".join(map(str, options))
            matches = difflib.get_close_matches(self.run_type, options, n=1, cutoff=0.6)
            if matches:
                raise ValueError(f"Invalid run type. Did you mean: {matches[0]}?")
            raise ValueError(f"Invalid redshift bin number {self.zbin_num}. Must be {valid_str}.")

        if self.run_type=='quick':
            self.zbin_num = 40
        else:
            self.zbin_num = 160
        
    def check_test_sample(self):

        options = ['source_catalog', 'intensity_map']
        if self.test_type not in options:
            valid_str = " or ".join(map(str, options))
            matches = difflib.get_close_matches(self.test_type, options, n=1, cutoff=0.6)
            if matches:
                raise ValueError(f"Invalid test data type. Did you mean: {matches[0]}?")
            raise ValueError(f'Invalid test data type: {self.test_type}. Must be {valid_str}.')

        # Check test_data_file
        file = Path(self.test_data_file)
        if not file.is_file():
            raise FileNotFoundError(f"Test data file does not exist.\n"
                                    f"Is your file located here: {file.resolve()}?\n"
                                    f"If not, try specify the absolute path, or double-check the spelling.\n")

        # Check test_random_file
        if self.test_random_file:
            file = Path(self.test_random_file)
            if not file.is_file():
                raise FileNotFoundError(f"Test random file does not exist.\n"
                                    f"Is your file located here: {file.resolve()}?\n"
                                    f"If not, try specify the absolute path, or double-check the spelling.\n")

        # Check test_map_ordering
        if self.test_type == 'intensity_map':
            self.test_map_ordering = self.config.get('Test Sample', 'test_map_ordering').upper()
            options = ['NEST', 'RING']
            if self.test_map_ordering not in options:
                valid_str = " or ".join(map(str, options))
                matches = difflib.get_close_matches(self.test_map_ordering, options, n=1, cutoff=0.6)
                if matches:
                    raise ValueError(f"Invalid map ordering. Did you mean: {matches[0]}?")
                raise ValueError(f"Invalid map ordering {self.test_map_ordering}. Must be {valid_str}.")

        # Check beam_fwhm_arcmin
        if self.test_type == 'intensity_map':
            val = self.config.getfloat('Test Sample', 'beam_fwhm_arcmin')
            if val>=1.24:
                if (self.run_type=='quick') and (val>=3.):
                    raise FileNotFoundError('Quick run does not provide files with physical scale larger than 0.5 Mpc.')
                self.beam_fwhm_arcmin = val
            else:
                logging.warning(f'Beam FWHM {val} increased to 1.24 arcmin, the effective NS2048 pixel size.')
        
        # Check footprint_definition
        options = ['user_defined', 'auto_detection', 'full_sky']
        if self.test_random_file: # If test random file exists, the footprint_definition is forced to be user_defined
            self.footprint_definition = 'user_defined'
            self.auto_footprint_detection = True
            
        if self.footprint_definition not in options:
            valid_str = ", ".join(map(str, options[:-1])) + f" or {options[-1]}"
            matches = difflib.get_close_matches(self.footprint_definition, options, n=1, cutoff=0.6)
            if matches:
                raise ValueError(f"Invalid footprint definition. Did you mean: {matches[0]}?")
            raise ValueError(f"Invalid footprint definition {self.footprint_definition}. Must be {valid_str}.")
        
        if self.footprint_definition=='user_defined':
            weighting = self.config['Test Sample']
            for weight in weighting:
                if 'spatial_weight_map_file' in weight: 
                    wpath = weighting.get(weight)
                    
                    if wpath: 
                        file = Path(wpath)
                        if not file.is_file():
                            raise FileNotFoundError(f"The weighting file {wpath} does not exist.")
                        else:
                            self.weight_path.append(wpath)
        
                        ordering = 'spatial_weight_map_ordering_' + weight.split('_')[-1]
                        word = weighting.get(ordering).upper()
                        options = ['NEST', 'RING']
                        if word not in options:
                            valid_str = " or ".join(map(str, options))
                            matches = difflib.get_close_matches(word, options, n=1, cutoff=0.6)
                            if matches:
                                raise ValueError(f"Invalid map ordering. Did you mean: {matches[0]}?")
                            raise ValueError(f"Invalid map ordering {word}. Must be {valid_str}.")
                        else:
                            self.weight_ord.append(word)
        
            if len(self.weight_path)==0 and self.test_random_file=='':
                raise ValueError(f"Must have either random file or spatial weighting map for user-defined mode.")
        
        elif self.footprint_definition=='auto_detection':
            self.auto_footprint_detection = True
        
        logging.info(f'Applying {self.footprint_definition} for footprint definition...')

        # Check output filename
        if self.output_filename=='':
            raise NameError('Please at least input an output filename.')
        else:
            self.output_filename = Path(self.output_filename)

    def check_test_sample_processing(self):
        # Check template_cleaning
        if self.template_cleaning:
            options = ['HI', 'CSFD']
            if self.template_cleaning not in options:
                valid_str = " or ".join(map(str, options))
                matches = difflib.get_close_matches(self.template_cleaning, options, n=1, cutoff=0.6)
                if matches:
                    raise ValueError(f"Invalid template cleaning. Did you mean: {matches[0]}?")
                raise ValueError(f"Invalid template cleaning {self.template_cleaning}. Must be {valid_str}.")

        # Check filter_fwhm_deg
        if self.filter_fwhm_deg in ['', 'inf']:
            self.filter_fwhm_deg = np.inf
            logging.warning('Set to no high-pass filtering.')
        else:
            options = [2.0, 4.0, 8.0, 16.0]
            valid_str = ", ".join(map(str, options))
            try:
                self.filter_fwhm_deg = float(self.filter_fwhm_deg)
                if self.filter_fwhm_deg not in options:
                    raise ValueError
            except (ValueError, TypeError):
                raise ValueError(f"Invalid degree of filtering FWHM {self.filter_fwhm_deg}. Must be {valid_str} or inf.")

    def check_error_estimation(self):
        if self.N_BS<10 or self.N_BS>1000: 
            raise ValueError(f'Number of Bootstrap sampling {self.N_BS} is out of range. Set between 10 and 1000.')

    def check_optional_input(self):
        try:
            self.N_prl = self.config.getint('Optional Inputs', 'N_prl')
            if self.N_prl<=0: 
                raise ValueError(f'Number of CPUs must exceed 0.')
        except:
            logging.warning('Set number of CPUs to 4.')
            
    
# ==== File Path Handling ====

# zbin_tag = '160log1pzbins'
# ref_sample = 'sdss_lss_plus'
# scale_tag = '0.5-10Mpch'
# gamma = -0.8
# max_theta_max = 15.0

class FilePathHandler():
    """
    Manage and validate the file path for precomputed data.
    """

    def __init__(self):
        self.home_dir = None
        self.need_download = {}

    def get_filepath(self, name, **kwargs):
        """
        Retrieve the file path for a given resource and verify its existence.
        If files do not exist, they will be downloaded from Zenodo: https://zenodo.org/records/20155554.
    
        Args:
            name (str): Identifier of the file to retrieve (e.g., "ref_footprint",
                "act_data", "mask"), corresponding to keys in the internal filepath dictionary.
            **kwargs: Additional parameters used to resolve or format the filepath.
    
        Returns:
            Path: The resolved file path.
    
        Raises:
            ValueError: If the template does not exist
            FileNotFoundError: If the requested file does not exist.

        Notes:
            ``self.home_dir`` is initialized for the path of precalculated_data during installation.
            A backup mechanism is provided through ``load_or_init_precaldata()``.
        """
        ref_dir = self.home_dir / 'references'
        ref_mask_dir = ref_dir / 'footprint'
        act_dir = self.home_dir / 'activation_maps'
        mask_dir = self.home_dir / 'masks'
        temp_dir = self.home_dir / 'templates'
        wm_dir = self.home_dir / 'cosmology' / 'planck_2018'
        BS_dir = self.home_dir / 'bootstrapping'
        
        filepath = {
            'ref_footprint' : ref_mask_dir / 'sdss_ref_footprint_nest2048.fits', # SDSS footprint
            'ref_valid_pix' : ref_mask_dir / 'valid_pixels_nest2048.fits', # Valid pixels = unmasked pixels
            'ref_sample' : ref_dir / '{zbin_tag}' / '{ref_sample}_data_z_maps.pddf.parquet', # The so-called 'zero_lag' maps, all z bins in one file, for calculating b_r, no random needed
            'act_data' : act_dir / '{zbin_tag}' / '{ref_sample}_data_{scale_tag}_gamma{gamma}.pddf.parquet',
            'act_rand' : act_dir / '{zbin_tag}' / '{ref_sample}_rand_{scale_tag}_gamma{gamma}.pddf.parquet',
            'wm': wm_dir / 'wmbar_prp_{scale_tag}_gamma{gamma}_{max_theta_max}degtmax_{zbin_tag}_multi_beam_filtering.fits',
            'pre_b_r' : ref_dir / '{zbin_tag}' / 'b_ref_{ref_sample}_{scale_tag}_gamma{gamma}.fits',
            'mask': mask_dir / "{mask_name}.fits",
            'bs': BS_dir / 'N_repeat_block_ID_{bs_no}BS_default_valid_pixels_flatten.npy',
            'CSFD': temp_dir / 'CSFD_ebv_ns2048_nest.fits',
            'HI': temp_dir / 'NHI_ns2048_nest.fits'
        }

        file_on_zenodo = {
            'b_ref_sdss_lss_plus_0.5-10Mpch_gamma-0.8.fits' : 'b67b5f44beeafcf37763a99febc49c2d',
            'b_ref_sdss_lss_plus_1-10Mpch_gamma-0.8.fits' : '76ed300305c40613ad62ae14773df832',
            'b_ref_sdss_lss_plus_1.5-10Mpch_gamma-0.8.fits' : '2808a730e33a136b1d84153759c998cb',
            'b_ref_sdss_lss_plus_2-10Mpch_gamma-0.8.fits' : 'ea1dce50ec55c63a0e45dfcaa3813906',
            'b_ref_sdss_lss_plus_2.5-10Mpch_gamma-0.8.fits' : 'e3c25c91443b8d0ce5e0a526fd3868e0',
            'sdss_lss_plus_data_0.5-10Mpch_gamma-0.8.pddf.parquet' : '987374b6c4099703e8ca52c4c6aa84f4',
            'sdss_lss_plus_data_1-10Mpch_gamma-0.8.pddf.parquet' : 'ec05c270b9c7189dc6a196bef9ac6ad6',
            'sdss_lss_plus_data_1.5-10Mpch_gamma-0.8.pddf.parquet' : '90d034a45c94568ba83999b5b94505ed',
            'sdss_lss_plus_data_2-10Mpch_gamma-0.8.pddf.parquet' : '97bb8783a1f4832f8a7710dcfa621a9d',
            'sdss_lss_plus_data_2.5-10Mpch_gamma-0.8.pddf.parquet' : '8424730463d932bf73591a73a850febb',
            'sdss_lss_plus_data_z_maps.pddf.parquet' : 'dfde6dae726248401bbfe58a819351aa',
            'sdss_lss_plus_rand_0.5-10Mpch_gamma-0.8.pddf.parquet' : 'effb36c781a34a0024f39da74da85305',
            'sdss_lss_plus_rand_1-10Mpch_gamma-0.8.pddf.parquet' : 'f9b9afe0820747528bc5390fe55765c6',
            'sdss_lss_plus_rand_1.5-10Mpch_gamma-0.8.pddf.parquet' : '24acdd5442d4c66bc0365394faf505be',
            'sdss_lss_plus_rand_2-10Mpch_gamma-0.8.pddf.parquet' : 'acc8dc58a5ac0fc92ef77f0424a6b8be',
            'sdss_lss_plus_rand_2.5-10Mpch_gamma-0.8.pddf.parquet' : '31a7c34704118dd987b74489bdc31514',
            'wmbar_prp_0.5-10Mpch_gamma-0.8_15.0degtmax_160log1pzbins_multi_beam_filtering.fits' : '940b182a0e837575527f89c0643dbc7c',
            'wmbar_prp_1-10Mpch_gamma-0.8_15.0degtmax_160log1pzbins_multi_beam_filtering.fits' : 'bd6858be78546fd3abd60364d6fc3412',
            'wmbar_prp_1.5-10Mpch_gamma-0.8_15.0degtmax_160log1pzbins_multi_beam_filtering.fits' : '02a564ab7e191c05632230fff64ab1f3',
            'wmbar_prp_2-10Mpch_gamma-0.8_15.0degtmax_160log1pzbins_multi_beam_filtering.fits' : '4d6995e9780eeaec7efdc249144b646a',
            'wmbar_prp_2.5-10Mpch_gamma-0.8_15.0degtmax_160log1pzbins_multi_beam_filtering.fits' : '4f9dd5075893950241554cbf03c0dc83'        
        }
        
        path_template = filepath.get(name)
        if path_template is None:
            raise ValueError(f"Invalid mask name: {name}")
    
        # Convert to string and format with kwargs if placeholders exist
        path_str = str(path_template).format(**kwargs)
        path = Path(path_str)            

        zenodo_id = file_on_zenodo.get(path.name)
        if not path.is_file():
            if zenodo_id and (path.parent.name=='160log1pzbins' or '160log1pzbins' in path.name):
                logging.error(f'The file does not exist and will be downloaded: {path}')
                self.need_download[zenodo_id] = path
                return self
            raise FileNotFoundError(
                    f"No such file: {path.resolve()}\n"
                    f"Please first confirm if the file exists.\n"
                    f"To fix this, either move the file to that location or update the path in: "
                    f"{Path.home() / '.precal_data_path.toml'}"
                )

        if zenodo_id and (path.parent.name=='160log1pzbins' or '160log1pzbins' in path.name):
            with open(path, "rb") as f:
                file_md5 = hashlib.file_digest(f, "md5").hexdigest()
            if file_md5!=zenodo_id:
                logging.error(f'The file does not match that on Zenodo and will be downloaded: {path}')
                self.need_download[file_on_zenodo[path.name]] = path
                return self
            else: # If matches, remove from the download list to avoid repeating downloads
                self.need_download.pop(file_on_zenodo[path.name], None)
                
        return path

    ##### Downloading precalculated_data #####
    class DownloadError(Exception):
        """
        An exception that occurs while trying to download a file.
        """
        pass
    
    async def download_and_verify(self, session, expected_md5, file_path, main_pbar):
        MAX_CONCURRENT_DOWNLOADS = 5
        semaphore = asyncio.Semaphore(MAX_CONCURRENT_DOWNLOADS)

        async with semaphore:  # Ensures only 5 tasks run this block at a time
            os.makedirs(file_path.parent, exist_ok=True)
            filename = file_path.name
            
            logging.info(f"Downloading: {filename}")
            md5 = hashlib.md5()
            url = f'https://zenodo.org/records/20155554/files/{filename}'
    
            async with session.get(url, timeout=None) as response:
                # Get file size for the individual progress bar
                total_size = int(response.headers.get('content-length', 0))
            
                if response.status != 200:
                    raise self.DownloadError(f"Failed {filename}: Status {response.status}")
    
                # Create an individual progress bar for this file
                with tqdm(
                    total=total_size, 
                    unit='iB', 
                    unit_scale=True, 
                    desc=filename, 
                    leave=False # Closes the bar when finished to save screen space
                ) as file_pbar:
                    
                    with open(file_path, 'wb') as f:
                        async for chunk in response.content.iter_chunked(16384): # 16KB chunks
                            if chunk:
                                md5.update(chunk)
                                f.write(chunk)
                                file_pbar.update(len(chunk))
                                    
            clean_expected = expected_md5.lower()
            actual_md5 = md5.hexdigest().lower()
    
            if actual_md5 == clean_expected:
                main_pbar.update(1)
                return True
            raise self.DownloadError('The MD5 sum of the downloaded file is incorrect.\n'
                            + '  download: {}\n'.format(actual_md5)
                            + '  expected: {}\n'.format(clean_expected))

    async def async_download(self, file_data):
        # TCPConnector can help manage connection pooling
        connector = aiohttp.TCPConnector(limit=50)
        async with aiohttp.ClientSession(connector=connector) as session:
            with tqdm(total=len(file_data), desc="Downloading needed pre-calculated files...", unit="file") as main_pbar:
                tasks = [
                    self.download_and_verify(session, md5sum, filepath, main_pbar) 
                    for md5sum, filepath in file_data.items()
                ]
                # Run all tasks and wait for them to finish
                await asyncio.gather(*tasks)
        raise SystemExit('Downloaded sucessfully. Please re-run the code.')

    def download_all(self):
        # Download all the precalculated files on Zenodo and put them in the right directory
        zbin_tag = '160log1pzbins'
        ref_sample = 'sdss_lss_plus'
        gamma = -0.8
        for rpmin in [0.5, 1, 1.5, 2, 2.5]:
            scale_tag = f'{rpmin}-10Mpch'
            for act in ['act_data', 'act_rand']:
                self.get_filepath(act, zbin_tag=zbin_tag, ref_sample=ref_sample, scale_tag=scale_tag, gamma=gamma)
            self.get_filepath('wm', scale_tag=scale_tag, gamma=gamma, max_theta_max=15., zbin_tag=zbin_tag)
            self.get_filepath('pre_b_r', zbin_tag=zbin_tag, ref_sample=ref_sample, scale_tag=scale_tag, gamma=gamma)
        self.get_filepath('ref_sample', zbin_tag=zbin_tag, ref_sample=ref_sample)

        if len(self.need_download)!=0:
            try:
                asyncio.run(self.async_download(self.need_download))
            except:
                nest_asyncio.apply()
                asyncio.run(self.async_download(self.need_download))

            
    ##### Loading or initializing ``self.home_dir`` #####
    PRECAL_DATA_PATH = Path.home() / ".precal_data_path.toml" # Precalculated Data Path is stated at home directory
    def detect_precaldata_path(self):
        # auto-detect location
        candidates = [
            Path("./precalculated_data"),
            Path("../precalculated_data"),
            Path("../../precalculated_data"),
        ]
    
        for c in candidates:
            if c.is_dir():
                return c.resolve()
    
        raise FileNotFoundError("precalculated_data not found.")

    def init_precaldata_path(self): 
        # Initialize the precalculated file path
        data_path = self.detect_precaldata_path()
    
        config = {
            "data": {
                "precalculated_data_path": str(data_path)
            }
        }
    
        with open(FilePathHandler.PRECAL_DATA_PATH, "wb") as f:
            tomli_w.dump(config, f)
    
        logging.info(f"Initializing precalculated data path to {FilePathHandler.PRECAL_DATA_PATH}")

    def load_or_init_precaldata(self):
        if not FilePathHandler.PRECAL_DATA_PATH.exists():
            self.init_precaldata_path()
            
        with open(FilePathHandler.PRECAL_DATA_PATH, "rb") as f:
            logging.info(f"Reading precalculated data path at {FilePathHandler.PRECAL_DATA_PATH}")
            libload = tomllib.load(f)
            
            precal_data_path = Path(libload['data']['precalculated_data_path'])
            if precal_data_path != self.detect_precaldata_path():
                self.init_precaldata_path()

            self.home_dir = precal_data_path
            
# ==== Data Handling ====

def tbl2array(in_table): 
    """
    Get the 1st column of the input table, and homogenize to unitless numpy array. Properly convert bad values to nan
    
    Args:
        in_table (astropy Table): table with arbitrary number of columns
    Returns:
        array (numpy array): 1st column of the input with colname and unit stripped, properly nan-ed
    """
    return np.array(in_table[in_table.colnames[0]])


def safe_read_hmap(fname, in_ordering, out_ordering='NEST', out_nside=2048):
    """
    Read a HEALPix map of either the "proper" healpy format or a plain fits table. Data type should be float or int
    
    Args:
        fname (string): file name, including the full directory 
        in_ordering (string): HEALPix ordering of the input, "RING" or "NEST"
        out_ordering (string): HEALPix ordering of the output, "RING" or "NEST"
        out_nside (int): HEALPix nside of the output map
    Returns:
        array (numpy array): 1st column of the input with colname and unit stripped, properly nan-ed
    """
    
    # Handle HEALPix maps in both "proper" healpy map or plain fits table 
    data, hdr = hp.read_map(fname, h=True)
    hdr = dict(hdr)
    
    if 'ORDERING' in hdr:
        ftype = 'proper'
        hdr_in_ordering = hdr['ORDERING']
        in_nside = hdr['NSIDE']
    else:
        ftype = 'table' # In this case hp.read_map would have already taken only the 1st column and dropped the rest
    
    data = data.astype(float)
    data[data == hp.UNSEEN] = np.nan # Mask unseen value -1.6375e+30
    in_nside = hp.npix2nside(len(data))
    data = hp.ud_grade(data, out_nside, order_in=in_ordering, order_out=out_ordering) # includes reordering 
    
    if (ftype == 'proper'):
        if (in_ordering != hdr_in_ordering):
            warnings.warn("User specified ordering is different from the header")
    return data


def average_plus(data_lst, weights):
    """
    An error-proof weighted averaging function
    
    Args:
        data_lst (numpy array): the data vector to be taken the average
        weights (numpy array): the weights
        in_ordering (string): HEALPix ordering of the input, "RING" or "NEST"
    Returns:
        out (float): weighted average of the input data, return NaN if weights sum to zero
    """
    try:
        out = np.average(data_lst, weights=weights)
    except ZeroDivisionError: # Happens if weights sum to zero
        out = np.nan
    return out
    

# ==== Map Processing ====

def add_GLON_GLAT(cat):
    """
    Add Galactic coordinates to a catalog given the equtorial
    
    Args:
        cat (astropy Table): catalog containing ra, dec in degrees
    Returns:
        cat (astropy Table): the input catalog with Galactic l, b ('GLON', 'GLAT' in degree) added
    """
    if ('GLON' in cat.colnames) & ('GLAT' in cat.colnames):
        pass
    elif ('glon' in cat.colnames) & ('glat' in cat.colnames):
        cat.rename_column('glon', 'GLON')
        cat.rename_column('glat', 'GLAT')
    elif ('l' in cat.colnames) & ('b' in cat.colnames):
        cat.rename_column('l', 'GLON')
        cat.rename_column('b', 'GLAT')
    else:
        if 'dec' in cat.colnames:
            c_icrs = SkyCoord(ra=np.array(cat['ra'])*u.degree, dec=np.array(cat['dec'])*u.degree, frame='icrs')
        elif 'Dec' in cat.colnames:
            c_icrs = SkyCoord(ra=np.array(cat['RA'])*u.degree, dec=np.array(cat['Dec'])*u.degree, frame='icrs')
        elif 'DEC' in cat.colnames:
            c_icrs = SkyCoord(ra=np.array(cat['RA'])*u.degree, dec=np.array(cat['DEC'])*u.degree, frame='icrs')
        
        c_galactic = c_icrs.transform_to('galactic')
        cat['GLON'] = c_galactic.l.degree
        cat['GLAT'] = c_galactic.b.degree
    return 0

def add_hpid(cat, nside = 2048, nest=True):
    """
    Add HEALPix ID to a catalog, use Galactic coordinates
    
    Args:
        cat (astropy Table): catalog containing Galactic coordinate 'GLON', 'GLAT' in degrees
        nside (int): HEALPix nside
        nest (bool): use nested ording or not
    Returns:
        cat (astropy Table): the input catalog with HEALPix ID added
    """
    hpid = hp.ang2pix(nside = nside,  theta=cat['GLON'], phi =cat['GLAT'], nest=nest, lonlat=True)
    cat['HP'+str(nside)+'ID'] = hpid
    return 0

def make_healpix_map(cat, nside=2048, area_norm=False, weights=1):
    """
    Make a HEALPix count map given a source catalog with HEALPix ID for each source
    
    Args:
        cat (astropy Table): source catalog containing 'HP'+str(nside)+'ID'
        nside (int): HEALPix nside
        area_norm (bool): whether to divide the count by pixel solid angle or not
        weights (int, float, or numpy array): 1 for count map, flux for flux map, etc. 
    Returns:
        hmap (numpy array of float): HEALPix count map, flux map, or intensity map of the input catalog. If area-normalized, the unite is per sq deg
    """
    hmap   = np.zeros((hp.nside2npix(nside)))
    np.add.at(hmap, cat['HP'+str(nside)+'ID'], weights)
    if area_norm:
        pix_area = hp.nside2pixarea(nside, degrees=True)
        hmap = hmap/pix_area
    return hmap

def hp_masked_smooth(U, sm_fwhm_deg, nest=False):
    """
    Smooth a masked HEALPix map 
    
    Args:
        U (numpy array): input HEALPix map, where masked area is set to nan
        sm_fwhm_deg (float): full width half max of the smoothing beam in deg
    Returns:
        U_smoothed (numpy array): smoothed HEALPix map, masked area is set to nan
    """
    V = U.copy()
    V[U!=U] = 0
    VV = hp.smoothing(V, fwhm = np.radians(sm_fwhm_deg), nest = nest)    
    W = 0*U.copy()+1.
    W[U!=U] = 0
    WW = hp.smoothing(W, fwhm = np.radians(sm_fwhm_deg), nest = nest)    
    VV[U!=U] = np.nan
    return VV/WW


def auto_footprint(in_count_map, med_non_zero_cpp_thres = 20, thres_ratio = 0, verbose = False):
    """
    Automatically guess the footprint using the data
    Iteratively lowering the HEALPix nside (from 2048) until median non-zero count per pixel >= threshold 
    
    Args:
        in_count_map (numpy array): input HEALPix map, where masked area is set to nan; should have unit of #, not # per area
        med_non_zero_cpp_thres (int): threshold median non-zero count per pixel, default 20
        thres_ratio (float): threshold relative to ``med_non_zero_cpp_thres``; pixels with values above this threshold are discarded.
        verbose (boolean): if intermediate maps and stats during the iterations need to be shown
    Returns:
        footprint_map (boolean numpy array): footprint in HEALPix Nested nside=2048, True for pixels within the footprint 
    """
    
    nside = 2048
    new_count_map = in_count_map.copy()
    footprint_map = new_count_map > 1e-10
    med_non_zero_cpp = np.round(np.median(new_count_map[new_count_map>0.]))
    footprint_area = hp.nside2pixarea(nside, degrees=True)*len(new_count_map[new_count_map > 0])
    
    if verbose:
        hp.mollview(new_count_map, cmap=cmap, nest=True, norm='hist', title='Input Count Map')
        plt.show()
        print(f'nside = {nside} ; median non-zero count per pixel = {med_non_zero_cpp}')
        print(f'                  footprint area = {footprint_area:.1f} deg')
        if med_non_zero_cpp >= med_non_zero_cpp_thres:
            print('nside = 2048 is sufficient to make the footprint map')
    
    while med_non_zero_cpp < med_non_zero_cpp_thres:
        nside = int(nside/2)
        new_count_map = 4*hp.ud_grade(new_count_map, nside, pess=False, order_in='NEST', order_out='NEST')
        med_non_zero_cpp = np.round(np.median(new_count_map[new_count_map>0.]))
        footprint_map = new_count_map > (thres_ratio * med_non_zero_cpp) 
        footprint_area = hp.nside2pixarea(nside, degrees=True)*len(new_count_map[footprint_map])
        new_count_map_viz = new_count_map.copy()
        new_count_map_viz[~footprint_map] = np.nan
        if verbose:
            print(f'nside = {nside} ; median non-zero count per pixel = {med_non_zero_cpp}')
            print(f'                  footprint area = {footprint_area:.1f} deg')
            hp.mollview(new_count_map_viz, cmap=cmap, nest=True, norm='hist', title='Count Map Iteration')
            plt.show()
    
    # back to nside 2048
    if nside < 2048:
        footprint_map_final = footprint_map
        footprint_map = hp.ud_grade(footprint_map, 2048, order_in='NEST', order_out='NEST')
        footprint_map[in_count_map!=in_count_map] = False

    if verbose:
        hp.mollview(footprint_map, cmap=cmap, nest=True, norm='hist', title=f'Auto Footprint (Based On NS={nside} Source Count)')
        plt.show()

    return footprint_map, footprint_area
    

# ==== Redshift Binning ====

def z_binning_log1pz(z_min, z_max, n_bins):
    """
    Redshift binning scheme with the bin width scaled with log10(1+z)
    
    Args:
        z_min (float): minimum redshift
        z_max (float): maximum redshift
        n_bins (float): number of bins
    
    Returns:
        z_edges (numpy array): bin edges
        z_ctrs (numpy array): bin centers, geometric mean of the lower and upper z bounds for each bin
    """
        
    log_1pz_min = np.log10(1.+z_min)
    log_1pz_max = np.log10(1.+z_max)
    
    log_1pz_edges = np.linspace(log_1pz_min, log_1pz_max, n_bins+1)
    log_1pz_ctrs = (log_1pz_edges[1:]+log_1pz_edges[:-1])/2.
    
    z_edges = (10.**log_1pz_edges)-1.
    z_ctrs = (10.**log_1pz_ctrs)-1.
    
    return z_edges, z_ctrs


def z_rebin_tomo_out(w, new_z_edges, new_z_ctrs):
    """
    Rebin Tomographer output to new redshift binining. Within each new bin, nan in the old bins will be ignored

    Args:
        w (astropy table): raw tomographer output of the default redshift binning
        new_z_edges (numpy array): new redshift bin edges
        new_z_ctrs (numpy array): new redshift bin centers
    Returns:
        rebinned_w (astropy table): redshift rebinned tomographer output. dN/dz b or dI/dz b are inverse covariance (if ivertable) weighted or inverse variance weighted.
    """

    if 'dNdz_b' in w.colnames:
        key_col = 'dNdz_b'
    elif 'dIdz_b' in w.colnames:
        key_col = 'dIdz_b'

    w_tmp = w.copy() # adjust to have the same columns of the rebinned one
    if 'w_rr_minus' in w_tmp.colnames: # b_r_triband_correction is on
        del w_tmp['w_rr_minus']
        del w_tmp['w_rr_plus']
    w_tmp['z_eff'] = w_tmp['z']
    w_tmp = w_tmp[np.concatenate([['z', 'z_eff'], w_tmp.colnames[1:-1]]).tolist()]

    rebinned_w_holder = [[] for i in range(len(w_tmp))]

    N_BS = len(w[key_col+'_BS'][-1])
    
    for z, z_low, z_up in zip(new_z_ctrs, new_z_edges[:-1], new_z_edges[1:]):
        this_w = w_tmp[(w_tmp['z']>=z_low) & (w_tmp['z']<z_up)]

        this_ddzb = np.array(this_w[key_col])
        this_w = this_w[this_ddzb == this_ddzb] # ignore those bins with NaN dN/dz b or dI/dz b
        if len(this_w) == 1:
            for i, col in enumerate(w_tmp.colnames):
                if col[-2:] == 'BS':
                    rebinned_w_holder[i].append(np.array(this_w[col][0]))
                else:
                    rebinned_w_holder[i].append(this_w[col][0])
        else:
            cov_mtx = np.cov(np.array(this_w[key_col+'_BS']))

            if det(cov_mtx) == 0 : #If covariance matrix not invertable, keep the diagonal variance only
                cov_mtx = np.diag(np.diag(cov_mtx))

            inv_cov = inv(cov_mtx)

            identity = np.ones(len(this_w))
            weights = np.dot(inv_cov, identity) / np.dot(identity, np.dot(inv_cov, identity))
            
            for i, col in enumerate(w_tmp.colnames):
                if col=='z':
                    rebinned_w_holder[i].append(z)
                elif col=='z_low':
                    rebinned_w_holder[i].append(z_low)
                elif col=='z_up':
                    rebinned_w_holder[i].append(z_up)
                elif col[-3:] == 'err':# ignore off-diagnal terms of the covariance and just use variance for now
                    if len(this_w) == 0:
                        rebinned_w_holder[i].append(np.nan)
                    else:
                        cov = np.diag(np.power(this_w[col], 2))
                        wgt_variance = np.dot(weights, np.dot(cov, weights))
                        rebinned_w_holder[i].append(wgt_variance**0.5)
                elif col[-2:] == 'BS':# Bootstrap samples, each z append a list instead of a number
                    if len(this_w) == 0:
                        rebinned_w_holder[i].append(np.array([np.nan]*N_BS))
                    else:
                        rebinned_w_holder[i].append(np.average(this_w[col], weights = weights, axis=0))
                elif col=='z_eff':
                    if len(this_w) == 0:
                        rebinned_w_holder[i].append(np.nan)
                    else:
                        rebinned_w_holder[i].append(np.average(this_w['z'], weights = weights))
                else:
                    if len(this_w) == 0:
                        rebinned_w_holder[i].append(np.nan)
                    else:
                        rebinned_w_holder[i].append(np.average(this_w[col], weights = weights))

    rebinned_w = Table()
    for i, col in enumerate(w_tmp.colnames):
        rebinned_w[col] = rebinned_w_holder[i]

    # For those with BS samples, change the errors to STD of the BS samples. This is the proper way to propogate covariance
    for col in rebinned_w.colnames:
        if col[-2:] == 'BS':
            col_err = col.split('BS')[0] + 'err'
            rebinned_w[col_err] = np.std(np.array(rebinned_w[col]), axis=1)

    return rebinned_w
