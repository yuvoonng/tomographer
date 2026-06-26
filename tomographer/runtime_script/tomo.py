#!/usr/bin/env python
# coding: utf-8

import numpy as np
import numexpr as ne
import healpy as hp
from scipy import interpolate
import pandas as pd
from sklearn.linear_model import LinearRegression
from pylab import cm
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, NullFormatter

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

from importlib import resources
from multiprocess import Pool
from pathlib import Path
import configparser
import argparse
import asyncio
import logging
import difflib
import shutil
import psutil
import time
import sys
import os
import gc

from . import tomo_utils as utils

logging.basicConfig(
        level=logging.INFO,           
        format="%(asctime)s - %(levelname)s - %(message)s"
)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--conf', type=str, default='conf.ini', help='Path to config file')
        
    # 2. Subparsers setup
    subparsers = parser.add_subparsers(dest="command", required=False)
    
    # 'init' subcommand
    subparsers.add_parser("init", help="Initialize a new project")
    
    # 'rebin' subcommand
    rebin_parser = subparsers.add_parser("rebin", help="Rebin measurement files")
    rebin_parser.add_argument("measurement_file", type=str, help="Path to the measurement file")
    rebin_parser.add_argument("bin_number", type=int, help="The rebinning number")
    
    args = parser.parse_args()
    
    if args.command and '--conf' in sys.argv:
        parser.error("Argument '--conf' cannot be used with subcommands (init, rebin).")

    if args.command == "init":
        init(args)
    elif args.command == "rebin":
        rebin(args)
    else:
        tomo(args)

def init(args):
    filepath_handler = utils.FilePathHandler()
    filepath_handler.init_precaldata_path()
    filepath_handler.check_file_default()

    logging.info("Initialization succeeded!")

def rebin(args):
    measurement_file, bin_number = args.measurement_file, args.bin_number
    w = Table.read(measurement_file)

    if (bin_number<=0) or (bin_number>len(w)):
        raise ValueError('Rebin number must be between 0 and the original bin number.')

    rebinned_z_edges, rebinned_z_ctrs = utils.z_binning_log1pz(0., 4.2, bin_number)
    rebinned_w = utils.z_rebin_tomo_out(w, rebinned_z_edges, rebinned_z_ctrs)

    hdr = fits.Header()
    primary_hdu = fits.PrimaryHDU(header=hdr)      
    hdu_w = fits.BinTableHDU(rebinned_w) 
    hdul = fits.HDUList([primary_hdu, hdu_w])
    
    output_filename = measurement_file.removesuffix('.fits') + f'_rebin{bin_number}.fits'
    hdul.writeto(output_filename, overwrite=True)

    logging.info(f'Rebin sucessfully. Check: {output_filename}')

def tomo(args):

    start = time.time()
    
    process = psutil.Process(os.getpid())
    def log_memory():
        mem_mb = process.memory_info().rss / 1024**2
        logging.info("Memory usage: %.2f MB", mem_mb)
    
    
    from matplotlib import rcParams
    rcParams['figure.dpi'] = 150
    
    cmap = cm.RdYlBu_r
    cmap.set_under('w') 
    cmap.set_bad('gray') 
    DPI = 300
    
    
    # # Read configuration file

    config_filename = args.conf
    if not Path(config_filename).exists():
        if config_filename == 'conf.ini':
            hint = "'--conf /path/to/your/config/file' "
        else:
            hint = ""
        raise FileNotFoundError(f"Config file {config_filename} does not exist.\n"
                                f"Try {hint}to specify a different location, or double-check the spelling.")
    logging.info(f"Loading configuration file {config_filename}...")

    # ==== Initiate Config ==== 
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(f'{config_filename}')
    config_par = utils.ConfigChecker(config).check_all() # Check the compatibility of the configuration file
    
    
    if config_par.test_type == 'source_catalog':
        ddz_unit = r'$\rm[\#/deg^2]$'
        ddz_label = '(dN/dz)b '+ddz_unit
    elif config_par.test_type == 'intensity_map':
        ddz_unit = r'$\rm[map\ unit]$'
        ddz_label = '(dI/dz)b '+ddz_unit
    
    
    filepath_handler = utils.FilePathHandler()
    filepath_handler.load_or_init_precaldata() # Load/initialize precalculated data file path
    get_filepath = filepath_handler.get_filepath
    
    
    # ==== Scale ==== 
    # rp_max is fixed to 10 Mpc/h 
    # rp_min can be [0.5, 1, 1.5, 2, 2.5] Mpc/h. Needs to be larger than the beam to avoid a strong 1-halo term signal
    
    zbin_tag = '{}log1pzbins'.format(config_par.zbin_num) 
    ref_sample_filepath = get_filepath('ref_sample', zbin_tag=zbin_tag, ref_sample=config_par.ref_sample)
    
    rp_min = (config_par.beam_fwhm_arcmin//1.5)/2. 
    rp_min = min(max(rp_min, 0.5), 2.5) # (min, max) = (0.5, 2.5)
    
    if np.isclose(rp_min%1, 0):
        scale_tag = str(int(rp_min))+'-10Mpch'
    else:
        scale_tag = str(rp_min)+'-10Mpch'
    print('(rp_min, rp_max) = {}'.format(scale_tag))
    
    max_theta_max = 15. # Deg, maximum angular scale used, which overwrites rp_max at the lowest z
    gamma = -0.8 # Scale weighting power index
    
    act_data_filepath = get_filepath('act_data', zbin_tag=zbin_tag, ref_sample=config_par.ref_sample, scale_tag=scale_tag, gamma=gamma)
    act_rand_filepath = get_filepath('act_rand', zbin_tag=zbin_tag, ref_sample=config_par.ref_sample, scale_tag=scale_tag, gamma=gamma)
    
    
    if config_par.template_cleaning:
        template_filepath = get_filepath(config_par.template_cleaning)
    
    if config_par.N_BS<=100: max_N_BS = 100
    elif config_par.N_BS<=1000: max_N_BS = 1000
    bootstrap_filepath = get_filepath('bs', bs_no=max_N_BS)
    
    ref_valid_pix_filepath = get_filepath('ref_valid_pix')
    ref_footprint_filepath = get_filepath('ref_footprint')
    
    matter_filepath = get_filepath('wm', scale_tag=scale_tag, gamma=gamma, max_theta_max=max_theta_max, zbin_tag=zbin_tag)
    
    pre_b_r_filepath = get_filepath('pre_b_r', zbin_tag=zbin_tag, ref_sample=config_par.ref_sample, scale_tag=scale_tag, gamma=gamma)
    
    
    if len(filepath_handler.need_download)!=0:
        asyncio.run(filepath_handler.async_download(filepath_handler.need_download))
    
    
    # # Test data
    
    from astropy.cosmology import FlatLambdaCDM
    
    h      = 0.6737
    rp_min = 0.5/h # Mpc
    rp_max = 10./h # Mpc
    max_theta_max = 15. # deg
    cosmo = FlatLambdaCDM(H0=h*100., Om0=0.3132, Tcmb0=2.725) # Planck18
    
    z_edges, z_ctrs = utils.z_binning_log1pz(0., 4.2, config_par.zbin_num)
    delta_z = z_edges[1:] - z_edges[:-1]
    print('dense z bin centers [{} log(1+z) bins in 0<z<4.2]:\n'.format(config_par.zbin_num), z_ctrs, '\n')
    
    if config_par.zbin_num==160:
        rebinned_z_edges, rebinned_z_ctrs = utils.z_binning_log1pz(0., 4.2, 40)
        print('rebinned z bin centers [40 log(1+z) bins in 0<z<4.2]:\n', rebinned_z_ctrs, '\n')
    
    
    logging.info('Reading test data file: {}...'.format(config_par.test_data_file))

    try:
        if config_par.test_type == 'source_catalog':    
            test_data = Table.read(config_par.test_data_file)
            logging.info('Converting test data to HEALPix map...')
            utils.add_GLON_GLAT(test_data)
            utils.add_hpid(test_data, 2048)
            if config_par.weight_colname:
                test_data_map = utils.make_healpix_map(test_data, 2048, weights=test_data[config_par.weight_colname]/np.nanmedian(test_data[config_par.weight_colname]))
            else:
                test_data_map = utils.make_healpix_map(test_data, 2048)
            
            if config_par.test_random_file != '':
                logging.info('Reading test random file: {}...'.format(config_par.test_random_file))
                if config_par.test_type == 'source_catalog':
                    test_rand = Table.read(config_par.test_random_file)
                else:
                    test_rand_map = utils.safe_read_hmap(config_par.test_random_file, config_par.test_map_ordering)
            
                logging.info('Converting test random to HEALPix map...')
                utils.add_GLON_GLAT(test_rand)
                utils.add_hpid(test_rand, 2048)
                
                if config_par.weight_colname:
                    test_rand_map = utils.make_healpix_map(test_rand, 2048, weights=test_rand[config_par.weight_colname]/np.nanmedian(test_rand[config_par.weight_colname]))
                else:
                    test_rand_map = utils.make_healpix_map(test_rand, 2048)
            
                del test_rand
            else: logging.warning('No random file.')
        
        else:
            test_data_map = utils.safe_read_hmap(config_par.test_data_file, config_par.test_map_ordering)
            if config_par.test_random_file != '':
                logging.info('Reading test random file: {}...'.format(config_par.test_random_file))
                test_rand_map = utils.safe_read_hmap(config_par.test_random_file, config_par.test_map_ordering)

        pass

    except Exception as e:
        raise ValueError("Failed to read the test file. Please ensure that the file is in the correct format and matches the `test_type`.") from e

    logging.info('Successfully read test file(s).')
    
    
    # ## Footprint selection
    
    med_non_zero_cpp_thres = 10 # 5-20 recommended, threshold for non zero count per pixel in the utils.auto_footprint function  
    
    if config_par.test_type == 'source_catalog':
        if config_par.auto_footprint_detection:
            if config_par.test_random_file:
                footprint_map = test_rand_map
            else:
                footprint_map = test_data_map
                logging.info('Auto detecting footprint...')
            test_footprint_map, total_area = utils.auto_footprint(footprint_map
                                           , med_non_zero_cpp_thres = med_non_zero_cpp_thres
                                           , thres_ratio=.0, verbose = False)
        else:
            test_footprint_map = np.ones(12*2048*2048).astype(bool)
    
    elif config_par.test_type == 'intensity_map': # Treat NaN or hp.UNSEEN (-1.6375e+30) as out of footprint        
        test_footprint_map = (test_data_map == test_data_map) & (test_data_map != hp.UNSEEN)
        if config_par.test_random_file != '': # Join non-NaN area in both the random and data
            test_footprint_map = test_footprint_map & ((test_rand_map == test_rand_map) & (test_rand_map != hp.UNSEEN))
    
        if config_par.auto_footprint_detection:
            uniq_test = len(np.unique(test_data_map[test_data_map>=0]))/len(test_data_map[test_data_map>=0])<0.05
            zero_test = len(test_data_map[test_data_map==0])/len(test_data_map[test_data_map==test_data_map])>0.01
            # If an intensity map pass either uniq_test or zero_test, it is treated as a discrete count map.
            # Its footprint can then be defined with auto_detection.
            # Otherwise it can only be applied with no_selection
            if uniq_test or zero_test:
                logging.info('It is a discrete count map. Auto detecting footprint...')
                test_footprint_map, total_area = utils.auto_footprint(test_data_map
                                   , med_non_zero_cpp_thres = med_non_zero_cpp_thres
                                   , thres_ratio=.0, verbose = False)
            else:
                logging.info('It is not a discrete count map. Applying no footprint selection...')
                footprint_definition = 'no_selection'
    
    
    # Apply to data/random
    test_data_map[~test_footprint_map] = np.nan
    if config_par.test_random_file != '':
        test_rand_map[~test_footprint_map] = np.nan
        
    
    if config_par.test_type == 'source_catalog':
        pix_area = hp.nside2pixarea(2048, degrees=True)
    
        test_data_map = test_data_map/pix_area # to unit of count per square degree
        if config_par.test_random_file != '':
            test_rand_map = test_rand_map/pix_area # to unit of count per square degree
    
    
    
    # ## Masking & Weighting
    
    default_mask_option = ['CSFD_cosmology_area'
        , 'fsky_25_low_dust'
        , 'fsky_50_low_dust'
        , 'fsky_75_low_dust'
        , 'globular_clusters_veto'
        , 'LMC_SMC_veto'
        , 'nearby_galaxy_clusters_veto'
        , 'Planck_point_source_veto'
        , 'SDSS_veto'
        , 'rec_cut']
    
    mask_map = []
    for mask_name in default_mask_option:
        optional_input = config.getboolean('Optional Inputs', mask_name)
        if optional_input: 
            mask_map.append(mask_name)
    
    
    test_weights_map = np.ones(12*2048*2048)
    
    for map_name in mask_map:
    
        if map_name == 'rec_cut':
            logging.info(f'Applying rectangular cut...')
    
            rec_min_lat = config.getfloat('Optional Inputs', 'rec_min_lat')
            rec_max_lat = config.getfloat('Optional Inputs', 'rec_max_lat')
            rec_min_lon = config.getfloat('Optional Inputs', 'rec_min_lon')
            rec_max_lon = config.getfloat('Optional Inputs', 'rec_max_lon')
    
            rec_cut_coord = config.get('Optional Inputs', 'rec_cut_coord').lower()
            options = ['galactic', 'ecliptic', 'equatorial']
            if rec_cut_coord not in options:
                matches = difflib.get_close_matches(rec_cut_coord, options, n=1, cutoff=0.6)
                if matches:
                    raise ValueError(f"Invalid coordinate frame. Did you mean: {matches[0]}?")
                raise ValueError(f'Invalid coordinate frame: {rec_cut_coord}')
                
            if rec_cut_coord == 'galactic':
                rec_frame = 'galactic'
            elif rec_cut_coord == 'ecliptic':
                rec_frame = 'geocentricmeanecliptic'
            elif rec_cut_coord == 'equatorial':
                rec_frame = 'icrs'
    
            coord_min = (rec_min_lon, rec_min_lat)
            coord_max = (rec_max_lon, rec_max_lat)
            if coord_min>coord_max: 
                raise ValueError("Minimum coordinate(s) cannot be greater than maximum coordinate(s).")
            coord_min, coord_max = [SkyCoord(x*u.deg, y*u.deg, frame=rec_frame).transform_to('galactic')
                  for x, y in [coord_min, coord_max]]
            
            theta_top = np.radians(90 - coord_max.b.degree)
            theta_bottom = np.radians(90 - coord_min.b.degree)
            pix_in_strip = hp.query_strip(2048, theta_top, theta_bottom, nest=True)
        
            _, phi = hp.pix2ang(2048, pix_in_strip, nest=True)
            lon_deg = np.degrees(phi)
            lon_mask = (lon_deg >= coord_min.l.degree) & (lon_deg <= coord_max.l.degree)
            
            rec_pix = pix_in_strip[lon_mask]
    
            rec_include = config.getboolean('Optional Inputs', 'rec_include')
            if rec_include == True:
                rec_mask_map = np.zeros(12*2048*2048)
                rec_mask_map[rec_pix] = 1.
            else:
                rec_mask_map = np.ones(12*2048*2048)
                rec_mask_map[rec_pix] = 0.
            
            test_weights_map *= rec_mask_map
            
        else:
            logging.info(f'Reading {map_name} for masking...')
            path = get_filepath('mask', mask_name = map_name)
    
            mask_data = Table.read(path)
            mask_data = utils.tbl2array(mask_data)
            
            test_weights_map *= mask_data
            del mask_data
    
    
    if len(config_par.weight_path)!=0:
        for wpath, word in zip(config_par.weight_path, config_par.weight_ord):
            logging.info(f'Reading spatial weighting file {wpath} ...')
            weight_data = utils.safe_read_hmap(wpath, word)
            test_weights_map *= weight_data
            del weight_data
    
    
    # ## Footprint
    
    ref_valid_pix = Table.read(ref_valid_pix_filepath) # Unmasked, valid pixels in SDSS footpint(Healpix NS2048)
    
    ref_footprint = Table.read(ref_footprint_filepath)
    ref_footprint = utils.tbl2array(ref_footprint)
    
    
    test_weights_map *= ref_footprint.astype(float)
    test_weights_map *= test_footprint_map.astype(float)
    
    # set zero for NaN weights, otherwise booling it would get True
    test_weights_map[test_weights_map!=test_weights_map] = 0.
    
    
    # ==== Remove Zero-Weight and Out of Reference/Test Footprint Areas From the Data Vector ====
    joint_footprint = test_weights_map.astype(bool)
    sel_joint_in_ref_vpix_space = joint_footprint[ref_valid_pix['HP2048ID']]
    joint_valid_pix = ref_valid_pix[sel_joint_in_ref_vpix_space]
    
    test_data_map[~joint_footprint] = np.nan
    test_data_map_vpix = test_data_map[joint_valid_pix['HP2048ID']]
    if config_par.test_random_file != '':
        test_rand_map[~joint_footprint] = np.nan
        test_rand_map_vpix = test_rand_map[joint_valid_pix['HP2048ID']]
        
    test_weights_map_vpix = test_weights_map[joint_valid_pix['HP2048ID']]
    if len(test_weights_map_vpix)==0:
        raise ValueError('No overlapping area between test data and reference sample.')
    
    
    # Summary tag (Boolean) to indicate if spatial weights are contentious for valid pixels, not just all 1's for 'unmasked'
    spatial_weighting = ~((np.min(test_weights_map_vpix) == 1.) & (np.max(test_weights_map_vpix) == 1.))
    
    
    # Mean match the test random to the data within unmasked area (include spatial weights, if given)
    if config_par.test_random_file != '':
        
        test_rand_map *= (test_data_map_vpix @ test_weights_map_vpix)/(test_rand_map_vpix @ test_weights_map_vpix)
        test_rand_map_vpix = test_rand_map[joint_valid_pix['HP2048ID']]
    
    
    total_area  = len(test_data_map_vpix) * hp.nside2pixarea(2048, degrees=True)
    logging.info('Overlapping pixels with the reference sample: {}/{}'.format(len(joint_valid_pix), 12*2048**2))
    logging.info('Overlapping area with the reference sample: {:.2f} deg^2'.format(total_area))
    
    if config_par.test_type == 'source_catalog':
        mean_density = np.nansum(test_data_map_vpix) * pix_area / total_area
        logging.info('Overlapping sources with the reference sample: {}/{}'.format(int(np.nansum(test_data_map_vpix) * pix_area), len(test_data)))
        logging.info('Overlapping density with the reference sample: {:.2f} /deg^2'.format(mean_density))    
    
    
    # ## Template-Based Foreground Mitigation
    
    # ==== HI or CSFD Cleaning ===
    # Regress out test data that can be linear predicted by HI column density in HI4PI or E(B-V) in CSFD, in unmaksed areas with non-zero weights
    if config_par.template_cleaning != '':
        logging.info(f'Reading {config_par.template_cleaning} for foreground mitigation')
        temp = Table.read(template_filepath)
        temp = temp[temp.colnames[0]] 
        
        x = np.array(temp)
        y = test_data_map
        
        x = x[y==y]# Apply join mask to temp 
        y = y[y==y]
        x = x.reshape((-1, 1))
        
        model = LinearRegression()
        model.fit(x, y)
        
        xx = np.linspace(np.min(x), np.max(x), 10)
        yy = model.intercept_ + model.coef_ * xx
        
        # Plot temp vs test_map values and the best fit linear model
        x_low = np.nanpercentile(x, 0.1)/1.4
        x_up = np.nanpercentile(x, 99.9)*1.4
        y_low = np.nanpercentile(y, 0.1)/1.4
        y_up = np.nanpercentile(y, 99.9)*1.4
        
        # Residual fluctuation after correction
        y_pred = model.predict(x)
        y_corr = y - y_pred
        
        # Put correction back to test_data_map
        temp_pred_test_data_map = np.array([np.nan]*len(test_data_map))
        temp_pred_test_data_map[test_data_map == test_data_map] = y_pred
        test_data_map = test_data_map - temp_pred_test_data_map
        test_data_map_vpix = test_data_map[joint_valid_pix['HP2048ID']]
        
        if config_par.test_random_file != '':
            test_rand_map = test_rand_map - temp_pred_test_data_map
            test_rand_map_vpix = test_rand_map[joint_valid_pix['HP2048ID']]
            
        # delete unused objects
        del x, y, model, y_pred, y_corr, temp_pred_test_data_map
    
    
    # ## High-pass Filtering
    
    # ==== High-Pass Filtering for the Test Sample ====
    if config_par.filter_fwhm_deg != np.inf:
        logging.info('Smoothing test data map with {} degree'.format(config_par.filter_fwhm_deg))
        test_data_map_smoothed = utils.hp_masked_smooth(test_data_map, config_par.filter_fwhm_deg, nest=True)
        test_data_map -= test_data_map_smoothed
        test_data_map[(test_data_map!=test_data_map) | (test_data_map == -1.6375e+30)] = np.nan
        test_data_map_vpix = test_data_map[joint_valid_pix['HP2048ID']]
        print(f'test data mean = {np.nanmean(test_data_map):.5f}; masked test data mean = {np.nanmean(test_data_map_vpix):.5f}')
        
        del test_data_map_smoothed
    
        if config_par.test_random_file != '':
            logging.info('Smoothing test random map')
            test_rand_map_smoothed = utils.hp_masked_smooth(test_rand_map, config_par.filter_fwhm_deg, nest=True)
            test_rand_map -= test_rand_map_smoothed
            test_rand_map[(test_rand_map!=test_rand_map) | (test_rand_map == -1.6375e+30)] = np.nan
            test_rand_map_vpix = test_rand_map[joint_valid_pix['HP2048ID']]
            print(f'test random mean = {np.nanmean(test_rand_map):.5f}; masked test random mean = {np.nanmean(test_rand_map_vpix):.5f}')
            
            del test_rand_map_smoothed
    
    
    # # Preparation for correlation analysis
    
    # ## Bootstrapping
    
    N_repeat_block_ID_BS = np.load(bootstrap_filepath) # Needs to be reshaped
    
    # ==== Slice Pre-Calculated Block-Boostrapping Pixel Realizations ====
    
    max_blk_id = int(len(N_repeat_block_ID_BS)/max_N_BS)
    N_repeat_block_ID_BS = N_repeat_block_ID_BS.reshape(max_N_BS, max_blk_id)
    N_repeat_block_ID_BS = N_repeat_block_ID_BS[:config_par.N_BS] # Chop to the config_par.N_BS specified
    N_repeat_block_ID_BS = N_repeat_block_ID_BS.astype(np.ushort) # BS weights take a lot of memory, attempt to reduce it
    
    BS_weights_BS = [] # BS weights for each final valid pixel (masked referece, test footprints + joint veto)
    for i_BS, N_repeat_block_ID in enumerate(N_repeat_block_ID_BS):    
        block = N_repeat_block_ID[joint_valid_pix['HP16ID']]
        BS_weights_BS.append(block)
    BS_weights_BS = np.array(BS_weights_BS)
    
    del N_repeat_block_ID_BS, N_repeat_block_ID
    print('Bootstrapping matrix shape: {}'.format(BS_weights_BS.shape))
    
    
    gc.collect()
    log_memory()
    
    
    # # Correlation analysis
    
    # ## Auto- & cross-correlation
    
    def get_w(i_z, count_pair, spatial_weighting, cal_err=False, **kwargs):
        """
        Compute the correlation function (w), with optional spatial weighting.
        
        Args:
            i_z (int): Index of the redshift slice.
            count_pair (dict): Pair counts, expected keys include
                "ref_data", "test_data", "test_random", "act_data", and "act_random".
            spatial_weighting (bool): If True, apply spatial weights.
            cal_err (bool): If True, estimate errors via bootstrapping.
            **kwargs: Additional parameters (e.g., weights).
        
        Returns:
            tuple:
                For auto-correlation: (w, dd, dr)
                For cross-correlation: (w, dd, dr, rd, rr)
        """
    
        def auto_est(dd, dr, rd=None, rr=None): 
            """
            Estimate cross-correlation using Landy–Szalay (if rd, rr provided)
            or Davis–Peebles (otherwise).
            """
            if rd is not None and rr is not None:
                return (dd - dr - rd + rr) / rr
            return dd/dr - 1.
    
        def cross_est(dd, dr, rd=None, rr=None):
            """
            Estimate cross-correlation using Landy–Szalay-like (if rd, rr provided)
            or Davis–Peebles-like (otherwise).
            """
            if rd is not None and rr is not None:
                return (dd - dr - rd + rr)
            return (dd - dr)
    
        corr_type = 'auto' if any('ref' in k for k in count_pair) else 'cross' # if 'ref' in count_pair, then it is calculation of auto-corr
        est_func = cross_est if corr_type == 'cross' else auto_est
    
        act_data = count_pair.get('act_data')
        act_rand = count_pair.get('act_rand')
        data = count_pair.get('ref_data', count_pair.get('test_data'))
        rand = count_pair.get('ref_rand', count_pair.get('test_rand'))
    
        dd = ne.evaluate("data * act_data")
        dr = ne.evaluate("data * act_rand")
        dd_nonzeroterm = np.where(dd!=0)[0] # select non-zero term
        dr_nonzeroterm = np.where(dr!=0)[0]
    
        dd_nonzero = dd[dd_nonzeroterm]
        dr_nonzero = dr[dr_nonzeroterm]
        
        if cal_err:
            bs_nonzero = BS_weights_BS[:, dd_nonzeroterm]
            dd_nonzero = ne.evaluate("dd_nonzero * bs_nonzero") # put "weight" for bootstrapping
            bs_nonzero = BS_weights_BS[:, dr_nonzeroterm]
            dr_nonzero = ne.evaluate("dr_nonzero * bs_nonzero")
        else:
            dd_nonzero = np.atleast_2d(dd_nonzero) # Convert to 2D format for compatibility with bootstrapping
            dr_nonzero = np.atleast_2d(dr_nonzero)
        
        if rand is not None:
            rd = ne.evaluate("rand * act_data")
            rr = ne.evaluate("rand * act_rand")
            rd_nonzeroterm = np.where(rd!=0)[0]
            rr_nonzeroterm = np.where(rr!=0)[0]
    
            rd_nonzero = rd[rd_nonzeroterm]
            rr_nonzero = rr[rr_nonzeroterm]
    
            if cal_err:
                bs_nonzero = BS_weights_BS[:, rd_nonzeroterm]
                rd_nonzero = ne.evaluate("rd_nonzero * bs_nonzero")
                bs_nonzero = BS_weights_BS[:, rr_nonzeroterm]
                rr_nonzero = ne.evaluate("rr_nonzero * bs_nonzero")
                del bs_nonzero
            else:
                rd_nonzero = np.atleast_2d(rd_nonzero)
                rr_nonzero = np.atleast_2d(rr_nonzero)
                
        # Compute dd, dr, rd and rr
        dd0 = np.sum(dd_nonzero, axis=1, dtype=np.float32)
        dr0 = np.sum(dr_nonzero, axis=1, dtype=np.float32)
    
        if corr_type == 'cross': 
            ave_data_weight = np.sum(act_data)
            ave_rand_weight = np.sum(act_rand)
            dd0 /= ave_data_weight
            dr0 /= ave_rand_weight
            
        return_pair = (dd0.squeeze(), dr0.squeeze())
        
        if rand is not None:
            rd0 = np.sum(rd_nonzero, axis=1, dtype=np.float32)
            rr0 = np.sum(rr_nonzero, axis=1, dtype=np.float32)
            if corr_type == 'cross': 
                rd0 /= ave_data_weight
                rr0 /= ave_rand_weight
            return_pair += (rd0.squeeze(), rr0.squeeze())
        else:
            rd0 = None
            rr0 = None
    
        w_no_wgt = est_func(dd0, dr0, rd0, rr0)
        
        if not spatial_weighting:
            gc.collect()
            return w_no_wgt.squeeze(), return_pair
    
    
        # If spatial_weighting is True, apply weighting to dd, dr, rd and rr
        weights = kwargs.get('weights')
        w_nonzero = weights[dd_nonzeroterm]
        dd_w = ne.evaluate("dd_nonzero * w_nonzero")
        w_nonzero = weights[dr_nonzeroterm]
        dr_w = ne.evaluate('dr_nonzero * w_nonzero')
        if corr_type == 'cross':
            ave_data_weight = np.sum(act_data * weights)
            ave_rand_weight = np.sum(act_rand * weights)
            dd_w /= ave_data_weight
            dr_w /= ave_rand_weight
    
        dd0 = np.sum(dd_w, axis=1, dtype=np.float32)
        dr0 = np.sum(dr_w, axis=1, dtype=np.float32)
        return_pair = (dd0.squeeze(), dr0.squeeze())
    
        if rand is not None:
            w_nonzero = weights[rd_nonzeroterm]
            rd_w = ne.evaluate('rd_nonzero * w_nonzero')
            w_nonzero = weights[rr_nonzeroterm]
            rr_w = ne.evaluate('rr_nonzero * w_nonzero')
            if corr_type == 'cross': 
                rd_w /= ave_data_weight
                rr_w /= ave_rand_weight
            
            rd0 = np.sum(rd_w, axis=1, dtype=np.float32)
            rr0 = np.sum(rr_w, axis=1, dtype=np.float32)
            return_pair += (rd0.squeeze(), rr0.squeeze())
        else:
            rd0 = None
            rr0 = None
    
        w_wgt = est_func(dd0, dr0, rd0, rr0)
    
        if corr_type == "auto":
            w = (w_wgt**2 / w_no_wgt) 
        else:
            w = w_wgt
    
        return (w.squeeze(), return_pair)
    
    
    def auto_corr(i_z, count_pair, **kwargs):
    
        # ==== Reference Auto Correlation to Get b_r ====
        ref_data_z_maps = pd.read_parquet(ref_sample_filepath, columns=[str(i_z)])
        ref_data_z_maps = ref_data_z_maps[sel_joint_in_ref_vpix_space][f'{i_z}'].to_numpy()
    
        count_pair_auto = count_pair.copy()
        count_pair_auto['ref_data'] = ref_data_z_maps
        extra_args = kwargs
        w_auto, w_pair = get_w(i_z, count_pair_auto, spatial_weighting, **extra_args)
    
        if calculate_br_error:
            w_BS_lst, _ = get_w(i_z, count_pair_auto, spatial_weighting, cal_err=True, **extra_args)
            w_auto_err = np.std(w_BS_lst)
            
        # For b_r triband correction (i.e., cross-amplitudes leaked to nearby z bins, matter when z bins are narrow)
        if b_r_triband_correction:
            # i-1 zbin cross i zbin (note: activation stays at the ith zbin)
            if i_z == 0:
                w_auto_minus = np.nan
            else:
                ref_data_z_maps = pd.read_parquet(ref_sample_filepath, columns=[str(i_z-1)])
                ref_data_z_maps = ref_data_z_maps[sel_joint_in_ref_vpix_space][f'{i_z-1}'].to_numpy()
                count_pair_auto['ref_data'] = ref_data_z_maps
                w_auto_minus, _ = get_w(i_z, count_pair_auto, spatial_weighting, **extra_args)
            
            # i+1 zbin cross i zbin (note: activation stays at the ith zbin)
            if i_z == (len(z_ctrs)-1):
                w_auto_plus = np.nan
            else:
                ref_data_z_maps = pd.read_parquet(ref_sample_filepath, columns=[str(i_z+1)])
                ref_data_z_maps = ref_data_z_maps[sel_joint_in_ref_vpix_space][f'{i_z+1}'].to_numpy()
                count_pair_auto['ref_data'] = ref_data_z_maps
                w_auto_plus, _ = get_w(i_z, count_pair_auto, spatial_weighting, **extra_args)
    
        return_data = (w_auto, )
        if b_r_triband_correction:
            return_data += (w_auto_minus, w_auto_plus, )
        if calculate_br_error:
            return_data += (w_auto_err, )
    
        del ref_data_z_maps, count_pair_auto
        gc.collect()
        return return_data
    
    
    def cross_corr(i_z, count_pair, **kwargs):    
    
        count_pair_cross = count_pair.copy()
        count_pair_cross['test_data'] = test_data_map_vpix
        if config_par.test_random_file:
            count_pair_cross['test_rand'] = test_rand_map_vpix
        extra_args = kwargs
    
        w_cross, w_pair = get_w(i_z, count_pair_cross, spatial_weighting, **extra_args)
        return_data = (w_cross, ) + w_pair
    
        w_BS_lst, _ = get_w(i_z, count_pair_cross, spatial_weighting, cal_err=True, **extra_args)
        w_cross_err = np.std(w_BS_lst)
        return_data += (w_cross_err,)
    
        del count_pair_cross
        gc.collect()
        return return_data, w_BS_lst
    
    
    def calc_corr(i_z):
        act_ref_data_z = pd.read_parquet(act_data_filepath, columns=[str(i_z)])
        act_ref_data_z = act_ref_data_z[sel_joint_in_ref_vpix_space][f'{i_z}'].to_numpy()
    
        act_ref_rand_z = pd.read_parquet(act_rand_filepath, columns=[str(i_z)])
        act_ref_rand_z = act_ref_rand_z[sel_joint_in_ref_vpix_space][f'{i_z}'].to_numpy()
    
        count_pair = {'act_data': act_ref_data_z, 'act_rand': act_ref_rand_z}
        extra_args = {"weights": test_weights_map_vpix} if spatial_weighting else {}
    
        ac = auto_corr(i_z, count_pair, **extra_args)
        cc = cross_corr(i_z, count_pair, **extra_args)
    
        del act_ref_data_z, act_ref_rand_z
        gc.collect()
        return ac, *cc
    
    
    # ==== Clustering Redshift Options ====
    def w_tr_rr_each_z(i_z):
        #  ==== Preparation ====
        z = z_ctrs[i_z]
        z_low = z_edges[i_z]
        z_up = z_edges[i_z+1]
        dz = delta_z[i_z]
        print(f'i = {i_z:<3} | z = {z:8.6f}')
        return_data = (z, z_low, z_up, dz)
    
        # Combine marked bootstrapping weights with pre-calculated activation maps
        auto_data, cross_data, w_BS_lst = calc_corr(i_z)
        return_data += (auto_data + cross_data)
        return (return_data, w_BS_lst)
    
    
    b_r_triband_correction = False
    if config_par.zbin_num==160:
        b_r_triband_correction = True # When dz for redshift binning is small, ref auto correlation leaks to the z bins in front and behind
    calculate_br_error = True
    
    
    start = time.time()
    i_z_lst = range(len(z_ctrs))
    
    with Pool(config_par.N_prl) as p:
        results = p.map(w_tr_rr_each_z, i_z_lst)
    w_tr_rr_data, dN_tr_BS = map(list, zip(*results))
    
    print('Take {:.2f} minutes to calculate the correlation functions.'.format((time.time()-start)/60))
    
    
    # # Save files
    
    # ## Matter
    
    wm = Table.read(matter_filepath)
    
    
    # ==== Beam and Filtering Info for w_m Correction ====
    beam_tag = 'NS2048' # The minimum would be the pixel window
    w_m_beam_fwhm_arcmin_lst = np.array([2., 3., 4., 4.3, 5., 5.5, 6., 7., 7.1, 8., 9., 10., 12., 15., 20., 25., 30.]) # implemented in w_m
    if config_par.beam_fwhm_arcmin > 1.5:
        i_pick = np.argmin(np.fabs(config_par.beam_fwhm_arcmin-w_m_beam_fwhm_arcmin_lst))
        w_m_beam = w_m_beam_fwhm_arcmin_lst[i_pick]
        beam_tag = f'{w_m_beam}_arcmin_beam'
    filtering_tag = f'{config_par.filter_fwhm_deg}_deg_filtered'
    
    
    # ==== Matter Clustering Correcttion Preparation ====
    # Ref-ref auto is effectively pixalized reference map cross reference catalog (pixel-window x point sources)
    f_wm_ref_ref = interpolate.interp1d(wm['z'], wm['wmbar_NS2048']) # For ref-ref auto and b_ref calculation
    # Ref-test cross is effectively pixalized test map (sometimes beam-smeared) cross reference catalog (pixel or beam window x point sources)
    f_wm_ref_test = interpolate.interp1d(wm['z'], wm[f'wmbar_{beam_tag}_{filtering_tag}'])
    
    
    w_tr_rr_data = np.array(w_tr_rr_data).transpose()
    
    i_next = 5
    
    w = Table()
    w['z'] = list(w_tr_rr_data[0])
    w['z_low'] = list(w_tr_rr_data[1])
    w['z_up'] = list(w_tr_rr_data[2])
    w['dz'] = list(w_tr_rr_data[3])
    
    w['w_rr'] = list(w_tr_rr_data[4])
    if b_r_triband_correction:
        w['w_rr_minus'] = list(w_tr_rr_data[i_next])
        w['w_rr_plus'] = list(w_tr_rr_data[i_next+1])
        i_next += 2
    if calculate_br_error:
        w['w_rr_err'] = list(w_tr_rr_data[i_next])
        i_next += 1
    
    w['w_tr'] = list(w_tr_rr_data[i_next]) 
    w['DD_tr']  = list(w_tr_rr_data[i_next+1])
    w['DR_tr']  = list(w_tr_rr_data[i_next+2])
    i_next += 2
    if config_par.test_random_file != '':
        w['RD_tr']  = list(w_tr_rr_data[i_next+1])
        w['RR_tr']  = list(w_tr_rr_data[i_next+2])
        i_next += 2
    w['w_tr_err'] = list(w_tr_rr_data[i_next+1])
    w['w_tr_BS'] = dN_tr_BS
    
    w['w_m_ref_auto'] = list(f_wm_ref_ref(w['z']))
    w['w_m_cross'] = list(f_wm_ref_test(w['z']))
    
    # # Set w_m to be nan if it is lower than 0.0003 (some are negative at very low z due to high-pass filtering)
    w['w_m_cross'][w['w_m_cross']<0.0003] = np.nan
    
    w.sort('z')
    
    
    # ## Bias correction
    
    # Correct for the formalism difference (see Appendix A of Chiang+26)
    
    deg_per_Mpc = cosmo.arcsec_per_kpc_proper(z_ctrs).to(u.deg/u.Mpc).value
    theta_min = deg_per_Mpc * rp_min
    theta_max = deg_per_Mpc * rp_max
    theta_max[theta_max>max_theta_max] = max_theta_max
    rmax, rmin = np.radians(theta_max), np.radians(theta_min)
    X = (rmax**(gamma+1) - rmin**(gamma+1)) / (gamma+1)
    
    w['w_rr'] *= X
    if b_r_triband_correction:
        w['w_rr_minus'] *= X
        w['w_rr_plus'] *= X
    if calculate_br_error: w['w_rr_err'] *= X 
        
    w['w_tr'] *= X
    w['w_tr_err'] *= X
    w['w_tr_BS'] = np.array(w['w_tr_BS']) * X[:,None]
    
    
    pre_b_r = Table.read(pre_b_r_filepath)
    
    def bias_corr(w_rr, w_rr_err=None):
        b_r_err = None
        b_r_SNR = None
    
        if b_r_triband_correction:
            tmp = w_rr * w['dz'] / w['w_m_ref_auto']
    
            tmp_minus = np.zeros_like(tmp)
            tmp_minus[1:] = (w['w_rr_minus'][1:] * w['dz'][:-1] * X[:-1]/ w['w_m_ref_auto'][:-1])
            tmp_minus[0] = tmp_minus[1]  # approximate 0th bin
    
            tmp_plus = np.zeros_like(tmp)
            tmp_plus[:-1] = (w['w_rr_plus'][:-1] * w['dz'][1:] * X[1:] / w['w_m_ref_auto'][1:])
            tmp_plus[-1] = tmp_plus[-2]
    
            b_r = np.sqrt(tmp + tmp_minus + tmp_plus)
    
        else:
            b_r = np.sqrt(w['dz'] * w_rr / w['w_m_ref_auto'])
    
        # --- error handling ---
        mask_bad = (w_rr <= 0) | (~np.isfinite(b_r))
        b_r = b_r.astype(float)
        b_r[mask_bad] = np.nan
    
        b_r_final = pre_b_r['b_r_w_highz_approx'].copy()
    
        if calculate_br_error==True:
            b_r_err = np.full_like(b_r, np.nan)
    
            valid = (~mask_bad) & (w_rr_err > 0)
            b_r_err[valid] = 0.5 * b_r[valid] * (w_rr_err[valid] / w_rr[valid])
            b_r_SNR = b_r / b_r_err
    
            use_on_the_fly = b_r_SNR >= 5
    
        else:
            frac_difference = np.abs(b_r - pre_b_r['b_r_w_highz_approx']) / pre_b_r['b_r_w_highz_approx']
            use_on_the_fly = frac_difference < 0.3
    
        b_r_final[use_on_the_fly] = b_r[use_on_the_fly]
        
        return b_r_final, b_r_err, b_r_SNR
    
    
    if calculate_br_error:
        b_r_final, b_r_err, b_r_SNR = bias_corr(w['w_rr'], w_rr_err=w['w_rr_err'])
        w['b_r_final'] = b_r_final
        w['b_r_err'] = b_r_err
        w['b_r_SNR'] = b_r_SNR
    else:
        b_r_final, b_r_err, b_r_SNR = bias_corr(w['w_rr'])
        w['b_r_final'] = b_r_final
    
    
    # ## Result & save
    
    def z_dist(w_tr, b_r_final, w_tr_err, w_tr_BS):
    
        dNdz_b = w_tr/b_r_final/w['w_m_cross']
        dNdz_b_err = w_tr_err/b_r_final/w['w_m_cross']
        dNdz_b_BS = (np.array(w_tr_BS).transpose()/b_r_final/w['w_m_cross']).transpose()
    
        return dNdz_b, dNdz_b_err, dNdz_b_BS
    
    
    dNdz_b, dNdz_b_err, dNdz_b_BS = z_dist(w['w_tr'], w['b_r_final'], w['w_tr_err'], w['w_tr_BS'])
    w['dNdz_b'] = dNdz_b
    w['dNdz_b_err'] = dNdz_b_err
    w['dNdz_b_BS'] = dNdz_b_BS
    
    
    # Correct for spatial resolution (see Appendix B of Chiang+26)
    
    def sigmoid(z):
        if config_par.zbin_num==160:
            L, U, k, z0 = 1.01, 1.56, 18.29, 1.05
        else:
            L, U, k, z0 = 1.01, 1.40, 20.67, 1.05
        return L + (U-L)/(1+np.exp(-k * (z-z0)))
        
    zless = np.argmin(np.abs(w['z']-3))
    zmore = np.argmin(np.abs(w['z']-.05))
    
    spatcorr = np.sqrt(np.array(sigmoid(w['z'][zmore:zless])))
    
    if rp_min<=0.5/h:
        w['dNdz_b'][zmore:zless] /= spatcorr
        w['dNdz_b_err'][zmore:zless] /= spatcorr
        w['dNdz_b_BS'][zmore:zless] /= spatcorr.reshape(len(spatcorr),1)
    
    
    if config_par.zbin_num==160:
        rebinned_w = utils.z_rebin_tomo_out(w, rebinned_z_edges, rebinned_z_ctrs)
    
    
    # ### Saving...
    
    # Change all 'dN' to 'dI' for intensity map case
    if config_par.test_type == 'intensity_map': 
        for col in w.colnames:
            if col[:2] == 'dN':
                col_new = 'dI' + col.split('dN')[1]
                w.rename_column(col, col_new)
    
        if config_par.zbin_num==160:
            for col in rebinned_w.colnames:
                if col[:2] == 'dN':
                    col_new = 'dI' + col.split('dN')[1]
                    rebinned_w.rename_column(col, col_new)
    

    import warnings
    from astropy.io.fits.verify import VerifyWarning
    warnings.filterwarnings('ignore', category=VerifyWarning)
    
    # Saving effective input parameters to 0th header
    hdr = fits.Header()

    hdr['REF_DATA'] = config_par.ref_sample
    hdr['TEST_TYPE'] = config_par.test_type
    hdr['TEST_DATA'] = config_par.test_data_file
    hdr['TEST_RANDOM'] = config_par.test_random_file
    if config_par.test_type=='intensity_map':
        hdr['TEST_MAP_ORDERING'] = config_par.test_map_ordering
    hdr['BEAM_FWHM_ARCMIN'] = config_par.beam_fwhm_arcmin
    hdr['FOOTPRINT_DEF'] = config_par.footprint_definition
    hdr['WEIGHT_MAP'] = ','.join(config_par.weight_path)
    hdr['TEMPLATE_CLEANING'] = config_par.template_cleaning
    hdr['FILTER_FWHM'] = config_par.filter_fwhm_deg if np.isfinite(config_par.filter_fwhm_deg) else 'inf'
    hdr['N_BOOTSTRAP'] = config_par.N_BS
    hdr['MASK_MAP'] = ', '.join(mask_map)
    
    
    primary_hdu = fits.PrimaryHDU(header=hdr)      
    hdu_w = fits.BinTableHDU(w)        
    if config_par.zbin_num==160:
        hdu_w_rebinned = fits.BinTableHDU(rebinned_w, name='REBINNED') 
        hdul = fits.HDUList([primary_hdu, hdu_w, hdu_w_rebinned])
    else:
        hdul = fits.HDUList([primary_hdu, hdu_w])
    
    # Write to file
    output_filename = config_par.output_filename

    fits_path = Path(output_filename + "_ddz.fits")
    ini_path = Path(output_filename + "_conf.ini")
    fig_path = Path(output_filename + "_ddz.pdf")
    
    # Create a copy of the configuration file using the same output filename
    fits_path.parent.mkdir(parents=True, exist_ok=True)
    
    hdul.writeto(fits_path, overwrite=True)
    shutil.copy2(f"{config_filename}", ini_path)

    
    # # Result visualisation
    
    default_font = rcParams['font.family']
    rcParams['font.family'] = 'monospace'
    plt.figure(figsize=(8, 3))
    
    if config_par.test_type=='source_catalog':
        key_col = 'dNdz_b'
        ylabel = r'b\ \mathrm{d}N/\mathrm{d}z\ \mathrm{[\#/deg^2]}'
    else:
        key_col = 'dIdz_b'
        ylabel = r'b\ \mathrm{d}I/\mathrm{d}z\ \mathrm{[map\ unit/deg^2]}'
    
    y1 = w[key_col]  
    y1err = w['{}_err'.format(key_col)]
    
    styles = {160: {'color': 'lightblue', 'ms': 3.5, 'capsize': 2}, 40: {'color': 'black', 'ms': 6, 'capsize': 3.5}}
    plt.errorbar(w['z']+1, y1, yerr=y1err, fmt='.', lw=0.8, **styles.get(len(w)))
    if config_par.zbin_num==160:
        plt.errorbar(rebinned_w['z']+1, rebinned_w[key_col], yerr=rebinned_w['{}_err'.format(key_col)], fmt='.', lw=0.8, zorder=100, **styles.get(len(rebinned_w)))
        
    plt.axhline(0, color='grey', ls='--', zorder=-10)
    
    plt.gca().set_xscale('log')
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    plt.gca().xaxis.set_minor_formatter(NullFormatter())
    zticks = np.array([0, 0.5, 1,1.5,  2, 2.5, 3, 3.5, 4.])
    zticks_str = ['0', '0.5', '1', '', '2', '', '3', '', '4']
    plt.gca().set_xticks(1.+zticks)
    plt.xlim([.97, 5.3])
    
    w_z03_pos = w[(w['z']<3.) & (w[key_col]>0.) ]
    w_z03_neg = w[(w['z']<3.) & (w[key_col]<0.) ]
    y_min, y_max = 0., 0.
    if len(w_z03_pos) > 0:
        w_z03_pos['SNR'] = w_z03_pos[key_col]/w_z03_pos[f'{key_col}_err']
        w_z03_pos.sort('SNR', reverse=True)
        w_z03_pos_top3 = w_z03_pos[:3]
        y_max = np.max(w_z03_pos_top3[key_col])
    if len(w_z03_neg) > 0:
        w_z03_neg['SNR'] = w_z03_neg[key_col]/w_z03_neg[f'{key_col}_err']
        w_z03_neg.sort('SNR')
        w_z03_neg_top3 = w_z03_neg[:3]
        y_min = np.min(w_z03_neg_top3[key_col])
    # enlarge the y range by 20% 
    y_mid = np.mean([y_min, y_max])
    y_range = y_max-y_min
    y_min = y_mid - (1.2*y_range)/2.
    y_max = y_mid + (1.2*y_range)/2.
    plt.ylim(y_min, y_max)
        
    plt.gca().set_xticklabels(zticks_str)
    plt.ylabel(r'${}$'.format(ylabel))
    plt.xlabel(r'$z$')
    
    plt.savefig(f"{fig_path}", dpi=DPI, bbox_inches='tight', pad_inches=0.02)

    logging.critical("Files saved:\n"
                    f"{fits_path}\n"
                    f"{ini_path}\n"
                    f"{fig_path}")
    
    end = time.time()
    logging.critical('\n' + '='*40 + '\nDURATION: {:.2f} minutes to finish.\n'.format((end-start)/60) + '='*40)
