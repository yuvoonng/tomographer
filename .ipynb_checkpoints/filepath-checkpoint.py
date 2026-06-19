from pathlib import Path

# zbin_tag = '160log1pzbins'
# ref_sample = 'sdss_lss_plus'
# scale_tag = '0.5-10Mpch'
# gamma = -0.8
# max_theta_max = 15.0

def get_filepath(fileloc, name, **kwargs):
    fileloc = '../Data'
    home_dir = Path("{}".format(fileloc))

    ref_dir = home_dir / 'references'
    ref_mask_dir = ref_dir / 'footprint'
    act_dir = Path("/home/ykchiang/workdir/Projects/Tomographer_2.1_Online/code/data") / 'activation_maps'
    mask_dir = home_dir / 'masks'
    temp_dir = home_dir / 'templates'
    wm_dir = home_dir / 'cosmology' / 'planck_2018'
    BS_dir = home_dir / 'bootstrapping'
    
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
    
    path_template = filepath.get(name)
    if path_template is None:
        raise ValueError(f"Invalid mask name: {name}")

    # Convert to string and format with kwargs if placeholders exist
    path_str = str(path_template).format(**kwargs)
    path = Path(path_str)

    if not path.exists():
        raise ValueError(f'The file does not exist: {path}')

    return path
