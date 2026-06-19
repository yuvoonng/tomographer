import numpy as np
import numexpr as ne
import healpy as hp
from scipy import interpolate
import pandas as pd
from sklearn.linear_model import LinearRegression

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

from multiprocessing import Pool
from importlib import resources
from pathlib import Path
import configparser
import argparse
import logging
import shutil
import psutil
import time
import os
import gc
import asyncio

from . import tomo_utils as utils

def main():

    logging.basicConfig(
    level=logging.INFO,           
    format="%(asctime)s - %(levelname)s - %(message)s"
    )
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--conf', type=str, default='conf.ini', help='Path to the config file')
    args = parser.parse_args()

    filepath_handler = utils.FilePathHandler()
    try:
        data_path = resources.files("precalculated_data")
        filepath_handler.home_dir = data_path
    except:
        filepath_handler.load_or_init_precaldata()

    config_filename = args.conf
    if not Path(config_filename).exists():
        raise FileNotFoundError(f"Config file {config_filename} does not exist.")
    
    logging.info(f"Loading configuration file {config_filename}...")
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(f'{config_filename}')
    config_par = utils.ConfigChecker(config).check_all()
    get_filepath = filepath_handler.get_filepath
    
    # ==== Scale ==== 
    # rp_max is fixed to 10 Mpc/h 
    # rp_min can be [0.5, 1, 1.5, 2, 2.5] Mpc/h. Needs to be larger than the beam to avoid a strong 1-halo term signal
    
    zbin_tag = '{}log1pzbins'.format(config_par.zbin_num) 
    # ref_sample_filepath = get_filepath('ref_sample', zbin_tag=zbin_tag, ref_sample=config_par.ref_sample)

    rp_min = (config_par.beam_fwhm_arcmin//1.5)/2. 
    rp_min = min(max(rp_min, 0.5), 2.5) # (min, max) = (0.5, 2.5)
    
    if np.isclose(rp_min%1, 0):
        scale_tag = str(int(rp_min))+'-10Mpch'
    else:
        scale_tag = str(rp_min)+'-10Mpch'
    print('(rp_min, rp_max) = {}'.format(scale_tag))
    
    max_theta_max = 15. # Deg, maximum angular scale used, which overwrites rp_max at the lowest z
    gamma = -0.8 # Scale weighting power index
    
    # act_data_filepath = get_filepath('act_data', zbin_tag=zbin_tag, ref_sample=config_par.ref_sample, scale_tag=scale_tag, gamma=gamma)
    # act_rand_filepath = get_filepath('act_rand', zbin_tag=zbin_tag, ref_sample=config_par.ref_sample, scale_tag=scale_tag, gamma=gamma)

    matter = get_filepath('wm', scale_tag='0.5-10Mpch', gamma=-0.8, max_theta_max=15.0, zbin_tag='160log1pzbins')

    # ref_data_z_maps = pd.read_parquet(ref_sample_filepath, columns=['0'])  
    if len(filepath_handler.need_download)!=0:
        asyncio.run(filepath_handler.async_download(filepath_handler.need_download))
