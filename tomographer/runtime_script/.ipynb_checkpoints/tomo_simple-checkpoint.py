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
import shutil
import psutil
import time
import os
import gc

from . import tomo_utils as utils
    
def main():

    start = time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--conf', type=str, default='conf.ini', help='Path to the config file')
    args = parser.parse_args()
    
    logging.basicConfig(
        level=logging.INFO,           
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    
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
        raise FileNotFoundError(f"Config file {config_filename} does not exist.")
    logging.info(f"Loading configuration file {config_filename}...")

    # ==== Initiate Config ==== 
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(f'{config_filename}')
    config_par = utils.ConfigChecker(config).check_all() # Check the compatibility of the configuration file
    