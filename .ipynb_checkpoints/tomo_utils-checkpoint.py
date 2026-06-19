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

cmap = cm.RdYlBu_r

# ==== Version ====
# dev_v06
# 01/25/2024

# ==== Data Handeling ====

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
