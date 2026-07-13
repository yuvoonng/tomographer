<h1><span style="font-variant: small-caps;">Tomographer-2</span></h1>

Tomographer is an integrated, end-to-end clustering-redshift estimator designed to make clustering-redshift analysis fast, robust, and accessible to non-experts, while retaining science-grade accuracy for astrophysical and cosmological applications.

It performs redshift inference on arbitrary extragalactic datasets using spatial cross-correlations with a fixed set of spectroscopic reference samples. It accepts either a <ins>source catalog</ins> (a list of sky coordinates) or a HEALPix <ins>intensity map</ins>, and returns the bias-weighted redshift distribution of the input, $b(z)\ {\rm d}N/{\rm d}z$ or $b(z)\ {\rm d}I/{\rm d}z$, respectively, in fine redshift bins over $0 < z < 4.2$.
 
For details, see Chiang et al. (2026).

---

## Installation
0. Git Large File Storage (LFS): 
Skip if installed. Otherwise, install [Git LFS](https://git-lfs.com/) to enable cloning large files. Example using `conda`:
```bash
conda install -c conda-forge git-lfs
git lfs install
```  

1. Clone [Tomographer](https://github.com/yuvoonng/tomographer) from Github (~3 GB).
```bash
git clone https://github.com/yuvoonng/tomographer.git
```   

2. Install the package.
```bash
cd tomographer
pip install .
```

3. Initialize the code.   
```bash
tomo init
```

---

## Project Structure

```text
tomographer/
├── tomographer/
│   ├── runtime_script/
│   │   ├── conf.ini             # Configuration file as user interface
│   │   ├── tomo.py              # Main code script
│   │   └── tomo.ipynb           # Main code equivalent Jupyter notebook
│   ├── runtime_demo/
│   │   ├── source_catalog/      # Example run for source catalog
│   │   └── intensity_map/       # Example run for intensity map
│   └── precalculated_data/      # Precomputed data for fast, pair-less correlations
├── pyproject.toml               # Build configuration and package metadata for pip install
└── README.md
```

The `precalculated_data` directory initially includes a lightweight 40-redshift-bin configuration for quick installation and testing. For science analyses, we recommend the full 160-bin configuration. The additional files are hosted on [Zenodo](https://zenodo.org/records/20155554) and are downloaded automatically when needed; no manual download is required.

---

## Quick Start

Users will mainly interact with the [configuration file](#configuration-file) `conf.ini`.

1. Edit [configuration file](#configuration-file) `conf.ini` as needed. <br>
   First-time users can start in `runtime_demo`, with ready-to-run examples in `source_catalog/` and `intensity_map/`.

2. Run tomographer from the directory containing `conf.ini` using one of the two equivalent modes:
<table>
  <tr>
    <th width="50%">Command Line Mode</th>
    <th width="50%">Jupyter Notebook Mode</th>
  </tr>
  <tr>
    <td width="50%" valign="top">
      <div>To run a full end-to-end analysis, simply execute:</div>
      <br>
      <div><code>tomo</code></div>
      <br>
      <div>If the configuration is in another directory or has a different file name, run <code>tomo --conf /path/to/config</code></div>
    </td>
    <td width="50%" valign="top">
      <div>To run interactively with step-by-step visualization, use Jupyter Notebook:</div>
      <br>
      <div><code>tomo.ipynb</code></div>
      <br>
      <div>Modify <code>config_filename</code> therein if using a different configuration </div>
    </td>
  </tr>
</table>

---

## Configuration File

We here provide descriptions for the `Basic Settings` section in the configuration file `conf.ini`. First-time users are recommended to start with this section.

| Key | Required? | Value & Description *<small> (Click &#9654; to expand more) </small>* |
 :---: | :---: | :--- |
| `run_type` | Yes | Options: <ul><li>`quick`: 40 redshift bins <li>`fiducial`: 160 redshift bins
| `test_type` | Yes | Options<sup>[1]</sup>:  <details> <summary>`source_catalog` (Type C) </summary>  Table containing coordinates <br> (RA/Dec or Glon/Glat or l/b; case-insensitive column names). </details> <details> <summary> `intensity_map` (Type M) </summary> HEALPix map (NSIDE=2048). <br> Out-of-footprint or masked area set to NaN; otherwise, provide a weight/footprint map.  </details> |
| `test_data_file` | Yes | <details> <summary>`/path/to/test_data_file`</summary>  <ul><li>Data source catalog ([supported file formats](https://docs.astropy.org/en/latest/io/unified_table.html#built-in-table-readers-writers): FITS, CSV, Parquet, HDF5, etc.) </li> <li> Data intensity field (HEALPix 2048)<sup> [2] </sup></li></ul></details> |
| `test_map_ordering` | Yes <br> (Type M only) | Options: <ul><li>`NEST` <li>`RING`|
| `beam_fwhm_arcmin` | Yes <br> (Type M only) | <details> <summary> Beam FWHM of the test intensity map </summary> For undersampled, unsmoothed maps, use 1.24 arcmin as a Gaussian approximation to the NSIDE=2048 pixel window. <br>To avoid 1-halo clustering signal, `rp_min` is derived from the beam FWHM, quantized to a 0.5 Mpc/h grid over [0.5, 2.5]. <br>`rp_max` is fixed at 10 Mpc/h. </details> |
|`footprint_definition`<sup>[3]</sup> | Yes | Options: <details> <summary> `user_defined` </summary> **Recommended** for science runs. <br> Provide a random catalog (Type C) and/or a spatial weight map (Type C/M). <br> Selected automatically if a test random file is provided. </details> <details>  <summary> `full_sky` </summary> Treat the full sky (Type C) or all non-NaN pixels (Type M) as valid.<br>Useful for truly full-sky data, pre-masked intensity maps, or when using built-in masks in `Advanced Settings`.</details> <details> <summary> `auto_detection` </summary> Use with caution.<br>Estimate a coarse footprint using [auto_footprint](https://github.com/yuvoonng/tomographer/blob/442954f4ee2cc4ade43db0d3c9e305c9bb9bbdb9/tomographer/runtime_script/tomo_utils.py#L611) from the non-zero areas in the test data; susceptible to boundary effects.</details>|
| `test_random_file` | No | <details> <summary>`/path/to/test_random_file`</summary>  Random source catalog or random HEALPix intensity map (same format as test data)<br>Used only when <code>footprint_definition = user_defined</code>. </details> |
| `spatial_weight_map_file` | No | <details> <summary> `/path/to/spatial_weight_map_file` <sup> [2] </sup> </summary> A simple footprint map (0/1 for masked/unmasked areas) or continuous spatial weights.<br>Used only when <code>footprint_definition = user_defined</code>. </details>  |
| `spatial_weight_map_ordering` | No | Options: <ul><li>`NEST` <li>`RING` |
| `per_source_weight_colname` | No (Type C only) | <details> <summary> `column_name` </summary> The column name specifying additional weights for each input source in data and random file. </details> |
| `output_prefix` | Yes | <details> <summary>`/path/to/output_prefix` </summary> Prefix of output filenames. <ul><li> `<prefix>_ddz.fits`: result table </li> <li> `<prefix>_ddz.pdf`: plot </li> <li> `<prefix>_conf.ini`: copy of configuration </li></ul>Can also accept a custom path. 

**[1]** A discrete count map (e.g., a source catalog read onto HEALPix) is also treated as an intensity map. <br>
**[2]** The map files are read with [`hp.read_map()`](https://healpy.readthedocs.io/en/latest/generated/healpy.fitsfunc.read_map.html), so please follow the expected format accordingly. <br>
**[3]** Footprint, veto masks, and spatial weighting are all represented as full-sky weight maps and multiplied together into a single selection function. Additional masks, rectangular cuts, and weight maps can be specified in `Advanced Settings`.

<br>

For advanced configurations, users can modify the settings in the `Advanced Settings` section. This includes:
<details> <summary> <strong> Template cleaning </strong>: Linearly regress out a foreground template from the test field. </summary>
<ul>
<li> <code>HI</code>: Galactic HI from HI4PI (<a href="https://doi.org/10.1051/0004-6361/201629178">HI4PI Collaboration 2016</a>). </li>
<li> <code>CSFD</code>: Large-scale-structure-free Galactic dust map (<a href="https://doi.org/10.3847/1538-4357/acf4a1">Chiang 2023</a>).</li>
</ul>
</details>
<details> <summary> <strong> High-pass filtering</strong>: High-pass filter the test field to suppress large-scale foregrounds and systematics. </summary>
<ul>
 <li> 2&deg; </li>
 <li> 4&deg; </li>
 <li> 8&deg; </li>
 <li> 16&deg; </li>
 <li> inf (no filtering) </li>
</ul>
</details>
<details> <summary> <strong> Masking </strong>: Rectangular cut and built-in masks. </summary>
<ul> 
<li> <code>CSFD_cosmology_area</code>: UV-optical-IR extragalactic sky outside high dust and WISE veto</li>
<li> <code>fsky_[25,50,75]_low_dust</code>: Lowest [25,50,75]% dust sky by CSFD E(B-V) </li>
<li> <code>globular_clusters_veto</code>: Globular clusters </li>
<li> <code>LMC_SMC_veto</code>: Magellanic Clouds </li>
<li> <code>nearby_galaxy_clusters_veto</code>: Nearby galaxy clusters </li>
<li> <code>Planck_point_source_veto</code>: Planck point sources </li>
<li> <code>SDSS_veto</code>: Veto bad imaging stripes/scans in the original SDSS imaging </li>
</ul>
</details>
<details> <summary> <strong> Additional spatial weight maps </strong>:  Same format as <code>spatial_weight_map_file</code>. All masks and weight maps are multiplied together. </summary>
To add more maps, replace <code>int</code> in <code>additional_weight_map_file_[int]</code> with different integer values.
</details>


---

## Output

Tomographer produces:

- Measurement file: `OUTPUT_PREFIX_ddz.fits`

- Saved configuration file: `OUTPUT_PREFIX_conf.ini`

- Result plot: `OUTPUT_PREFIX_ddz.pdf`

---

## FAQ

<details> <summary> <strong> Can I download all precalculated files from Zenodo at once? </strong> </summary>

Yes. To download all precalculated files from Zenodo in one step and automatically organize them into the correct directory structure:

```python
from tomographer.runtime_script import tomo_utils as utils
filepath_handler = utils.FilePathHandler()
filepath_handler.load_or_init_precaldata()
filepath_handler.download_all()
```
</details>

<details> <summary> <strong> Why do I get <code>FileNotFoundError: precalculated_data not found.</code>? </strong> </summary>

 During initialization, the code searches for the <code>precalculated_data</code> directory within the downloaded package and stores its path in 
<code>~/.tomographer_path.toml</code>

First, verify that the initialization completed successfully by checking whether this file exists. If it does not, run the following command from within the downloaded package:
 ```bash
 tomo init
 ```

If the file exists, ensure that the stored path correctly points to the location of the <code>precalculated_data</code> directory. <br>
If the <code>precalculated_data</code> directory is later moved to a different location, either update the path stored in <code>~/.tomographer_path.toml</code> or rerun the initialization. Otherwise, Tomographer will not be able to locate the precalculated data files.

</details>

<details>
<summary><strong> Can I rebin my measurements into fewer bins to improve the signal-to-noise ratio?</strong></summary>

Yes. To rebin a measurement into a smaller number of bins, run:
```bash
tomo rebin YOUR_MEASUREMENT_FILE.fits BIN_NUMBER
```

e.g. To rebin the measurements in <code>out.fits</code> into 20 bins:
```bash
tomo rebin out.fits 20
```

The rebinned measurements will be saved as <code>ORIGINAL_FILENAME_rebin{BIN_NUMBER}.fits</code> (e.g. <code>out_rebin20.fits</code>)

</details>

<details>
<summary><strong>What if I am interested in this framework but need a more advanced analysis pipeline?</strong></summary>
We welcome collaborations based on the Tomographer framework. If your scientific application requires a more customized or sophisticated analysis workflow, please feel free to reach out using the contact information below.
</details>

<details> <summary> <strong> What should I do if I encounter an issue that is not listed here? </strong> </summary>
Feel free to join the discussion group and ask your question there. Your question may help other users facing the same issue.
</details>


---

## Contact

Yi-Kuan Chiang: ykchiang (at) asiaa.sinica.edu.tw\
Yu Voon Ng: yvng (at) asiaa.sinica.edu.tw\
Discussion group: 
