# Tomographer

Tomographer is an integrated, end-to-end clustering-redshift estimator designed to make clustering-redshift analysis fast, robust, and accessible to non-experts, while retaining science-grade accuracy for astrophysical and cosmological applications.

It performs redshift inference on arbitrary extragalactic datasets using spatial cross-correlations with a fixed set of spectroscopic reference samples. It accepts either a <ins>source catalog</ins> (a list of sky coordinates) or a HEALPix <ins>intensity map</ins>, and returns the bias-weighted redshift distribution of the input, $b(z)\ {\rm d}N/{\rm d}z$ or $b(z)\ {\rm d}I/{\rm d}z$, in fine redshift bins over $0 < z \lesssim 4$.
 
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
│   │   ├── conf.ini             # General configuration file
│   │   ├── tomo.py              # Main code
│   │   └── tomo.ipynb           # Interactive notebook
│   ├── runtime_test/
│   │   ├── source_catalog/      # Example for running source catalog
│   │   └── intensity_map/       # Example for running intensity map
│   └── precalculated_data/      # Precomputed data for fast, pair-less correlations
├── pyproject.toml               # Build configuration and package metadata (pip install entry point)
└── README.md
```

The `precalculated_data` directory includes a lightweight configuration with 40 redshift bins for quick installation and testing. For science analyses, we recommend the extended 160-bin configuration. The additional data files are hosted on [Zenodo](https://zenodo.org/records/20155554) and downloaded automatically when needed, so no manual download is required.

---

## Quick Start

Users will mainly interact with the [configuration file](#configuration-file) `conf.ini`.

1. Edit [configuration file](#configuration-file) `conf.ini` as needed. <br>
   First-time users should start with `source_catalog/` or `intensity_map/`, which include complete example configurations and test data.

2. Run the pipeline from the directory containing the configuration file:
<table>
  <tr>
    <th width="50%">Command mode</th>
    <th width="50%">Interactive mode</th>
  </tr>
  <tr>
    <td valign="top">
      <pre><code> tomo </code></pre>
   If the configuration file is located in another directory or has a different name:
   <pre><code> tomo --conf /path/to/config </code></pre>
    </td>
    <td valign="top">
      For users who would like to explore the code in more detail or visualize intermediate results, we also provide the notebook:
<pre><code>tomo.ipynb</code></pre>
    </td>
  </tr>
</table>

---

## Configuration File

We here provide descriptions for the `Basic Settings` section in the configuration file `conf.ini`. First-time users are recommended to start with this section.

| Key | Required? | Value & Description *<small> (Click &#9654; to expand more) </small>* |
 :---: | :---: | :--- |
| `run_type` | Yes | Options: <ul><li>`quick`: 40 redshift bins <li>`fiducial`: 160 redshift bins
| `test_type` | Yes | Options<sup>[1]</sup>:  <details> <summary>`source_catalog` (Type C) </summary>  Table containing coordinates <br> (RA/Dec or Glon/Glat or l/b; case-insensitive column names). </details> <details> <summary> `intensity_map` (Type M) </summary> HEALPix map (NSIDE=2048). <br> Out-of-footprint or masked area set to NaN, otherwise, provide weight/footprint map  </details> |
| `test_data_file` | Yes | <details> <summary>`/path/to/test_data_file`</summary>  <ul><li>Data source catalog ([supported file formats](https://docs.astropy.org/en/latest/io/unified_table.html#built-in-table-readers-writers): FITS, CSV, Parquet, HDF5, etc.) </li> <li> Data intensity field (HEALPix 2048)<sup> [2] </sup></li></ul></details> |
| `test_map_ordering` | Yes <br> (Type M only) | Options: <ul><li>`NEST` <li>`RING`|
| `beam_fwhm_arcmin` | Yes <br> (Type M only) | <details> <summary> &geq; 1.24  </summary> For source catalogs, the smallest $r_\mathrm{p, min} = 0.5/h\mathrm{Mpc}$ is used. <br> For intensity maps, the beam FWHM corresponds to $$r_{\mathrm{p,min}}/h\mathrm{Mpc} = {\lfloor \mathrm{FWHM}/1.5 \rfloor}/2$$ with constraints within the range $[0.5, 2.5]/h \mathrm{Mpc}$. </details> |
|`footprint_definition` | Yes | Options: <details> <summary> `user_defined` </summary> **RECOMMENDED** if the selection function is well understood. <br> A test random file and/or spatial weight map must be provided. <br> If a test random file is provided, this mode is selected automatically. </details> <details>  <summary> `full_sky` </summary> Apply no selection; use all-sky pixels. <br> Any provided spatial weight map is ignored. <br> A rectangular cut can be applied in Advanced Settings after selecting this mode. </details> <details> <summary> `auto_detection` </summary> Use with caution; provide only a coarse footprint estimate. <br> Run [auto_footprint](https://github.com/yuvoonng/tomographer/blob/442954f4ee2cc4ade43db0d3c9e305c9bb9bbdb9/tomographer/runtime_script/tomo_utils.py#L611) on the test data file. <br> Any provided spatial weight map is ignored.</details>|
| `test_random_file` | No | <details> <summary>`/path/to/test_random_file`</summary>  To capture the selection function <br> Random source catalog, or random intensity field (same format as test data) </details> |
| `spatial_weight_map_file_1`<sup> [3] </sup> | No | <details> <summary> `/path/to/spatial_weight_map_file` <sup> [2] </sup> </summary> A simple footprint map (0/1 for masked/unmasked areas) or continuous spatial weights. </details>  |
| `spatial_weight_map_ordering_1`<sup> [3] </sup> | No | Options: <ul><li>`NEST` <li>`RING` |
| `per_source_weight_colname` | No (Type C only) | <details> <summary> `column_name` </summary> The column name specifying additional weights for each input source in data and random file (fits or csv) </details> |
| `output_prefix` | Yes | <details> <summary>`/path/to/output_prefix` </summary> Prefix of output filenames. <ul><li> `<prefix>_ddz.fits`: result table </li> <li> `<prefix>_ddz.pdf`: plot </li> <li> `<prefix>_conf.ini`: copy of configuration </li></ul>Can also accept a custom path. 

**[1]** A discrete count map (e.g., a source catalog read onto HEALPix) is also treated as an intensity map. <br>
**[2]** The map files are read with [`hp.read_map()`](https://healpy.readthedocs.io/en/latest/generated/healpy.fitsfunc.read_map.html), so please follow the expected format accordingly. <br>
**[3]** Replace `1` with other integers to specify additional spatial weight map files.

<br>

For advanced configurations, users can modify the settings in the `Advanced Settings` section. This includes:
<details> <summary> <strong> Template cleaning </strong>: Regress out a list of foreground templates for intensity maps </summary>
<ul>
<li> HI </li>
<li> CSFD </li>
</ul>
</details>
<details> <summary> <strong> High-pass filtering</strong>: Suppress large-scale Galactic foregrounds and survey systematics </summary>
<ul>
 <li> 2&deg; </li>
 <li> 4&deg; </li>
 <li> 8&deg; </li>
 <li> 16&deg; </li>
 <li> inf (no filtering) </li>
</ul>
</details>
<details> <summary> <strong> Masking </strong>: Rectangular cut and built-in masks </summary>
<ul> 
<li> <code>CSFD_cosmology_area</code>: Dust mask based on the corrected SFD map of <a href="https://doi.org/10.3847/1538-4357/acf4a1">Chiang (2023)</a></li>
<li> <code>fsky_[25,50,75]_low_dust</code>: Galactic-latitude cuts (|b| > 25°, 50°, 75°) </li>
<li> <code>globular_clusters_veto</code>: Globular clusters </li>
<li> <code>LMC_SMC_veto</code>: Magellanic Clouds </li>
<li> <code>nearby_galaxy_clusters_veto</code>: Nearby galaxy clusters </li>
<li> <code>Planck_point_source_veto</code>: Planck point sources </li>
<li> <code>SDSS_veto</code>: SDSS footprints </li>
</ul>
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
