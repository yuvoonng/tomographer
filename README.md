# Tomographer

Tomographer is an an integrated clustering-redshift estimator designed to be fast, robust, and usable end-to-end by non-experts, while
retaining science-grade accuracy for survey collaborations. \
It accepts either a <ins>source catalog</ins> (a list of sky coordinates) or a HEALPix <ins>intensity map</ins>, and returns the bias-weighted redshift distribution of the input, $b(z)\ {\rm d}N/{\rm d}z$ or $b(z)\ {\rm d}I/{\rm d}z$, in fine redshift bins over $0 < z \lesssim 4$.\
For details of the methodology, see Chiang et al. (2026).

---

## Installation

1. Fetch the repository from [Github](https://github.com/yuvoonng/tomographer).
```bash
git clone https://github.com/yuvoonng/tomographer.git
```   

2. Install the package
```bash
cd tomographer
pip install .
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
│   └── precalculated_data/      # Downloaded precalculated files
├── pyproject.toml               # Build configuration and package metadata (pip install entry point)
└── README.md
```

The `precalculated_data` directory includes the default configuration with 40 tomographic bins.
We also provide an extended version with 160 bins. Using the 160-bin configuration requires additional data downloads from [Zenodo](https://zenodo.org/records/20155554). These downloads are triggered automatically by the code based on the selected configuration, so users do not need to download them manually.

---

## Quick Start
    
<!-- <table style="table-layout: fixed; width: 100%; border-collapse: collapse;">
<tr>
<td width="50%" style="vertical-align: top; font-size: 14px; line-height: 1.4; padding: 0 10px 0 0; word-break: break-all; max-width: 50%;">
<strong>First-time users</strong>
<hr style="margin: 5px 0;">
1. Go to one of the runtime directories within the downloaded package:<br>
<code>runtime_script/</code><br>
<code>runtime_test/source_catalog/</code><br>
<code>runtime_test/intensity_map/</code><br><br>
2. Modify the <a href="#configuration-file">configuration file</a> <code>conf.ini</code>. It is recommended to start from one of the example directories (<code>source_catalog/</code> or <code>intensity_map/</code>), which include complete example configurations and test data.<br><br>
3. Run the pipeline from the directory to initialize the precalculated data path:
<pre style="max-width: 100%; white-space: pre-wrap; word-break: break-all;"><code style="white-space: pre-wrap; word-break: break-all;">tomo</code></pre>
</td>
<td width="50%" style="vertical-align: top; font-size: 14px; line-height: 1.4; padding: 0 0 0 10px; word-break: break-all; max-width: 50%;">
<strong>Advanced users</strong>
<hr style="margin: 5px 0;">
Users will mainly interact with the <a href="#configuration-file">configuration file</a> <code>conf.ini</code>.<br><br>
1. Modify the <a href="#configuration-file">configuration file</a> <code>conf.ini</code> according to your needs.<br><br>
2. Run the pipeline from the directory containing the configuration file:
<pre style="max-width: 100%; white-space: pre-wrap; word-break: break-all;"><code style="white-space: pre-wrap; word-break: break-all;">tomo</code></pre>
If the configuration file is located in another directory or has a different name:
<pre style="max-width: 100%; white-space: pre-wrap; word-break: break-all;"><code style="white-space: pre-wrap; word-break: break-all;">tomo --conf /path/to/your/config/file</code></pre>
</td>
</tr>
</table> -->

<table>
  <tr>
    <th width="50%">First-time users</th>
    <th width="50%">Advanced users</th>
  </tr>
  <tr>
    <td valign="top">
      <ol>
      <li> Go to one of the runtime directories within the downloaded package: <br>
        <code>runtime_script/</code> <br>
        <code>runtime_test/source_catalog/</code> <br>
        <code>runtime_test/intensity_map/</code> </li> <br>
      <li> Modify the <a href="#configuration-file">configuration file</a> <code>conf.ini</code>. It is recommended to start from one of the example directories (<code>source_catalog/</code> or <code>intensity_map/</code>), which include complete example configurations and test data. </li>  <br>
      <li> Run the pipeline from the directory to initialize the precalculated data path: <br>
      <pre><code>tomo</code></pre> </li> </ol>
    </td>
    <td valign="top">
      Users will mainly interact with the <a href="#configuration-file">configuration file</a> <code>conf.ini</code>. <br><br>
    <ol>
      <li> Modify the <a href="#configuration-file">configuration file</a> <code>conf.ini</code> according to your needs. </li> <br>
      <li> Run the pipeline from the directory containing the configuration file: <br>
      <pre><code>tomo</code></pre> 
      If the configuration file is located in another directory or has a different name: <br>
      <pre><code>tomo --conf /path/to/config</code></pre> </li> </ol>
    </td>
  </tr>
</table>

At the *first* run, the code searches for the `precalculated_data` directory within the downloaded package and stores its path in 
`~/.precal_data_path.toml`. If the `precalculated_data` directory is later moved to another location, please update the path file accordingly. Otherwise, the code will not be able to locate the directory.

For users who would like to explore the code in more detail or visualize intermediate results, we also provide the notebook:
```text
tomo.ipynb
```

---

## Configuration File

We here provide descriptions for the `Basic Settings` section in the configuration file `conf.ini`. First-time users are recommended to start with this section.

| Key | Required? | Value & Description *<small> (Click &#9654; to expand more) </small>* |
 :---: | :---: | :--- |
| `run_type` | Yes | Options: <ul><li>`quick`: 40 redshift bins <li>`fiducial`: 160 redshift bins
| `test_type` | Yes | Options<sup>[1]</sup>:  <details> <summary>`source_catalog` (Type C) </summary>  Table containing coordinates <br> (RA/Dec or Glon/Glat or l/b; case-insensitive column names). </details> <details> <summary> `intensity_map` (Type M) </summary> HEALPix map (NSIDE=2048). <br> Out-of-footprint or masked area set to NaN, otherwise, provide weight/footprint map  </details> |
| `test_data_file` | Yes | <details> <summary>`/path/to/test_data_file`</summary>  <ul><li>Data source catalog ([supported file formats](https://docs.astropy.org/en/latest/io/unified_table.html#built-in-table-readers-writers): FITS, CSV, Parquet, HDF5, etc.) </li> <li> Data intensity field (HEALPix 2048)<sup> [2] </sup></li></ul></details> |
| `test_map_ordering` | Yes <br> (Type M only) | Options: <ul><li>`NEST` <li>`RING`|
| `beam_fwhm_arcmin` | Yes <br> (Type M only) | <details> <summary> &geq; 1.24  </summary> For source catalogs, the smallest $r_\mathrm{p, min} = 0.5/h\mathrm{Mpc}$ is used. <br> For intensity maps, the beam FWHM correponds to $$r_{\mathrm{p,min}}/h\mathrm{Mpc} = {\lfloor \mathrm{FWHM}/1.5 \rfloor}/2$$ </details> |
|`footprint_definition` | Yes | Options: <details> <summary> `user_defined` </summary> **RECOMMENDED** if the selection function is well understood. <br> If a test random file is provided, this mode is selected automatically. <br> A spatial weight map may also be provided (not mutually exclusive), but is required if no random file is given. </details> <details>  <summary> `full_sky` </summary> Apply no selection; use all-sky pixels. <br> Any provided spatial weight map is ignored. <br> A rectangular cut can be applied in Advanced Settings after selecting this mode. </details> <details> <summary> `auto_detection` </summary> Use with caution; provide only a coarse footprint estimate. <br> Run auto_footprint_detection on the test data file. <br> Any provided spatial weight map is ignored.</details>|
| `test_random_file` | No | <details> <summary>`/path/to/test_random_file`</summary>  To capture the selection function <br> Random source catalog, or random intensity field (same format as test data) </details> |
| `spatial_weight_map_file_1`<sup> [3] </sup> | No | <details> <summary> `/path/to/spatial_weight_map_file` <sup> [2] </sup> </summary> A simple footprint map (0/1 for masked/unmasked areas) or continuous spatial weights. </details>  |
| `spatial_weight_map_ordering_1`<sup> [3] </sup> | No | Options: <ul><li>`NEST` <li>`RING` |
| `per_source_weight_colname` | No (Type C only) | <details> <summary> `column_name` </summary> The column name specifying additional weights for each input source in data and random file (fits or csv) </details> |
| `filename` | Yes | <details> <summary>`/path/to/output_files` </summary> Your output filename. Can also accept a custom path. </details>

**[1]** A discrete count map (e.g., a source catalog read onto HEALPix) is also treated as an intensity map. <br>
**[2]** The map files are read with [`hp.read_map()`](https://healpy.readthedocs.io/en/latest/generated/healpy.fitsfunc.read_map.html), so please follow the expected format accordingly. <br>
**[3]** Replace `1` with other integers to specify additional spatial weight map files.

<br>

For advanced configurations, users can modify the settings in the `Advanced Settings` section. This includes:
- **Template cleaning**: Regress out a list of foreground templates for intensity maps
- **High-pass filtering**: Suppress large-scale Galactic foregrounds and survey systematics
- **Masking**: Rectangular cut and built-in masks


---

## Output

Tomographer produces:

- Measurement file: `OUTPUT_FILENAME.fits`

- Saved configuration file: `OUTPUT_FILENAME.ini`

- Result plot: `OUTPUT_FILENAME.pdf`

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

Note: If the precalculated data path is not initialized first, run the above code from a runtime directory for initialization.
</details>

<details> <summary> <strong> Why do I get <code>FileNotFoundError: precalculated_data not found.</code>? </strong> </summary>
Please check whether the precalculated data path has been initialized by verifying the existence of the file <code>~/.precal_data_path.toml</code>. 

If the file does not exist, refer to <strong> First-time users</strong> in [Quick Start](#quick-start).
<br> If it exists, ensure that the path stored inside it correctly points to the location of the precalculated data files.
</details>

<details> <summary> <strong> What built-in masks does Tomographer provide? </strong> </summary>
We provide the following built-in masks that users can optionally apply:
<ul><li> Galactic-latitude cuts (|b| > 25°, 50°, 75°) </li>
<li> Magellanic Clouds </li>
<li> Globular clusters </li>
<li> Nearby galaxy clusters </li>
<li> Planck point sources </li>
<li> SDSS footprints </li>
<li> Dust mask based on the corrected SFD map of <a href="https://doi.org/10.3847/1538-4357/acf4a1">Chiang (2023)</a> </li>
</ul>
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
