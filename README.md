# ch4ne: TROPOMI methane screening with false-positive control over onshore oil and gas pipelines in Northeast Brazil

Code (version 4.1.0) for the article:

> Cardoso, J. F., et al. *Detectability of pipeline-scale methane enhancements from TROPOMI under false-positive control: a screening protocol for onshore oil and gas networks in Northeast Brazil.* Manuscript submitted for publication (2026).

The pipeline screens six years (2020–2025) of TROPOMI/Sentinel-5P column methane (XCH₄) over the nine states of Northeast Brazil. It removes a monthly regional background and a local moving-window background, applies minimum-sampling rules, flags robust z-scores and estimates false discoveries with a mirror (knockoff+) test. Cells near the onshore pipelines of the National Agency of Petroleum, Natural Gas and Biofuels (ANP) are then compared with distant cells using spatial block-bootstrap, regression, matching and coastal stratification.

**Authors:** Jean Firmino Cardoso, Hugo Emanoel Bezerra de Lima, Bruno Felipe de Carvalho, Vitória Regina da Silva Bulhões, Rebeca Carolina de Albuquerque Oliveira Lancha, Daiane Francisca do Nascimento Silva, Hugo Barcelar de Araújo Filho, José Roberto da Silva, Ivan García-Fornaris, Alexandre Carlos Araújo de Santana, Yaicel Ge Proenza, Abel Gámez Rodríguez, Raquel Milani, Yelina González Pérez, Severino Martins dos Santos Neto, Allan Almeida Albuquerque, Clériston Moura Vieira Júnior, Emmanuel Damilano Dutra, Marcelo Metri Corrêa, Marcus Metri Corrêa, João Paulo Cerquinho Cajueiro, Manoel Henriques de Sá Campos Filho, José Miguel Reichert, Daniel Milian Pérez, Antonio Celso Dantas Antonino. ORCIDs and affiliations are in [CITATION.cff](CITATION.cff).

Contact: Jean Firmino Cardoso (jean.firmino@ufpe.br), Departamento de Energia Nuclear, Universidade Federal de Pernambuco (UFPE), Recife, Brazil.

A Portuguese version of this file is in [README.pt-BR.md](README.pt-BR.md). The notes in `docs/` are in Portuguese.

> **Data policy:** this repository neither contains nor generates synthetic or invented data. The analysis runs only on rasters downloaded from Google Earth Engine, with a provenance record (`PROVENANCE.json`, with SHA-256 checksums) that the code checks before use. The tests use real inputs only: the ANP pipeline layers, the state boundaries and the tables published with the earlier version (v3).

## Supersedes v1.0.0 (v3 notebooks)

Version 4 replaces the v2/v3 notebooks archived in Zenodo as version 1.0.0 (https://doi.org/10.5281/zenodo.18717652). An audit found errors in that version, documented in [docs/AUDITORIA_v3.md](docs/AUDITORIA_v3.md) (Portuguese). **Do not use v1.0.0.**

| v3 | v4 |
|---|---|
| State names rebuilt by index (6 of 9 swapped) | State always taken from the polygon attribute, with an independent point-in-polygon check |
| Cluster–pipeline distances always 0 (read from the wrong object) | Real geodesic distances (SIRGAS 2000/UTM zone per cell; azimuthal equidistant per cluster) |
| Six-year mean without removing the trend | Per-orbit anomaly = XCH₄ minus the regional median of the month |
| Z against a single regional mean and SD | Local background (moving median, about 220 km) and robust z (median/MAD) |
| No minimum number of observations | At least 50 retrievals and 24 months with data per cell |
| `reproject` (sampling) to 7 km | `reduceResolution` (area mean) on a fixed 0.0625° grid |
| "Morphological cleaning" that did not erode, at 10 km | Connected components on the native grid and a false-discovery q-value (mirror test, knockoff+) |
| Nested zones; spatial SD read as temporal | Exclusive rings; real temporal variance; frequency of positive and negative deviations; matched controls (coast, albedo, water, elevation); spatial block bootstrap; regressions with block-robust errors |
| Pipelines with a 1° buffer (including Espírito Santo and offshore) | Classified inventory (onshore Northeast / offshore / other state); distances to the on-land parts only |
| Figures from thumbnails with embedded titles | Figures from the GeoTIFFs, 600 dpi plus PDF, geodesic scale bar, north arrow, °S/°W graticule, accessible palettes |

## Repository layout

```
ch4ne/                  Python package
  config.py             ALL parameters (the Methods can be written from here)
  gee_export.py         Earth Engine stage: orbits -> anomaly -> per-cell statistics -> assets -> download
  provenance.py         provenance record of the rasters and its verification
  boundaries.py         state boundaries (UF/name attributes) and coastline
  infra.py              harmonized and classified ANP inventory
  raster.py             GeoTIFF reading, rasterization, distances
  detection.py          local background, robust z, clusters, q-values, attribution
  zones.py              distance rings, block bootstrap, matching, regressions
  states.py             per-state table and independent check
  timeseries.py         regional trend and cluster time series
  figures.py            publication figures
  analysis.py           runs the local analysis and writes tables, figures and numbers
scripts/run_pipeline.py command-line interface
notebooks/              Google Colab notebook
tests/                  tests on real data (pytest)
data/anp/               public ANP pipeline layers used in the article
data/limites/           state boundaries (pinned copy; checksums in data/SHA256SUMS.txt)
validacao_offline/      reproduction of the v3 errors from the published tables
validacao_v4_1/         checks of v4.1 on real data (new vs old version, small area)
docs/                   audit, verified references, method notes, results (Portuguese)
outputs_v4/             outputs of the main run (tables, figures, numbers)
outputs_v4_P1/, outputs_v4_P2/   outputs before (P1) and after (P2) the 17 July 2022 product change
```

The Earth Engine rasters (`outputs_v4*/gee/*.tif`), `outputs_v4*/derived_fields.tif` and `outputs_v4*/pixels_strict.csv.gz` are not in Git. They are archived with the code in Zenodo (see *Data and code archive*).

## How to run

### 1. Environment

```bash
python -m pip install -r requirements.txt
```

`requirements.txt` gives minimum versions. The exact environment that produced the published numbers is recorded in `outputs_v4*/manuscript_numbers.json` (key `run`), and `requirements-lock.txt` is a `pip freeze` of the development environment (Python 3.14, Windows 11, 23 September 2026).

### 2. Earth Engine authentication (once per machine)

```bash
python -c "import ee; ee.Authenticate(auth_mode='localhost', force=True)"
```

The browser opens; sign in with the Google account registered for Earth Engine. The Google Cloud project is set in `ch4ne/config.py` and can be changed with `--project <id>` in any stage or with the environment variable `CH4NE_GEE_PROJECT`. The project must be registered at https://code.earthengine.google.com/register.

### 3. Processing (from the repository root)

```bash
python scripts/run_pipeline.py check
```

```bash
python scripts/run_pipeline.py all
```

- `check` (about 1 minute): bands, orbits, processor versions, one prepared orbit and the background of one month.
- `all` exports in two phases. **Phase A** writes the monthly backgrounds to a table asset. **Phase B** uses that asset and writes the image assets:
  - `stats`: all per-cell statistics, in a single pass over the orbits;
  - `monthly_2020` … `monthly_2025`: monthly anomaly and count, one asset per year;
  - `legacy`: the v3 composite, kept for comparison.

  Then come `wait`, `download` and `analyze`. `download` computes the static covariates (MODIS albedo, SRTM elevation and JRC water) directly on the grid, in blocks, with `ee.data.computePixels` (`covariates.tif`). It then assembles `stats.tif`, `monthly_anom.tif` and `monthly_count.tif`, plus the CSV and `PROVENANCE.json`.
- The export can be resumed. An existing asset with the same run identity is reused and a running task is awaited; running again never redoes finished work. `--jobs stats monthly` exports only the listed assets.
- Existing assets are never deleted without `--overwrite`. Even then, the old asset is replaced only when the new task finishes.
- If only the monthly series fail, the analysis still runs and omits only the cluster time series. The number of months with data already comes from `stats`.
- Stages can be run separately. `export` writes `run_config.json`, and later stages stop if the command line asks for different parameters.
- **Quota (Community Tier).** The cost of every task is logged in `outputs_v4*/gee/gee_tasks.csv`. In the main run (72 months) the successful tasks cost 8.2 EECU-h: `stats` 4.6 and the six `monthly_*` together 3.5. P1 cost 3.7 EECU-h and P2 8.7; the cost of `stats` varies considerably between executions. The `*_tile_px` parameters in `config.py` only control memory, never results. Measured failures: as a batch export, the JRC water covariate (30 m) exceeded the 2^31-pixel limit even with 32-px tiles (hence `computePixels`), and a single 144-band monthly asset ran out of memory.

**Robustness runs** by processor version (the label is mandatory for dates other than the default):

```bash
python scripts/run_pipeline.py all --start 2020-01-01 --end 2022-07-16 --label P1 --min-obs 22
```

```bash
python scripts/run_pipeline.py all --start 2022-07-17 --end 2025-12-31 --label P2 --min-obs 29
```

```bash
python scripts/run_pipeline.py compare --compare P1 P2
```

`--min-obs` keeps the main rule's density (50 retrievals in 72 months) in the shorter periods: P1 has 31 months and P2 42. With the absolute minimum of 50, only 6 cells near pipelines remain in P1. Re-running the main analysis removes Tables S12 and S19 of the comparison, so run `compare` again afterwards.

### 4. Tests

```bash
python -m pytest tests -q
```

### Google Colab

Open `notebooks/CH4_NE_v4_colab.ipynb`. It clones or mounts this repository, installs the dependencies, authenticates and runs the same stages.

## Outputs (`outputs_v4/`)

| File | Content |
|---|---|
| `tables/Table1_infrastructure_inventory.csv` | ANP segments by class and location |
| `tables/Table2_states_strict.csv` | per-state statistics (strict screen) |
| `tables/Table3_clusters.csv` | clusters: peak, centroid, state, z, q-value, distances, nearest asset |
| `tables/Table4_attribution_enrichment.csv` | observed vs expected peaks within 5/10/20 km of pipelines |
| `tables/Table5_proximity_rings_*.csv` | statistics by distance ring |
| `tables/Table6_near_vs_far_bootstrap_*.csv` | near vs far contrasts, 95% block-bootstrap intervals |
| `tables/Table7_near_vs_matched_controls.csv` | contrasts with matched controls |
| `tables/Table8_regressions_cluster_robust.csv` | regressions (mean and dispersion) with block-robust errors |
| `tables/TableS*.csv` | per-segment inventory, v3 reproduction, states (relaxed screen), matching balance, monthly background, cluster series, sensitivities (thresholds, window, block size), pipeline–coast collocation, P1 vs P2 comparison |
| `figures/Fig0*.png/.pdf` | figures (no titles; captions are in the article) |
| `manuscript_numbers.json` | every number cited in the text, with provenance, parameters and software versions of the run |
| `checks.json` | results of the automatic checks |
| `clusters_v4.geojson`, `derived_fields.tif` | GIS layers |

The table numbers in the file names are those of the pipeline, not of the article. Fig. 1 of the article is drawn without the locator inset (`figures.fig_study_area(..., inset=False)`).

## Data sources

- **TROPOMI/Sentinel-5P:** Google Earth Engine collection `COPERNICUS/S5P/OFFL/L3_CH4`, band `CH4_column_volume_mixing_ratio_dry_air_bias_corrected` (Earth Engine keeps only qa_value > 0.5). Contains modified Copernicus Sentinel data (2020–2025).
- **Pipelines:** public layers of the Agência Nacional do Petróleo, Gás Natural e Biocombustíveis (ANP): gas gathering, transport and distribution pipelines, and oil pipelines (https://geomaps.anp.gov.br/geoanp/). A copy is in `data/anp/`.
- **State boundaries:** `click_that_hood/brazil-states.geojson`, the same version used in v3 (pinned copy with checksum). It can be replaced by the official IBGE boundaries via `config.states_file`.
- **Covariates** (in Earth Engine): MODIS MCD43A3 (albedo, band 7), CGIAR SRTM90 (elevation), JRC Global Surface Water (water).

## Data and code archive

The code and the derived data of the three runs (main, P1, P2), including the Earth Engine rasters that are not in Git, are archived in Zenodo as **version 2.0.0** of the record that holds the superseded v1.0.0:

- all versions (concept DOI): https://doi.org/10.5281/zenodo.18717651
- version 2.0.0: <!-- add the version DOI reserved in Zenodo here, e.g. https://doi.org/10.5281/zenodo.NNNNNNNN -->

## How to cite

Please cite the article (reference above; it will be updated on publication) and the Zenodo archive, version 2.0.0. GitHub shows a ready-made citation under "Cite this repository", built from [CITATION.cff](CITATION.cff).

## License

- **Code:** MIT License (see [LICENSE](LICENSE)).
- **Derived data and figures** in `outputs_v4*/` and in the Zenodo archive: Creative Commons Attribution 4.0 International (CC BY 4.0). They contain modified Copernicus Sentinel data (2020–2025). The ANP pipeline layers in `data/anp/` are public data of the ANP and remain under the ANP's terms; cite the source when reusing them.

## Acknowledgments

This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior – Brasil (CAPES) – Finance Code 001; by the Human Resources Training Program of the ANP (PRH-ANP), managed by the São Paulo Research Foundation (FAPESP), grants 2025/17670-6, 2026/09019-6 and PRH-48/UFPE 2025/03556-7, with earlier support from PRH 48.1/UFPE (ANP/FINEP 48610.201019/2019-38; FAPESP 2024/10544-2); and by the Fundação de Amparo à Ciência e Tecnologia do Estado de Pernambuco (FACEPE), through the doctoral scholarship IBPG-1473-3.09/25 and the MAISSER-PE project (APQ-1428-3.01/25).
