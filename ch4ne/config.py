"""
Central configuration for the NE Brazil TROPOMI methane screening pipeline (v4).

Every number that shapes a result lives here, so the Methods section can be
written from this file and a sensitivity run only needs to override values.
"""
import math
import os
from dataclasses import asdict, dataclass, field
from datetime import date
from pathlib import Path

DEFAULT_START = "2020-01-01"
DEFAULT_END = "2025-12-31"


@dataclass
class Config:
    # ------------------------------------------------------------------ paths
    # Repository root; data are read from and outputs written below it.
    base_dir: Path = Path(".")
    infra_dir: str = "data/anp"
    # Output folder name, or an absolute path (e.g. a Google Drive folder).
    out_dir_name: str = "outputs_v4"
    # Label for alternative runs (e.g. "P2" = homogeneous processor period).
    # It separates the output folder and the Earth Engine asset folder, and is
    # mandatory whenever the dates differ from the defaults.
    run_label: str = ""

    # ------------------------------------------------------------ Earth Engine
    # Google Cloud project with Earth Engine enabled. Override with the
    # CH4NE_GEE_PROJECT environment variable or `--project` on the CLI.
    gee_project: str = field(default_factory=lambda: os.environ.get("CH4NE_GEE_PROJECT",
                                                                    "deep-presence-470415-u2"))
    collection: str = "COPERNICUS/S5P/OFFL/L3_CH4"
    band: str = "CH4_column_volume_mixing_ratio_dry_air_bias_corrected"
    uncertainty_band: str = "CH4_column_volume_mixing_ratio_dry_air_uncertainty"
    aod_band: str = "aerosol_optical_depth"
    # Inclusive analysis window. filterDate() is end-exclusive, so the code
    # adds one day to `end_date` (the v3 notebook silently dropped 31 Dec 2025).
    start_date: str = DEFAULT_START
    end_date: str = DEFAULT_END
    # Retrievals with reported uncertainty above this value (in ppb, NOT %)
    # are discarded. Same rule as v3 (`lte(15)`); it removes almost nothing,
    # because GEE keeps only qa_value > 0.5 at ingestion (noise <= 10 ppb).
    max_uncertainty_ppb: float = 15.0
    # Retrievals with SWIR aerosol optical depth above this are discarded
    # (Schuit et al., 2023, ACP 23, 9071).
    max_aod_swir: float = 0.13
    # Each L3 image is one orbit, and the footprints of all ~14 daily orbits
    # intersect the export window. Only orbits STARTING between these UTC
    # hours (inclusive) carry data over it: in 24 sample days spread over
    # 2020-2025 every orbit with data started at 14.96-17.11 UTC. The 2-h
    # margins keep all of them; dry_run re-checks this on sample days.
    overpass_start_hour_utc: int = 13
    overpass_end_hour_utc: int = 18

    # ------------------------------------------------------------ analysis grid
    # Regular lat/lon grid, ~6.9 km at the equator (close to the 5.5 x 7 km
    # TROPOMI footprint). Values are AGGREGATED onto it with reduceResolution
    # (v3 used reproject() alone: nearest-neighbour sampling of a pyramid level).
    grid_res_deg: float = 0.0625
    # Export window = NE Brazil + margin, so the local-background filter has
    # data beyond the state borders (ocean pixels are masked later).
    grid_west: float = -50.0
    grid_east: float = -34.0
    grid_south: float = -20.0
    grid_north: float = 0.5
    # A grid cell of one orbit is kept only if at least this fraction of the
    # native 0.01 deg L3 bins inside it hold a valid retrieval.
    min_cell_valid_fraction: float = 0.3
    # Months whose regional background rests on fewer valid NE cells than this
    # are dropped. Coverage is strongly seasonal (measured: ~1,000-4,700 cells
    # in Jan-Apr vs 21,000-26,000 in Jul-Oct of ~32,750), so the limit is low
    # enough to keep wet-season months; the latitude of the observed cells is
    # recorded per month to diagnose coverage bias.
    min_bg_cells_per_month: int = 500
    # Earth Engine memory controls: raw-value buffer of the per-cell
    # percentile reducer (exact below this many retrievals per cell),
    # ImageCollection.reduce parallelScale, and reduceRegion tileScale.
    percentile_max_raw: int = 3000
    parallel_scale: int = 4
    tile_scale: int = 4
    # Tile size (pixels) in which each asset export is computed, and block
    # size of the interactive covariate computation; they change memory use
    # only, never results. The per-cell statistics hold up to
    # percentile_max_raw values per cell; the inland-water covariate reads
    # ~54,000 JRC 30 m pixels per grid cell (a 32-cell block takes ~20 s).
    stats_tile_px: int = 128
    covariates_tile_px: int = 32
    monthly_tile_px: int = 128

    # ----------------------------------------------------- boundaries & infra
    states_url: str = ("https://raw.githubusercontent.com/codeforamerica/click_that_hood/"
                       "master/public/data/brazil-states.geojson")
    # Pinned copy of the boundary file used in v3 (see data/SHA256SUMS.txt);
    # set to "" to download states_url, or point to an IBGE BR_UF shapefile.
    states_file: str = "data/limites/brazil-states.geojson"
    ne_states: tuple = ("MA", "PI", "CE", "RN", "PB", "PE", "AL", "SE", "BA")
    infra_files: dict = field(default_factory=lambda: {
        "gasoduto-escoamento.shp": ("Gathering gas pipeline", "Nome_Dut_1"),
        "gasoduto-distribuicao.shp": ("Distribution gas pipeline", "Distrib"),
        "gasoduto-transporte.shp": ("Transport gas pipeline", "Nome_Dut_1"),
        "oleoduto.shp": ("Oil pipeline", "Nome_Dut_1"),
    })
    # A segment is "onshore NE" if at least this share of its length lies on
    # NE land. Offshore and out-of-region segments (e.g. Espirito Santo) are
    # kept in the inventory table but excluded from proximity analysis.
    onshore_share_min: float = 0.5
    simplify_tolerance_deg: float = 1e-5

    # ------------------------------------------------------ pixel screening
    # Minimum number of valid orbit-level retrievals per cell (absolute; it
    # bounds the standard error of the cell mean), and minimum share of the
    # months in the window with at least one retrieval (24 of 72 months).
    min_obs: int = 50
    min_valid_month_fraction: float = 1 / 3
    # Cells whose centre lies closer than this to the coastline, or whose
    # land fraction / inland-water fraction fail the limits, are excluded from
    # the "strict" analysis AND from the local-background estimate.
    coast_buffer_km: float = 10.0
    min_land_fraction: float = 0.9
    max_water_fraction: float = 0.2

    # -------------------------------------------------------- anomaly detection
    # Local background: nan-median in a square window of ~2 deg (Schuit et al.
    # 2023 and related work); 1, 3 and 5 deg are run as sensitivity cases.
    local_bg_window_km: float = 220.0
    sensitivity_windows_km: tuple = (110.0, 220.0, 330.0, 550.0)
    z_threshold: float = 2.5
    min_cluster_pixels: int = 2
    eight_connected: bool = True
    q_threshold: float = 0.10

    # ------------------------------------------------------ proximity analysis
    ring_edges_km: tuple = (0.0, 5.0, 15.0, 50.0, float("inf"))
    near_ring: str = "0-5 km"
    far_ring: str = ">=50 km"
    attribution_radius_km: tuple = (5.0, 10.0, 20.0)
    # The positive-vs-mirror comparison of peak distances is also stratified by
    # the peak's distance to the coast (strata split at these values, km).
    attribution_coast_strata_km: tuple = (50.0,)
    # Spatial blocks for the bootstrap and the cluster-robust errors; at least
    # as wide as the background window (pixels in adjacent 1-deg blocks share
    # most of their background). 1 and 3 deg are reported as sensitivity.
    bootstrap_block_deg: float = 2.0
    sensitivity_blocks_deg: tuple = (1.0, 2.0, 3.0)
    n_boot: int = 2000
    n_matched_controls: int = 5
    # Primary matching: Mahalanobis distance on the covariates among the far
    # cells whose propensity-score logit (logistic regression of near vs far
    # on the covariates) lies within 0.2 SD of the near cell's (Rosenbaum &
    # Rubin 1985; Austin 2011, Pharm. Stat. 10:150). Balance is reported as
    # standardised mean differences (TableS4).
    match_method: str = "ps_mahalanobis"
    match_caliper_ps_sd: float = 0.2
    # Sensitivity: Euclidean distance in the standardised covariate space with
    # this caliper (in 9 dimensions it discarded 158 of 168 near cells).
    match_caliper: float = 0.5
    random_seed: int = 42

    # ----------------------------------------------------------------- figures
    dpi: int = 600
    fig_width_full_in: float = 6.85  # EMAS (Springer) double-column width, 174 mm
    fig_width_col_in: float = 3.4
    font_size: float = 10.0   # EMAS: lettering 8-12 pt (smallest labels = font_size - 2)

    # --------------------------------------------- legacy (v3) reproduction
    legacy_scale_m: float = 7000.0
    published_mu: float = 1880.14
    published_sigma: float = 9.39
    published_n: int = 32904

    # ------------------------------------------------------------ helpers
    @property
    def out_dir(self) -> Path:
        suffix = f"_{self.run_label}" if self.run_label else ""
        return Path(self.base_dir) / f"{self.out_dir_name}{suffix}"

    @property
    def n_months(self) -> int:
        a, b = date.fromisoformat(self.start_date), date.fromisoformat(self.end_date)
        return (b.year - a.year) * 12 + (b.month - a.month) + 1

    @property
    def min_valid_months(self) -> int:
        return int(math.ceil(self.min_valid_month_fraction * self.n_months - 1e-9))

    @property
    def is_default_period(self) -> bool:
        return (self.start_date, self.end_date) == (DEFAULT_START, DEFAULT_END)

    @property
    def period_label(self) -> str:
        return f"{self.start_date[:7]} to {self.end_date[:7]}"

    def validate(self):
        a, b = date.fromisoformat(self.start_date), date.fromisoformat(self.end_date)
        if a >= b:
            raise ValueError(f"start_date {a} must precede end_date {b}")
        if not self.is_default_period and not self.run_label:
            raise ValueError("A non-default period needs a run_label (CLI: --label), so it cannot "
                             "overwrite the main run's outputs or Earth Engine assets.")
        return self

    def gee_params(self) -> dict:
        """Settings that shape the Earth Engine outputs (stamped on every asset)."""
        keys = ("collection", "band", "uncertainty_band", "aod_band", "start_date", "end_date",
                "max_uncertainty_ppb", "max_aod_swir", "overpass_start_hour_utc", "overpass_end_hour_utc",
                "grid_res_deg", "grid_west", "grid_east",
                "grid_south", "grid_north", "min_cell_valid_fraction", "min_bg_cells_per_month",
                "percentile_max_raw", "run_label")
        d = asdict(self)
        return {k: (list(d[k]) if isinstance(d[k], tuple) else d[k]) for k in keys}

    def as_dict(self):
        d = asdict(self)
        d["base_dir"] = str(self.base_dir)
        d["ring_edges_km"] = [str(x) for x in self.ring_edges_km]
        d["min_valid_months"] = self.min_valid_months
        d["n_months"] = self.n_months
        return d


# Value written for masked pixels in every Earth Engine export (read back as NaN).
NODATA = -9999.0

# Band order of the per-pixel statistics GeoTIFF (stats.tif). Earth Engine
# writes two assets, the per-cell statistics of the retrievals
# (STATS_CORE_BANDS, one pass over the orbits) and the static surface
# covariates (COVARIATE_BANDS); download_all() stacks them in this order and
# the local reader relies on it when the GeoTIFF carries no band descriptions.
STATS_CORE_BANDS = (
    "n_obs",            # valid orbit-level retrievals per cell
    "xch4_mean",        # raw bias-corrected XCH4, mean of cell-aggregated retrievals (ppb)
    "anom_mean",        # de-trended anomaly (retrieval - monthly NE median), mean (ppb)
    "anom_median",
    "anom_sd",          # TEMPORAL standard deviation of the anomaly (ppb)
    "anom_p10",
    "anom_p90",
    "anom_p05",
    "anom_p95",
    # Tail widths of each cell's anomaly distribution. Intermittent emission
    # widens only the UPPER tail; retrieval noise and artefacts widen both.
    # (Computed from the percentiles in the same pass, so the prepared orbits
    # are not re-read; they replace v4-beta exceedance frequencies.)
    "tail_up_90",       # p90 - p50
    "tail_lo_10",       # p50 - p10
    "tail_up_95",       # p95 - p50
    "tail_lo_05",       # p50 - p05
    "mean_t_years",     # mean observation time, years since start (sampling diagnostic)
    "frac_janjun",      # share of retrievals in Jan-Jun (wet-season sampling diagnostic)
    "aod_mean",         # mean aerosol optical depth of retained retrievals
    "n_months",         # number of distinct calendar months with >= 1 retrieval
)
COVARIATE_BANDS = (
    "albedo_b7",        # MODIS MCD43A3 black-sky albedo, band 7 (2.1 um)
    "elev_mean",        # SRTM elevation mean (m)
    "elev_sd",          # SRTM elevation SD within cell (m), terrain roughness
    "water_frac",       # JRC GSW occurrence>=50% fraction (inland water)
)
STATS_BANDS = STATS_CORE_BANDS + COVARIATE_BANDS
