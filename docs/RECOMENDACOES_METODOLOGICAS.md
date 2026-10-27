# Recomendações metodológicas com base na literatura

Levantamento feito em 23/09/2026 para substituir o método da v3 (Z > 2 com média e desvio globais sobre a média bruta de 6 anos). Cada item traz os parâmetros sugeridos e as referências; o pipeline v4 já implementa os itens marcados em `ch4ne/config.py`. Confira as referências antes de citá-las.

## 1. 0. Input product: the GEE L3 OFFL CH4 collection is not homogeneous across processor versions (prerequisite for all other topics)

Describe COPERNICUS/S5P/OFFL/L3_CH4 as a screening product and state its limits in the Methods. GEE builds it with harpconvert bin_spatial on a 0.01° grid and applies the filter CH4_column_volume_mixing_ratio_dry_air_validity>50, which is the same as qa_value>0.5. The collection contains only XCH4, bias-corrected XCH4, precision (the 'uncertainty' band), SWIR aerosol optical depth, aerosol height and viewing angles. It has no SWIR surface albedo, no land_fraction or surface_classification, no chi-square and no destriped XCH4. It is OFFL only; the 2022 full-mission reprocessing (RPRO v02.04.00) is not listed.

Because of this, the 2020–2025 'bias_corrected' band mixes several algorithms:
- Jan 2020 to Jun 2021 (v01.03.02 and v01.04.00): the a-posteriori correction was referenced to GOSAT.
- From 1 Jul 2021 (v02.02.00): the correction uses TROPOMI only (the Lorente et al. 2021 B-spline in SWIR albedo).
- From 14 Nov 2021 (v02.03.01): ocean sun-glint retrievals are included. GEE also notes that ocean data start in November 2021.
- From 17 Jul 2022 (v02.04.00): the third-order surface-reflectance polynomial removes surface-spectral artefacts of 20–40 ppb (Lorente et al. 2023). Earlier data do not have this fix.
- From 12 Mar 2023 (v02.05.00): anomalously low values at land–water boundaries become especially common.

Required actions:
(i) Run a version-split robustness test. Build independent maps for P1 (2020-01-01 to 2022-07-16) and P2 (2022-07-17 to 2025-12-31) and keep only anomalies present in both. Optionally add a split at 2021-11-14.
(ii) For the final inference, or at least one sensitivity run, use a homogeneous Level-2 record. Either follow the product readme (RPRO v02.04.00 up to mid-July 2022 plus OFFL v02.04 or later afterwards, with qa_value>0.5, the land_fraction filter and the destriped key), or use the blended TROPOMI+GOSAT product, as Hancock et al. (2025) did for South America.
(iii) Compare regional TROPOMI with GOSAT over Northeast Brazil before interpreting anomalies, which Jacob et al. (2022) call good practice.
(iv) Aggregate the 0.01° bins to the analysis grid explicitly with reduceResolution(ee.Reducer.mean()) plus a count band, instead of relying on the implicit resampling of .reproject(scale=7000).

**Parâmetros:** L3 grid 0.01° (≈1113 m). Implicit filter qa_value>0.5. Precision filter <10 ppb: the manuscript's `uncertainty ≤ 15 ppb` has no effect, because qa>0.5 already rejects CH4 noise error >10 ppb. Processor cut dates: 2021-07-01 (v02.02.00), 2021-11-14 (v02.03.01, ocean glint), 2022-07-17 (v02.04.00, surface-spectral fix), 2023-03-12 (v02.05.00, land–water low bias), 2023-11-26 (v02.06.00), 2024-09-07 (v02.07.01, new destriping and Copernicus 90 m DEM), 2024-11-16 (v02.08.00), 2025-11-22 (v02.09.01, land_fraction key). Version-split periods: P1 2020-01-01 to 2022-07-16; P2 2022-07-17 to 2025-12-31. Analysis grid: 0.05°–0.1° (or 7 km) by mean aggregation, with a per-cell count band.

**Notas:** Sources checked: the PRF text (scratchpad/prf_ch4.txt, Table 2 and sections 3.1 and 4.1) and the GEE catalog page, fetched 2026-09-23. The GEE page states that qa_value filtering 'does not remove all pixels considered bad'. The PRF adds that pixels with too-low CH4 remain after qa>0.5; these low outliers inflate standard-deviation statistics, which argues for robust statistics (topic 5). The PRF says a total uncertainty estimate should be taken as 2× the reported single-sounding precision. Jacob et al. (2022) give Lorente et al.'s variable bias as 5.6 ppb against TCCON, but also show TROPOMI–GOSAT regional differences above 20 ppb in some regions. Both are comparable to the manuscript's threshold of about +19 ppb (μ+2σ with σ = 9.4 ppb). This check could not confirm whether GEE ingests any RPRO files; the image metadata (processor version) should be inspected when stratifying by version.

**Referências:**

- Landgraf, J., Lorente, A., Borsdorff, T., Mandal, S., Langerock, B., Kumar Sha, M. (2025). S5P Mission Performance Cluster Methane [L2__CH4___] Product Readme File, S5P-MPC-SRON-PRF-CH4, issue 2.9, product version V02.09.01, 27 Nov 2025 (product DOI). https://doi.org/10.5270/S5P-3lcdqiv
- Google Earth Engine Data Catalog: Sentinel-5P OFFL CH4: Offline Methane, COPERNICUS/S5P/OFFL/L3_CH4 (accessed 2026-09-23). https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S5P_OFFL_L3_CH4
- Lorente, A., Borsdorff, T., Butz, A., Hasekamp, O., aan de Brugh, J., Schneider, A., et al. (2021). Methane retrieved from TROPOMI: improvement of the data product and validation of the first 2 years of measurements. Atmos. Meas. Tech., 14, 665–684. https://doi.org/10.5194/amt-14-665-2021
- Lorente, A., Borsdorff, T., Martinez-Velarte, M. C., Landgraf, J. (2023). Accounting for surface reflectance spectral features in TROPOMI methane retrievals. Atmos. Meas. Tech., 16, 1597–1608. https://doi.org/10.5194/amt-16-1597-2023
- Balasus, N., Jacob, D. J., Lorente, A., Maasakkers, J. D., Parker, R. J., Boesch, H., et al. (2023). A blended TROPOMI+GOSAT satellite data product for atmospheric methane using machine learning to correct retrieval biases. Atmos. Meas. Tech., 16, 3787–3807. https://doi.org/10.5194/amt-16-3787-2023
- Jacob, D. J., Varon, D. J., Cusworth, D. H., Dennison, P. E., Frankenberg, C., Gautam, R., et al. (2022). Quantifying methane emissions from the global scale down to point sources using satellite observations of atmospheric methane. Atmos. Chem. Phys., 22, 9617–9646. https://doi.org/10.5194/acp-22-9617-2022
- Hancock, S. E., Jacob, D. J., Chen, Z., Nesser, H., Davitt, A., Varon, D. J., et al. (2025). Satellite quantification of methane emissions from South American countries: a high-resolution inversion of TROPOMI and GOSAT observations. Atmos. Chem. Phys., 25, 797–817. https://doi.org/10.5194/acp-25-797-2025

## 2. 1. Local background removal (replace the global μ/σ of the whole region)

Replace Z = (X − μ_NE)/σ_NE, computed on the 6-year mean, with a local anomaly map:

ΔXCH4(x,t) = XCH4_ec(x,t) − median{XCH4_ec in window W(x)}

XCH4_ec is XCH4 corrected for elevation. Compute ΔXCH4 on monthly (or per-orbit) maps before compositing (see topic 2). Use the median as the background statistic. A robust local spread s_W(x) comes from 1.4826·MAD within W, or, following Vanselow et al. (2024), from the SD within W after excluding values above the 95th percentile.

Window sizes in the literature:
- Lauvaux et al. (2022): 11×11-pixel moving patch. Background is the SExtractor-type mode, 2.5·median − 1.5·mean, or the median when (mean − median)/σ > 0.3. Candidates are contiguous pixels above background + 3σ_patch.
- Schuit et al. (2023): 32×32-pixel scenes, about 176×232 km at nadir. Each scene is normalized between mean − 1σ and mean + 100 ppb − 1σ. For the integrated mass enhancement (IME), background is the median of the scene outside the plume mask.
- Vanselow et al. (2024): high-pass filter areas of 1°–5°, median background, at least 25% of the area filled.
- Barré et al. (2021): Gaussian high-pass filter with a 2.0° length scale over 30-day windows, applied to observation-minus-forecast departures.
- Liu et al. (2021, 2024): background = mean of the lowest 10% of values within ±5 cells (0.25° grid) or ±3 cells (0.2° grid), roughly 140–275 km.
- de Gouw et al. (2020): daily 10th percentile within the basin box.
- Sadavarte et al. (2021): mean of a 0.5°×0.5° box centred 0.1° upwind, needing at least 20 observations; otherwise the domain median.
- Buchwitz et al. (2017): source region minus the mean of a surrounding region, on annual maps.

For a regional screen at 7–11 km, use a 2°×2° median window (about 220 km) as the primary setting and report sensitivity for 1°, 3° and 5°. The window must be much larger than the clusters (1–3 cells) and smaller than synoptic or latitudinal gradients.

Apply an elevation correction before background removal. Northeast Brazil has plateaus of 0.5–1.5 km (Borborema, Araripe, Ibiapaba, Chapada Diamantina). At 7–8.5 ppb/km this gives an artificial contrast of about 3.5–13 ppb between plateaus and coastal lowlands, i.e. 0.4–1.4 of the manuscript's σ. This alone can raise the global Z-scores of lowland coastal cells, where the pipelines are.

For screening, do not use low-percentile backgrounds (p10 or p25). They are meant for quantifying fluxes in source-dominated windows and shift every anomaly upward.

**Parâmetros:** Primary window: 2°×2° median, about 220 km. In GEE: ee.Kernel.square(radius=110000, units='meters') with reduceNeighborhood(ee.Reducer.median()). Minimum fill in the window: ≥25% valid cells. Sensitivity windows: 1°, 3°, 5°. Robust spread: 1.4826·MAD within the window, or SD after excluding values above p95. Elevation correction: +7 ppb per km above mean sea level (Buchwitz 2017; Sadavarte 2021) or +8.5 ppb/km (Vanselow 2024); alternatively a regression on surface pressure fitted to the 25th percentile (Veefkind 2023). Report the literature anchors: Lauvaux 11×11 px and 3σ; Schuit 32×32 px; Barré 2° Gaussian over 30 days; Liu ±3–5 cells and lowest 10%.

**Notas:** The Lauvaux SI equations were lost in text extraction. The coefficient α = 2.5, the skewness criterion 0.3 and k = 3 are given in the SI text; the '>25 ppb averaged over several pixels' acceptance rule is stated in the main text. Zhang et al. (2020) used an elevation correction with a third-order polynomial, following Kort et al. (2014), and also filtered pixels with elevation SD ≥80 m within 5 km. Vanselow et al. (2024) is the closest peer-reviewed analogue to what the manuscript attempts: monthly 0.1° TROPOMI/WFMD maps for 2018–2021, local median high-pass filtering, 2σ enhancement plus a persistence rule. Citing and adopting it gives reviewers a recognizable benchmark. The Liu et al. (2021) background rule is taken from a secondary summary; the Wiley full text returned 403.

**Referências:**

- Lauvaux, T., Giron, C., Mazzolini, M., d'Aspremont, A., Duren, R., Cusworth, D., Shindell, D., Ciais, P. (2022). Global assessment of oil and gas methane ultra-emitters. Science, 375, 557–561 (Supplementary Materials, section 2.1). https://doi.org/10.1126/science.abj4351
- Schuit, B. J., Maasakkers, J. D., Bijl, P., Mahapatra, G., van den Berg, A.-W., Pandey, S., et al. (2023). Automated detection and monitoring of methane super-emitters using satellite data. Atmos. Chem. Phys., 23, 9071–9098. https://doi.org/10.5194/acp-23-9071-2023
- Vanselow, S., Schneising, O., Buchwitz, M., Reuter, M., Bovensmann, H., Boesch, H., Burrows, J. P. (2024). Automated detection of regions with persistently enhanced methane concentrations using Sentinel-5 Precursor satellite data. Atmos. Chem. Phys., 24, 10441–10473. https://doi.org/10.5194/acp-24-10441-2024
- Barré, J., Aben, I., Agustí-Panareda, A., Balsamo, G., Bousserez, N., Dueben, P., et al. (2021). Systematic detection of local CH4 anomalies by combining satellite measurements with high-resolution forecasts. Atmos. Chem. Phys., 21, 5117–5136. https://doi.org/10.5194/acp-21-5117-2021
- de Gouw, J. A., Veefkind, J. P., Roosenbrand, E., Dix, B., Lin, J. C., Landgraf, J., Levelt, P. F. (2020). Daily satellite observations of methane from oil and gas production regions in the United States. Sci. Rep., 10, 1379. https://doi.org/10.1038/s41598-020-57678-4
- Zhang, Y., Gautam, R., Pandey, S., Omara, M., Maasakkers, J. D., Sadavarte, P., et al. (2020). Quantifying methane emissions from the largest oil-producing basin in the United States from space. Sci. Adv., 6, eaaz5120. https://doi.org/10.1126/sciadv.aaz5120
- Sadavarte, P., Pandey, S., Maasakkers, J. D., Lorente, A., Borsdorff, T., Denier van der Gon, H., Houweling, S., Aben, I. (2021). Methane emissions from superemitting coal mines in Australia quantified using TROPOMI satellite observations. Environ. Sci. Technol., 55, 16573–16580. https://doi.org/10.1021/acs.est.1c03976
- Buchwitz, M., Schneising, O., Reuter, M., Heymann, J., Krautwurst, S., Bovensmann, H., et al. (2017). Satellite-derived methane hotspot emission estimates using a fast data-driven method. Atmos. Chem. Phys., 17, 5751–5774. https://doi.org/10.5194/acp-17-5751-2017
- Liu, M., van der A, R., van Weele, M., Eskes, H., Lu, X., Veefkind, P., et al. (2021). A new divergence method to quantify methane emissions using observations of Sentinel-5P TROPOMI. Geophys. Res. Lett., 48, e2021GL094151. https://doi.org/10.1029/2021GL094151
- Liu, M., van der A, R., van Weele, M., Bryan, L., Eskes, H., Veefkind, P., et al. (2024). Current potential of CH4 emission estimates using TROPOMI in the Middle East. Atmos. Meas. Tech., 17, 5261–5277. https://doi.org/10.5194/amt-17-5261-2024
- Kort, E. A., Frankenberg, C., Costigan, K. R., Lindenmaier, R., Dubey, M. K., Wunch, D. (2014). Four corners: The largest US methane anomaly viewed from space. Geophys. Res. Lett., 41, 6898–6903. https://doi.org/10.1002/2014GL061503
- Veefkind, J. P., Serrano-Calvo, R., de Gouw, J., Dix, B., Schneising, O., Buchwitz, M., et al. (2023). Widespread frequent methane emissions from the oil and gas industry in the Permian Basin. J. Geophys. Res. Atmos., 128, e2022JD037479. https://doi.org/10.1029/2022JD037479

## 3. 2. Remove the temporal trend and seasonality before compositing

Do not average raw retrievals over 2020–2025. The NOAA global annual mean rose from 1878.7 ppb (2020) to 1935.9 ppb (2025), about +57 ppb, with growth rates of 14.8, 17.7, 13.0, 8.3, 7.9 and 5.1 ppb/yr for 2020–2025. That is about six times the manuscript's regional σ of 9.4 ppb.

The number of valid retrievals per cell varies by year and season, and differently from place to place: wet-season clouds, ocean-glint data only after November 2021, and version changes. A per-pixel mean over all images is therefore weighted by sampling. A cell sampled on average one year later than its neighbours gets about +8–15 ppb, roughly 1σ, from the trend alone.

Compute anomalies against a background from the same time step, then composite with equal weight per time unit:
(a) Default, following Vanselow et al. (2024): build monthly maps. Keep a cell-month only if it has more than 3 days of data. Apply the local median high-pass of topic 1 to each month. The final map is the median of monthly ΔXCH4. Also report the persistence frequency F_enh = N_enh/N_meas.
(b) Sensitivity, per orbit or day: subtract the median of the same day within the 2° window (as in Lauvaux's per-orbit tiles and Liu et al.'s daily background, which needs at least 10 retrievals), or use 30-day windows as in Barré et al. (2021).

In GEE, a cheap first step is to subtract each image's domain median (Northeast Brazil plus a 1° buffer) before any temporal reduction. This removes the growth rate, the seasonal cycle and orbit-level offsets. Then apply the spatial high-pass. Also publish annual anomaly maps and year-to-year persistence, as Veefkind et al. (2023) did for 2019 and 2020. Episodic sources raise F_enh in specific years; artefacts tied to a processor version switch on or off at the version dates (topic 0).

**Parâmetros:** Time step: monthly (primary); daily/orbit (sensitivity). Monthly validity: N_days > 3 per cell-month. Composite: median of monthly ΔXCH4 with equal weight per month, not a mean over all retrievals. Persistence: F_enh ≥ 0.5 of valid months. Minimum months: N_meas ≥ 24 of 72, scaled from Vanselow's 16 of 48. Detrending shortcut: subtract the per-image domain median (Northeast Brazil plus 1°) before temporal reduction. NOAA global means for context: 2020 1878.72; 2021 1894.79; 2022 1910.95; 2023 1921.40; 2024 1929.41; 2025 1935.94 ppb.

**Notas:** The NOAA values are surface flask global means, taken from gml.noaa.gov/webdata/ccgg/trends/ch4/ch4_annmean_gl.txt and ch4_gr_gl.txt on 2026-09-23. The XCH4 trend has a similar magnitude, but the NOAA numbers are there only to show scale. Jacob et al. (2022) document seasonal sampling bias in TROPOMI and sparse coverage in the wet tropics. Hancock et al. (2025) note that TROPOMI data are almost absent over the Amazon because of clouds and dark surfaces, which matters for the Maranhão/Pará frontier cells in the manuscript's top-20. The manuscript's monthly time series (Figure 8) uses the monthly median of the domain. Reusing that same monthly domain median as the per-month offset gives a consistent, documented detrending step.

**Referências:**

- Lan, X., Thoning, K. W., Dlugokencky, E. J. Trends in globally-averaged CH4, N2O, and SF6 determined from NOAA Global Monitoring Laboratory measurements (global annual means and growth rates; file created 5 Sep 2026, accessed 23 Sep 2026). https://doi.org/10.15138/P8XG-AA10
- Vanselow, S., Schneising, O., Buchwitz, M., Reuter, M., Bovensmann, H., Boesch, H., Burrows, J. P. (2024). Automated detection of regions with persistently enhanced methane concentrations using Sentinel-5 Precursor satellite data. Atmos. Chem. Phys., 24, 10441–10473. https://doi.org/10.5194/acp-24-10441-2024
- Barré, J., Aben, I., Agustí-Panareda, A., Balsamo, G., Bousserez, N., Dueben, P., et al. (2021). Systematic detection of local CH4 anomalies by combining satellite measurements with high-resolution forecasts. Atmos. Chem. Phys., 21, 5117–5136. https://doi.org/10.5194/acp-21-5117-2021
- Lauvaux, T., et al. (2022). Global assessment of oil and gas methane ultra-emitters. Science, 375, 557–561. https://doi.org/10.1126/science.abj4351
- Liu, M., van der A, R., van Weele, M., Bryan, L., Eskes, H., Veefkind, P., et al. (2024). Current potential of CH4 emission estimates using TROPOMI in the Middle East. Atmos. Meas. Tech., 17, 5261–5277. https://doi.org/10.5194/amt-17-5261-2024
- Veefkind, J. P., Serrano-Calvo, R., de Gouw, J., Dix, B., Schneising, O., Buchwitz, M., et al. (2023). Widespread frequent methane emissions from the oil and gas industry in the Permian Basin. J. Geophys. Res. Atmos., 128, e2022JD037479. https://doi.org/10.1029/2022JD037479
- Jacob, D. J., et al. (2022). Quantifying methane emissions from the global scale down to point sources using satellite observations of atmospheric methane. Atmos. Chem. Phys., 22, 9617–9646. https://doi.org/10.5194/acp-22-9617-2022

## 4. 3. Minimum number of valid retrievals per cell

Mask cells with insufficient sampling and publish a map of the count of valid observations.

Count distinct orbit-days with valid data per analysis cell, not 0.01° L3 bins. One 5.5×7 km TROPOMI footprint fills about 35 bins of 0.01°, so bin counts overstate the number of independent samples roughly 35-fold.

Thresholds in the literature:
- Vanselow et al. (2024): a monthly 0.1° cell needs N_days > 3; a region needs data in at least 16 of 48 months; a high-pass area must be at least 25% filled.
- Hancock et al. (2025): drop 0.25°×0.3125° cells with fewer than 30 individual TROPOMI retrievals in the year.
- Liu et al. (2024): at least 10 valid retrievals to build the daily background, and more than 7 valid days per source.
- Barré et al. (2021): N ≥ 2 observations per 0.1° bin per 30-day window.
- Schuit et al. (2023): scenes need more than 20% valid pixels.
- Lauvaux et al. (2022): coverage statistics on 0.05° pixels need at least 10 valid measurements per year; plumes need at least 5 pixels with qa>0.75.
- Sadavarte et al. (2021): orbits need more than 500 valid pixels in a 4°×4° domain, and the upwind background needs at least 20 observations.
- Schneising et al. (2020): at least 25% coverage of the background, plume and hotspot regions, and at least 60% for each cross-section.

For this manuscript, use:
(i) monthly cell valid if N_days > 3;
(ii) cell kept in the composite if it has at least 24 valid months and at least 30 valid retrieval-days, with data in at least 4 of the 6 years and in both the wet and dry seasons;
(iii) background window at least 25% valid.

Statistical rationale: single-sounding precision against TCCON is about 14.5 ppb (Balasus et al. 2023). The standard error of a cell mean is about σ/√N, so a standard error of 2.5 ppb or less needs N ≳ 34 if errors were independent and identically distributed. That assumption is optimistic (Jacob et al. 2022), and averaging never reduces the variable bias (5.6 ppb against TCCON; more than 20 ppb against GOSAT in some regions). So N sets only a lower bound on detectability.

**Parâmetros:** Monthly cell: N_days > 3. Composite cell: ≥ 24 valid months (of 72) and ≥ 30 valid orbit-days, coverage in ≥ 4 of 6 years and in both seasons. Background window: ≥ 25% valid cells, or ≥ 10–20 retrievals for daily or upwind backgrounds (Liu 2024; Sadavarte 2021). Cluster: ≥ 2 contiguous cells and ≥ 5 valid pixels. Count unit: distinct orbit-days per cell, not 0.01° bins. Expected precision: σ_single ≈ 12–15 ppb, so N ≳ 34 for SE ≤ 2.5 ppb under independent errors.

**Notas:** Items (ii) and (iii) of the recommendation are my synthesis: Vanselow's 16/48 months scaled to 72 months, and Hancock's ≥30 retrievals per 25-km cell per year used as a floor. They are not a single published rule and should be presented as a design choice, with a sensitivity analysis over N_min ∈ {15, 30, 60} orbit-days. In GEE, count per cell with ee.ImageCollection.map(img → img.select(band).mask().reduceResolution(ee.Reducer.max()).reproject(grid)).sum(). This counts days on which the cell has any valid bin.

**Referências:**

- Vanselow, S., et al. (2024). Automated detection of regions with persistently enhanced methane concentrations using Sentinel-5 Precursor satellite data. Atmos. Chem. Phys., 24, 10441–10473. https://doi.org/10.5194/acp-24-10441-2024
- Hancock, S. E., Jacob, D. J., Chen, Z., et al. (2025). Satellite quantification of methane emissions from South American countries: a high-resolution inversion of TROPOMI and GOSAT observations. Atmos. Chem. Phys., 25, 797–817. https://doi.org/10.5194/acp-25-797-2025
- Liu, M., et al. (2024). Current potential of CH4 emission estimates using TROPOMI in the Middle East. Atmos. Meas. Tech., 17, 5261–5277. https://doi.org/10.5194/amt-17-5261-2024
- Barré, J., et al. (2021). Systematic detection of local CH4 anomalies by combining satellite measurements with high-resolution forecasts. Atmos. Chem. Phys., 21, 5117–5136. https://doi.org/10.5194/acp-21-5117-2021
- Schuit, B. J., et al. (2023). Automated detection and monitoring of methane super-emitters using satellite data. Atmos. Chem. Phys., 23, 9071–9098. https://doi.org/10.5194/acp-23-9071-2023
- Lauvaux, T., et al. (2022). Global assessment of oil and gas methane ultra-emitters. Science, 375, 557–561. https://doi.org/10.1126/science.abj4351
- Sadavarte, P., et al. (2021). Methane emissions from superemitting coal mines in Australia quantified using TROPOMI satellite observations. Environ. Sci. Technol., 55, 16573–16580. https://doi.org/10.1021/acs.est.1c03976
- Schneising, O., Buchwitz, M., Reuter, M., Vanselow, S., Bovensmann, H., Burrows, J. P. (2020). Remote sensing of methane leakage from natural gas and petroleum systems revisited. Atmos. Chem. Phys., 20, 9169–9182. https://doi.org/10.5194/acp-20-9169-2020
- Balasus, N., et al. (2023). A blended TROPOMI+GOSAT satellite data product for atmospheric methane using machine learning to correct retrieval biases. Atmos. Meas. Tech., 16, 3787–3807. https://doi.org/10.5194/amt-16-3787-2023
- Jacob, D. J., et al. (2022). Quantifying methane emissions from the global scale down to point sources using satellite observations of atmospheric methane. Atmos. Chem. Phys., 22, 9617–9646. https://doi.org/10.5194/acp-22-9617-2022

## 5. 4a. Artefacts at coasts and water bodies: masks usable in GEE

Coasts and large water bodies are the dominant, documented artefact mechanism, and they also confound the pipeline proximity analysis.

Evidence from the literature:
- The PRF (issue 2.9) says qa>0.5 includes both land and ocean pixels since v2.3.1. Low CH4 values occur often at land–water boundaries, especially since v02.05.00 (12 Mar 2023). It recommends excluding pixels with 0 < land_fraction < 0.8.
- Balasus et al. (2023) find persistent coastal biases, also along lakes and large rivers, caused by sub-pixel albedo contrast. Removing surface_classification = 3 pixels, plus class 2 pixels with SWIR χ² > 20,000, largely removes them at a cost of 0.6% of the data.
- Schuit et al. (2023) list coastlines and cloud borders as common artefact locations.
- Vanselow et al. (2024) reject regions whose largest enhancements often occur over or next to water.
- GEE confirms that ocean-glint data enter only from November 2021, so coastal cells change their mix of observations over time.

GEE lacks land_fraction, so reproduce the filter from a static water mask. Compute the water fraction of each analysis cell within a 7×7 km footprint window from JRC Global Surface Water (occurrence ≥ 50%) or MOD44W. Exclude cells with water fraction > 0.2, which is equivalent to land_fraction < 0.8, and exclude ocean cells. Run sensitivity tests that also remove a 1-footprint (about 7 km) and a 2-footprint (about 14 km) buffer around:
- the coastline and estuaries, including São Marcos and São José bays and the Parnaíba delta;
- large reservoirs: Sobradinho, Itaparica, Xingó, Boa Esperança, Castanhão, Orós, Armando Ribeiro Gonçalves;
- the seasonal interdunal lagoons of the Lençóis Maranhenses.

The water fraction must be added as a matching covariate in the proximity test (topic 6).

**Parâmetros:** Exclude cells with 0 < land fraction < 0.8 (water fraction > 0.2), computed in a 7×7 km window from JRC GSW occurrence ≥ 50% (GEE asset JRC/GSW1_4/GlobalSurfaceWater) or MOD44W. Exclude pure-water and ocean cells to remove glint data (present only from Nov 2021). Sensitivity buffers: 7 km and 14 km from coast and large reservoirs. Level-2 alternative: drop surface_classification = 3, and class 2 with SWIR χ² > 20,000 (Balasus 2023). Diagnostic from this check: 78.0% of 7-km cells within 5 km of ANP pipelines lie within 50 km of the coast, against 18.2% of all Northeast cells. Within 20 km: 30.3% vs 5.9%. Within 10 km: 12.7% vs 2.7%. Median distance to coast: 28 km vs 114 km (n = 821 vs 31,712 cells).

**Notas:** I computed the coastal-share diagnostic for this report with a scratch script: (local working copy, not archived) It uses the four ANP pipeline shapefiles clipped to Northeast Brazil, IBGE-style state polygons simplified by 300 m, the Brazil outline as the coastline, and a 7-km point grid in the South America Albers projection. The numbers are approximate but the contrast is large. Pipeline-proximal cells are four times more likely than the region to lie within 50 km of the coast, so any coastal artefact appears as a 'pipeline effect'. Several top-20 peaks in Table 4 appear to lie on or next to water: rank 1 (−1.53, −45.91) on the Reentrâncias Maranhenses mangrove coast, rank 18 (−2.62, −44.47) on São Marcos Bay, rank 8 (−9.75, −41.68) on the Sobradinho shore near Sento Sé, and rank 19 (−8.79, −38.72) near the Itaparica reservoir around Rodelas. Confirm these with the JRC overlay before stating them in the paper.

**Referências:**

- Landgraf, J., Lorente, A., Borsdorff, T., et al. (2025). S5P MPC Methane Product Readme File, S5P-MPC-SRON-PRF-CH4 issue 2.9 (V02.09.01), sections 3.1 and 4.1. https://doi.org/10.5270/S5P-3lcdqiv
- Balasus, N., et al. (2023). A blended TROPOMI+GOSAT satellite data product for atmospheric methane using machine learning to correct retrieval biases. Atmos. Meas. Tech., 16, 3787–3807 (Appendix D: filtering coastal TROPOMI scenes). https://doi.org/10.5194/amt-16-3787-2023
- Schuit, B. J., et al. (2023). Automated detection and monitoring of methane super-emitters using satellite data. Atmos. Chem. Phys., 23, 9071–9098. https://doi.org/10.5194/acp-23-9071-2023
- Vanselow, S., et al. (2024). Automated detection of regions with persistently enhanced methane concentrations using Sentinel-5 Precursor satellite data. Atmos. Chem. Phys., 24, 10441–10473. https://doi.org/10.5194/acp-24-10441-2024
- Pekel, J.-F., Cottam, A., Gorelick, N., Belward, A. S. (2016). High-resolution mapping of global surface water and its long-term changes. Nature, 540, 418–422. https://doi.org/10.1038/nature20584
- Google Earth Engine Data Catalog: COPERNICUS/S5P/OFFL/L3_CH4 (statement that glint ocean observations were added from November 2021). https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S5P_OFFL_L3_CH4

## 6. 4b. Albedo, surface-spectral, aerosol, cloud-edge and smoke artefacts; what the bias-corrected band does and does not remove

State explicitly that the bias-corrected XCH4 removes only a global dependence on SWIR surface albedo. Lorente et al. (2021) fit a B-spline to the ratio of reference to TROPOMI XCH4 in SWIR-albedo bins of 0.05, referenced to albedo 0.2 through a 'small-area approximation'. The correction does not remove:
(a) surface-spectral artefacts, such as carbonate rocks, of up to 20–40 ppb; the third-order surface polynomial that fixes them exists only from v02.04.00 (Lorente et al. 2023), so GEE data before 17 Jul 2022 are uncorrected;
(b) biases from aerosol size and scattering, cirrus, coasts and striping (Balasus et al. 2023). Their largest residual corrections are associated with coarse aerosols, high SWIR albedo and across-track pixel index, and exceed 10 ppb over arid land and persistently cloudy regions. They also suggest the operational albedo correction pushes XCH4 too low over bright surfaces;
(c) albedo-shaped plumes: Lauvaux et al. (2022) state that albedo biases are 'not properly removed' and discarded candidates that correlate with SWIR albedo or match optical-imagery patterns.

Filters to apply in GEE:
- SWIR aerosol_optical_depth < 0.13 (Schuit 2023; Maasakkers 2022).
- Precision < 10 ppb.
- Surface albedo diagnostic, because L3 has no albedo band. Use MODIS MCD43A4 Band 7 (2105–2155 nm) NBAR as a SWIR-albedo proxy. Reject clusters where |r| between ΔXCH4 and the SWIR proxy exceeds 0.5 within the cluster dilated by 2 cells; apply the same rule to AOD. Sadavarte et al. (2021) rejected orbits with |R| > 0.5 against albedo or AOT. Zhang et al. (2020) reported R² = 0.00 with albedo and 0.07 with AOT as evidence of robustness; report the same statistics.
- de Gouw et al. (2020) used albedo > 0.05 as a sensitivity bound (qa already removes albedo < 0.02).
- Terrain: exclude cells whose DEM SD within the cell exceeds 150 m (Vanselow et al. 2024). qa already drops roughness above 80 m.
- Clouds: drop orbit images with less than 20% valid pixels in the background window (Schuit: scenes with more than 20% valid pixels). As a sensitivity test, erode the valid mask by one footprint around data gaps to suppress cloud-border artefacts.
- Smoke: mask cell-days with VIIRS 375 m active fires within about 50 km, or with a TROPOMI CO enhancement (COPERNICUS/S5P/OFFL/L3_CO), since smoke aerosols bias XCH4.
- Shape test: reject anomalies whose elongation follows surface features rather than the wind; about 30° acceptance in Lauvaux's labelling (see topic 8).

**Parâmetros:** AOD (SWIR) < 0.13. Precision < 10 ppb. Albedo-proxy and AOD correlation: reject candidates with |r| > 0.5 (cluster plus 2-cell dilation). Optional albedo > 0.05 sensitivity. Cell elevation SD ≤ 150 m. Orbit images: ≥ 20% valid pixels in the window. Cloud-edge erosion: 1 footprint (about 7 km). Fire mask: VIIRS 375 m detections within about 50 km on the same day, or TROPOMI CO anomaly above p90 of the local window. Full Schuit L2 filter set for reference (not all reproducible in GEE): QA ≥ 0.4, precision < 10 ppb, SWIR AOD < 0.13, NIR AOD < 0.30, SWIR albedo > 0.02, mixed albedo 2.4·NIR − 1.13·SWIR < 0.95, SWIR cloud fraction < 0.02.

**Notas:** The MODIS Band 7 proxy is my adaptation for GEE, not a published TROPOMI filter. MODIS Band 7 is at 2.1 µm, TROPOMI's SWIR albedo at 2.3 µm, so present the check as a correlation diagnostic rather than a threshold on absolute albedo. The 50 km smoke radius and the CO-p90 rule are also design choices and need sensitivity tests. Geological hypothesis worth testing: the Jandaíra Formation carbonates crop out widely in the onshore Potiguar Basin (Rio Grande do Norte and Ceará), where onshore O&G is concentrated. Lorente et al. (2023) showed that carbonate surfaces produced 20–40 ppb artefacts before v02.04.00. The corrected table in the author's notes gives Rio Grande do Norte the highest state mean (+6.74 ppb); the P1/P2 version split (topic 0) can test whether that mean is partly an artefact. Reference for the formation: Santos Filho, M.A.B. et al. (2015), Braz. J. Geol. 45, 23–34, doi:10.1590/23174889201500010002.

**Referências:**

- Lorente, A., et al. (2021). Methane retrieved from TROPOMI: improvement of the data product and validation of the first 2 years of measurements. Atmos. Meas. Tech., 14, 665–684. https://doi.org/10.5194/amt-14-665-2021
- Lorente, A., Borsdorff, T., Martinez-Velarte, M. C., Landgraf, J. (2023). Accounting for surface reflectance spectral features in TROPOMI methane retrievals. Atmos. Meas. Tech., 16, 1597–1608. https://doi.org/10.5194/amt-16-1597-2023
- Balasus, N., et al. (2023). A blended TROPOMI+GOSAT satellite data product for atmospheric methane using machine learning to correct retrieval biases. Atmos. Meas. Tech., 16, 3787–3807. https://doi.org/10.5194/amt-16-3787-2023
- Lauvaux, T., et al. (2022). Global assessment of oil and gas methane ultra-emitters. Science, 375, 557–561. https://doi.org/10.1126/science.abj4351
- Schuit, B. J., et al. (2023). Automated detection and monitoring of methane super-emitters using satellite data. Atmos. Chem. Phys., 23, 9071–9098. https://doi.org/10.5194/acp-23-9071-2023
- Maasakkers, J. D., Varon, D. J., Elfarsdóttir, A., McKeever, J., Jervis, D., Mahapatra, G., et al. (2022). Using satellites to uncover large methane emissions from landfills. Sci. Adv., 8, eabn9683. https://doi.org/10.1126/sciadv.abn9683
- Sadavarte, P., et al. (2021). Methane emissions from superemitting coal mines in Australia quantified using TROPOMI satellite observations. Environ. Sci. Technol., 55, 16573–16580. https://doi.org/10.1021/acs.est.1c03976
- Zhang, Y., et al. (2020). Quantifying methane emissions from the largest oil-producing basin in the United States from space. Sci. Adv., 6, eaaz5120. https://doi.org/10.1126/sciadv.aaz5120
- de Gouw, J. A., et al. (2020). Daily satellite observations of methane from oil and gas production regions in the United States. Sci. Rep., 10, 1379. https://doi.org/10.1038/s41598-020-57678-4
- Vanselow, S., et al. (2024). Automated detection of regions with persistently enhanced methane concentrations using Sentinel-5 Precursor satellite data. Atmos. Chem. Phys., 24, 10441–10473. https://doi.org/10.5194/acp-24-10441-2024
- Schaaf, C. B., Gao, F., Strahler, A. H., Lucht, W., Li, X., Tsang, T., et al. (2002). First operational BRDF, albedo nadir reflectance products from MODIS. Remote Sens. Environ., 83, 135–148. https://doi.org/10.1016/S0034-4257(02)00091-3
- Schroeder, W., Oliva, P., Giglio, L., Csiszar, I. A. (2014). The New VIIRS 375 m active fire detection data product: Algorithm description and initial assessment. Remote Sens. Environ., 143, 85–96. https://doi.org/10.1016/j.rse.2013.12.008

## 7. 5. Robust thresholds and false-discovery control under spatial autocorrelation

Replace the 'Z > 2 against global μ/σ' rule with robust local statistics plus a persistence rule and multiplicity control.

(a) Robust local score:
z_rob = (ΔXCH4 − median_W) / (1.4826·MAD_W)
Use z_rob ≥ 2.5 as the candidate threshold (Leys et al. 2013, 'moderately conservative') and ≥ 3.5 as the strict threshold (Iglewicz and Hoaglin's modified Z). Alternatively follow Vanselow et al. (2024): ΔXCH4 ≥ 2σ on monthly maps, with σ computed after excluding values above p95, and persistence F_enh ≥ 0.5 of valid months.

Literature anchors for single-overpass work:
- Lauvaux et al. (2022): 3σ_local, a mean enhancement above 25 ppb, and at least 5 pixels with qa>0.75.
- Barré et al. (2021): flag when |filtered departure| exceeds TROPOMI precision (1%, about 18 ppb).

(b) Cell-level inference: test each cell's series of monthly ΔXCH4 against zero, one-sided. Use a Wilcoxon signed-rank test, or a moving-block bootstrap with 3-month blocks to respect serial correlation. Then control the false discovery rate: Benjamini–Hochberg at q = 0.10, following Wilks (2016), who shows that α_FDR ≈ 2·α_global gives about the right global error rate under moderate to strong spatial correlation. Report Benjamini–Yekutieli as the dependence-robust, more conservative variant. If local Getis–Ord Gi* or LISA statistics are used, apply the FDR procedure of Caldas de Castro and Singer (2006).

(c) Clusters: require at least 2 contiguous cells and persistence in a 3×3 moving window (Vanselow). Drop the erosion/dilation step as a substitute for inference.

**Parâmetros:** Robust scale: 1.4826·MAD within the 2° window (or Rousseeuw–Croux Qn). Candidate: z_rob ≥ 2.5. Strict: modified Z ≥ 3.5. Vanselow alternative: Nσ = 2 with σ excluding values above p95, plus F_enh ≥ 0.5 and N_meas ≥ 24 months. Absolute floor: monthly ΔXCH4 ≥ 10 ppb, or ≥ 18 ppb (≈1% precision) for single-orbit screening. FDR: BH q = 0.10 (≈ 2 × 0.05, per Wilks 2016); BY q = 0.05 as sensitivity. Temporal block bootstrap: 3-month blocks, B = 9,999. Minimum cluster: 2 cells (≈100 km² at 7 km) with 3×3 persistence.

**Notas:** Key point for the manuscript. A one-sided Z > 2 tail holds 2.28% of a Gaussian field, i.e. about 750 of the 32,904 cells at 7 km would exceed it with no signal at all. The reported 18,000 km² is about 367 cells at 49 km² (1.1–1.2%), which is not an excess over that no-signal expectation. The code comment calling Z > 2 a '95% statistical confidence' threshold is therefore incorrect, because it is a per-cell descriptive cut with no multiplicity control. Pair the robust thresholds with the persistence and FDR steps, and state the expected null exceedance in the text. Rousseeuw and Croux's Qn is a more efficient robust scale than MAD when a window is heavily skewed by sources.

**Referências:**

- Leys, C., Ley, C., Klein, O., Bernard, P., Licata, L. (2013). Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median. J. Exp. Soc. Psychol., 49, 764–766. https://doi.org/10.1016/j.jesp.2013.03.013
- Rousseeuw, P. J., Croux, C. (1993). Alternatives to the median absolute deviation. J. Am. Stat. Assoc., 88, 1273–1283. https://doi.org/10.1080/01621459.1993.10476408
- Iglewicz, B., Hoaglin, D. C. (1993). How to Detect and Handle Outliers. ASQC Basic References in Quality Control, vol. 16. ASQC Quality Press, Milwaukee (no DOI).
- Benjamini, Y., Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. J. R. Stat. Soc. B, 57, 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Benjamini, Y., Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Ann. Stat., 29, 1165–1188. https://doi.org/10.1214/aos/1013699998
- Wilks, D. S. (2016). 'The stippling shows statistically significant grid points': How research results are routinely overstated and overinterpreted, and what to do about it. Bull. Am. Meteorol. Soc., 97, 2263–2273. https://doi.org/10.1175/BAMS-D-15-00267.1
- Ventura, V., Paciorek, C. J., Risbey, J. S. (2004). Controlling the proportion of falsely rejected hypotheses when conducting multiple tests with climatological data. J. Climate, 17, 4343–4356. https://doi.org/10.1175/3199.1
- Caldas de Castro, M., Singer, B. H. (2006). Controlling the false discovery rate: A new application to account for multiple and dependent tests in local statistics of spatial association. Geogr. Anal., 38, 180–208. https://doi.org/10.1111/j.0016-7363.2006.00682.x
- Ord, J. K., Getis, A. (1995). Local spatial autocorrelation statistics: Distributional issues and an application. Geogr. Anal., 27, 286–306. https://doi.org/10.1111/j.1538-4632.1995.tb00912.x
- Anselin, L. (1995). Local indicators of spatial association—LISA. Geogr. Anal., 27, 93–115. https://doi.org/10.1111/j.1538-4632.1995.tb00338.x
- Künsch, H. R. (1989). The jackknife and the bootstrap for general stationary observations. Ann. Stat., 17, 1217–1241. https://doi.org/10.1214/aos/1176347265
- Vanselow, S., et al. (2024). Automated detection of regions with persistently enhanced methane concentrations using Sentinel-5 Precursor satellite data. Atmos. Chem. Phys., 24, 10441–10473. https://doi.org/10.5194/acp-24-10441-2024
- Lauvaux, T., et al. (2022). Global assessment of oil and gas methane ultra-emitters. Science, 375, 557–561. https://doi.org/10.1126/science.abj4351
- Barré, J., et al. (2021). Systematic detection of local CH4 anomalies by combining satellite measurements with high-resolution forecasts. Atmos. Chem. Phys., 21, 5117–5136. https://doi.org/10.5194/acp-21-5117-2021

## 8. 6. Testing proximity effects (mean and 'variance signature') under spatial autocorrelation

The current comparison, the SD of 7-km cells of a 6-year mean inside a 5-km buffer (12.89 ppb) against the whole region (9.39 ppb), mixes spatial heterogeneity with emissions. It uses independent-sample logic on autocorrelated cells, and it is confounded by coastal location (topic 4a). Redesign as a matched, spatially aware test.

(1) Outcomes per cell:
(a) median of monthly ΔXCH4 (location);
(b) temporal robust dispersion, the MAD of monthly ΔXCH4, normalized by its expected noise σ_i/√n_i. The 'episodic emission' hypothesis concerns temporal variance, not the spatial SD of a multi-year mean.

(2) Groups: treated cells within 5 km of a pipeline; a 5–25 km exclusion ring; control cells more than 25 km from any O&G asset (pipelines, wells, stations, the 417 ANP assets).

(3) Matching covariates: distance to coast (log), water fraction, elevation and roughness, SWIR-albedo proxy, latitude, number of valid months, MapBiomas class (mangrove, wetland/flooded field, pasture, urban), VIIRS fire frequency, and optionally cattle density. Use 1:k nearest-neighbour matching (k = 3–5) with a caliper of 0.2 SD of the logit propensity score, and require all standardized mean differences below 0.1.

(4) Inference that respects autocorrelation, in one or more of these forms:
(a) spatial block bootstrap: resample 1°×1° tiles, or blocks at least as large as the empirical semivariogram range of ΔXCH4;
(b) placement null: randomly rotate (0–360°) and translate the pipeline network (torus or translation test), at least 999 placements, keeping only placements that stay on land and preserve the distance-to-coast distribution;
(c) Moran spectral randomization;
(d) effective sample size correction (Clifford et al. 1989; Dutilleul 1993) for correlations with distance.

For dispersion, compute Brown–Forsythe (median-centred Levene) or Fligner–Killeen statistics, but take their p-values from the spatial null, not from iid tables. Test the claim that 4 of 21 clusters (19%) intersect pipelines against the same placement null. Report effect sizes with bootstrap confidence intervals.

**Parâmetros:** Treated: ≤ 5 km. Exclusion ring: 5–25 km. Control: > 25 km from any asset. Matching: 1:3–5 nearest neighbour, caliper 0.2 SD of logit PS, balance SMD < 0.1. Block bootstrap: tiles of 1°×1° or ≥ variogram range (fit an exponential or spherical variogram to ΔXCH4; expect tens to a few hundred km), B = 9,999. Placement null: ≥ 999 (ideally 9,999) rotations or translations, constrained to land and to the distance-to-coast distribution; p = (1 + #{T* ≥ T_obs})/(1 + B). Dispersion statistic: Brown–Forsythe or Fligner–Killeen on the MAD of monthly anomalies, normalized by σ_i/√n_i.

**Notas:** Why this matters here: 78% of pipeline-proximal cells lie within 50 km of the coast, against 18% of the region (topic 4a). Coastal cells combine land–water boundary low biases, ocean glint present only from November 2021, and mangrove, estuary and urban sources. All of these raise spatial dispersion. The 37% higher SD may therefore be a coastal signature, and this should be tested before it is interpreted. The 5–25 km exclusion ring and the > 25 km control distance are my choices and need sensitivity checks (for example 15/50 km). For the placement null, rotating and translating the real network keeps its internal geometry, including network length and clustering in capitals, which random points would not.

**Referências:**

- Künsch, H. R. (1989). The jackknife and the bootstrap for general stationary observations. Ann. Stat., 17, 1217–1241. https://doi.org/10.1214/aos/1176347265
- Politis, D. N., Romano, J. P. (1994). The stationary bootstrap. J. Am. Stat. Assoc., 89, 1303–1313. https://doi.org/10.1080/01621459.1994.10476870
- Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer Series in Statistics, Springer, New York. https://doi.org/10.1007/978-1-4757-3803-2
- Harms, K. E., Condit, R., Hubbell, S. P., Foster, R. B. (2001). Habitat associations of trees and shrubs in a 50-ha neotropical forest plot. J. Ecol., 89, 947–959 (torus-translation test). https://doi.org/10.1111/j.1365-2745.2001.00615.x
- Wagner, H. H., Dray, S. (2015). Generating spatially constrained null models for irregularly spaced data using Moran spectral randomization methods. Methods Ecol. Evol., 6, 1169–1178. https://doi.org/10.1111/2041-210X.12407
- Viladomat, J., Mazumder, R., McInturff, A., McCauley, D. J., Hastie, T. (2014). Assessing the significance of global and local correlations under spatial autocorrelation: A nonparametric approach. Biometrics, 70, 409–418. https://doi.org/10.1111/biom.12139
- Clifford, P., Richardson, S., Hémon, D. (1989). Assessing the significance of the correlation between two spatial processes. Biometrics, 45, 123–134. https://doi.org/10.2307/2532039
- Dutilleul, P., Clifford, P., Richardson, S., Hémon, D. (1993). Modifying the t test for assessing the correlation between two spatial processes. Biometrics, 49, 305–314. https://doi.org/10.2307/2532625
- Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera-Arroita, G., et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography, 40, 913–929. https://doi.org/10.1111/ecog.02881
- Dormann, C. F., McPherson, J. M., Araújo, M. B., Bivand, R., Bolliger, J., Carl, G., et al. (2007). Methods to account for spatial autocorrelation in the analysis of species distributional data: a review. Ecography, 30, 609–628. https://doi.org/10.1111/j.2007.0906-7590.05171.x
- Stuart, E. A. (2010). Matching methods for causal inference: A review and a look forward. Stat. Sci., 25, 1–21. https://doi.org/10.1214/09-STS313
- Ho, D. E., Imai, K., King, G., Stuart, E. A. (2007). Matching as nonparametric preprocessing for reducing model dependence in parametric causal inference. Polit. Anal., 15, 199–236. https://doi.org/10.1093/pan/mpl013
- Austin, P. C. (2011). Optimal caliper widths for propensity-score matching when estimating differences in means and differences in proportions in observational studies. Pharm. Stat., 10, 150–161. https://doi.org/10.1002/pst.433
- Brown, M. B., Forsythe, A. B. (1974). Robust tests for the equality of variances. J. Am. Stat. Assoc., 69, 364–367. https://doi.org/10.1080/01621459.1974.10482955
- Conover, W. J., Johnson, M. E., Johnson, M. M. (1981). A comparative study of tests for homogeneity of variances, with applications to the outer continental shelf bidding data. Technometrics, 23, 351–361. https://doi.org/10.1080/00401706.1981.10487680
- Legendre, P. (1993). Spatial autocorrelation: trouble or new paradigm? Ecology, 74, 1659–1673. https://doi.org/10.2307/1939924

## 9. 7. Natural and non-O&G sources in Northeast Brazil to screen and attribute

Build an attribution overlay for every candidate, and check the seasonal signature, before linking any anomaly to O&G.

(1) Mangroves and estuaries. The Pará–Maranhão Amazon Macrotidal Mangrove Coast holds about 7,600 km² of continuous mangrove, among the largest in the world (Souza Filho 2005). Amazon-coast mangrove soils emit substantial CH4, with 1,439–3,312 µg C m⁻² h⁻¹ near the Amazon mouth, favoured by low salinity (Bernardino et al. 2024). CH4 from mangroves is significant globally (Rosentreter et al. 2018), and aquatic sources are highly variable (Rosentreter et al. 2021). The rank-1 peak sits on this coast (Reentrâncias Maranhenses).

(2) Floodplains and wetlands. The Baixada Maranhense pre-Amazonian floodplain lakes (Ibañez 2000) and other seasonally flooded fields. Tropical floodplain CH4 is well established (Melack et al. 2004). Eastern Amazonia shows large and rising emissions, with the largest positive trend in the December–March wet season (Wilson et al. 2021), and regional variation in the airborne budget (Basso et al. 2021). Vanselow et al. (2024) found that 46.5% of global persistent enhancements were not assigned to any source type, mostly in wetland regions where WetCHARTs gives less than 50 kt/yr (Bloom et al. 2017). Inventory absence therefore does not exclude a wetland origin.

(3) Reservoirs. Sobradinho, Brazil's largest reservoir by area, emits about 1.33 Mt CO2-eq/yr from all carbon, with diffusion 96% of the flux and littoral zones 40% (Rodriguez and Casper 2018). The global synthesis of Deemer et al. (2016) supports reservoir CH4. Reservoir shores are also land–water artefact zones (topic 4a). Ranks 8 and 19 appear to lie on the Sobradinho and Itaparica shores.

(4) Biomass burning. The MATOPIBA Cerrado and the eastern Amazon frontier in Maranhão/Pará have a dry-season fire regime (Silva et al. 2021). Use GFED4.1s (van der Werf et al. 2017), VIIRS 375 m active fires (Schroeder et al. 2014) and emission factors from Andreae (2019). Smoke also biases the retrieval (topic 4b). Ranks 5, 10 and 11 (about 4.4–4.9°S, 46.2–48.0°W) fall in this frontier.

(5) Livestock, landfills and rice. Livestock dominates the anthropogenic budget: 65% of South America's posterior anthropogenic total in Hancock et al. (2025). Landfills near state capitals can be detected by TROPOMI (Maasakkers et al. 2022). For O&G priors use GFEI v2 (Scarpelli et al. 2022) and ANP assets; for other sectors use EDGAR and the national SEEG estimates.

The discriminating test is phase. Wetland, mangrove and reservoir signals should rise in the flood or wet season (Feb–Jul, depending on the sub-region). Fire signals should appear in Aug–Nov together with CO. O&G signals should be aseasonal and should not switch on or off at processor-version dates.

**Parâmetros:** Overlay layers, all in or importable to GEE: Global Mangrove Watch v3 (Bunting et al. 2022); MapBiomas land use and land cover (Souza et al. 2020), including mangrove, flooded field and wetland classes; JRC GSW occurrence, seasonality and max_extent (Pekel et al. 2016); WetCHARTs (Bloom et al. 2017); GFED4.1s and VIIRS/MODIS FIRMS; TROPOMI CO L3; GFEI v2 and EDGAR; SEEG. Attribution rule after Vanselow et al. (2024): sum emissions per sector inside the cluster mask dilated by one cell, and label the cluster 'unknown' if the total is below 50 kt/yr. Seasonal-phase test: compare the monthly climatology of ΔXCH4 per cluster with flood (JRC monthly), fire (FIRMS counts) and CO anomaly series.

**Notas:** Rank coordinates come from Table4_Critical_Targets_Ranking.csv. The placements on water bodies and the frontier are my visual and geographic reading, to be confirmed with the overlays; do not state them in the paper before that. Vanselow et al. (2024) detected persistent enhancements over the Pantanal and Argentine wetlands, which shows that South American wetlands produce persistent TROPOMI anomalies of the kind the manuscript flags. No peer-reviewed CH4 flux study specific to the Baixada Maranhense was found; cite it as a floodplain system (Ibañez 2000) and use the WetCHARTs and JRC seasonality evidence. The wet-season timing (Feb–Jul) varies within Northeast Brazil: roughly Feb–May in the north (MA, PI, CE), Apr–Jul on the east coast, and Nov–Mar in western Bahia. Define the phase test per sub-region.

**Referências:**

- Souza Filho, P. W. M. (2005). Costa de manguezais de macromaré da Amazônia: cenários morfológicos, mapeamento e quantificação de áreas usando dados de sensores remotos. Rev. Bras. Geofís., 23, 427–435. https://doi.org/10.1590/S0102-261X2005000400006
- Bernardino, A. F., Queiroz, H. M., Nobrega, G. N., Coppo, G. C., Sanders, C. J., Silva, A. E., et al. (2024). Soil greenhouse gas fluxes partially reduce the net gains in carbon sequestration in mangroves of the Brazilian Amazon. Environ. Res., 263, 120102. https://doi.org/10.1016/j.envres.2024.120102
- Rosentreter, J. A., Maher, D. T., Erler, D. V., Murray, R. H., Eyre, B. D. (2018). Methane emissions partially offset 'blue carbon' burial in mangroves. Sci. Adv., 4, eaao4985. https://doi.org/10.1126/sciadv.aao4985
- Rosentreter, J. A., Borges, A. V., Deemer, B. R., Holgerson, M. A., Liu, S., Song, C., et al. (2021). Half of global methane emissions come from highly variable aquatic ecosystem sources. Nat. Geosci., 14, 225–230. https://doi.org/10.1038/s41561-021-00715-2
- Ibañez, M. S. R. (2000). Limnological characteristics of three aquatic systems of the pre-amazonian floodplain, Baixada Maranhense (Maranhão, Brazil). Aquat. Ecosyst. Health Manag., 3, 521–531. https://doi.org/10.1016/S1463-4988(00)00048-8
- Melack, J. M., Hess, L. L., Gastil, M., Forsberg, B. R., Hamilton, S. K., Lima, I. B. T., Novo, E. M. L. M. (2004). Regionalization of methane emissions in the Amazon Basin with microwave remote sensing. Glob. Change Biol., 10, 530–544. https://doi.org/10.1111/j.1365-2486.2004.00763.x
- Wilson, C., Chipperfield, M. P., Gloor, M., Parker, R. J., Boesch, H., McNorton, J., et al. (2021). Large and increasing methane emissions from eastern Amazonia derived from satellite data, 2010–2018. Atmos. Chem. Phys., 21, 10643–10669. https://doi.org/10.5194/acp-21-10643-2021
- Basso, L. S., Marani, L., Gatti, L. V., Miller, J. B., Gloor, M., Melack, J., et al. (2021). Amazon methane budget derived from multi-year airborne observations highlights regional variations in emissions. Commun. Earth Environ., 2, 246. https://doi.org/10.1038/s43247-021-00314-4
- Bloom, A. A., Bowman, K. W., Lee, M., Turner, A. J., Schroeder, R., Worden, J. R., et al. (2017). A global wetland methane emissions and uncertainty dataset for atmospheric chemical transport models (WetCHARTs version 1.0). Geosci. Model Dev., 10, 2141–2156. https://doi.org/10.5194/gmd-10-2141-2017
- Rodriguez, M., Casper, P. (2018). Greenhouse gas emissions from a semi-arid tropical reservoir in northeastern Brazil. Reg. Environ. Change, 18, 1901–1912. https://doi.org/10.1007/s10113-018-1289-7
- Deemer, B. R., Harrison, J. A., Li, S., Beaulieu, J. J., DelSontro, T., Barros, N., et al. (2016). Greenhouse gas emissions from reservoir water surfaces: A new global synthesis. BioScience, 66, 949–964. https://doi.org/10.1093/biosci/biw117
- Silva, P. S., Nogueira, J., Rodrigues, J. A., Santos, F. L. M., Pereira, J. M. C., DaCamara, C. C., et al. (2021). Putting fire on the map of Brazilian savanna ecoregions. J. Environ. Manage., 296, 113098. https://doi.org/10.1016/j.jenvman.2021.113098
- van der Werf, G. R., Randerson, J. T., Giglio, L., van Leeuwen, T. T., Chen, Y., Rogers, B. M., et al. (2017). Global fire emissions estimates during 1997–2016. Earth Syst. Sci. Data, 9, 697–720. https://doi.org/10.5194/essd-9-697-2017
- Andreae, M. O. (2019). Emission of trace gases and aerosols from biomass burning – an updated assessment. Atmos. Chem. Phys., 19, 8523–8546. https://doi.org/10.5194/acp-19-8523-2019
- Hancock, S. E., et al. (2025). Satellite quantification of methane emissions from South American countries: a high-resolution inversion of TROPOMI and GOSAT observations. Atmos. Chem. Phys., 25, 797–817. https://doi.org/10.5194/acp-25-797-2025
- Scarpelli, T. R., Jacob, D. J., Grossman, S., Lu, X., Qu, Z., Sulprizio, M. P., et al. (2022). Updated Global Fuel Exploitation Inventory (GFEI) for methane emissions from the oil, gas, and coal sectors: evaluation with inversions of atmospheric methane observations. Atmos. Chem. Phys., 22, 3235–3249. https://doi.org/10.5194/acp-22-3235-2022
- Maasakkers, J. D., et al. (2022). Using satellites to uncover large methane emissions from landfills. Sci. Adv., 8, eabn9683. https://doi.org/10.1126/sciadv.abn9683
- Bunting, P., Rosenqvist, A., Hilarides, L., Lucas, R. M., Thomas, N., Tadono, T., et al. (2022). Global mangrove extent change 1996–2020: Global Mangrove Watch version 3.0. Remote Sens., 14, 3657. https://doi.org/10.3390/rs14153657
- Souza, C. M., Shimbo, J. Z., Rosa, M. R., Parente, L. L., Alencar, A. A., Rudorff, B. F. T., et al. (2020). Reconstructing three decades of land use and land cover changes in Brazilian biomes with Landsat archive and Earth Engine. Remote Sens., 12, 2735. https://doi.org/10.3390/rs12172735
- Vanselow, S., et al. (2024). Automated detection of regions with persistently enhanced methane concentrations using Sentinel-5 Precursor satellite data. Atmos. Chem. Phys., 24, 10441–10473. https://doi.org/10.5194/acp-24-10441-2024

## 10. 8. Wind-based methods with ERA5 (future work) and external plume cross-checks

Move beyond concentration maps by using winds.

(1) Wind-rotated oversampling to locate sources (Maasakkers et al. 2022; Sadavarte et al. 2021, SI S3; origin in Fioletov et al. 2015):
- Rotate each day's L2 pixels around a candidate source so that the local wind points in one direction. Use the ERA5 hourly 10-m wind closest to the overpass; Sadavarte used GEOS-FP 10 m.
- Oversample at 0.01° with footprint weighting (Sun et al. 2018).
- Scan a 13×13 grid of candidate sources at 0.05° spacing, then refine at 0.01°.
- Keep the location that maximizes the downwind minus upwind enhancement.

(2) Rotated-frame mass balance (Schneising et al. 2020): daily 0.05° grids rotated to the pressure-weighted boundary-layer wind from ERA5, a fixed upwind background region, wind speed 2–10 m/s, at least 25% coverage of the background, plume and hotspot regions, at least 60% coverage per cross-section, and exclusion of boxes below the 10th percentile within 700 km (residual clouds).

(3) Divergence method (Beirle et al. 2019 for NOx; Liu et al. 2021, 2024 and Veefkind et al. 2023 for CH4):
- Liu et al. (2021): 0.25° grid; background = mean of the lowest 10% within ±5 cells.
- Liu et al. (2024): 0.2° grid; ±3 cells; at least 10 retrievals per daily background; ERA5 wind at half the boundary-layer height (boundary layer fixed at 500 m); wind 0–10 m/s; more than 7 valid days; detection threshold about 3 kg km⁻² h⁻¹; non-divergent wind correction for orography; temporal filter keeping sources where more than 50% of days exceed background.
- Veefkind et al. (2023): daily 0.05° grid; ERA5 winds averaged from the surface to 500 m, 14 levels between 600 and 1000 hPa; background from a regression on surface pressure using the 25th percentile; about 10 km effective resolution; about 25% uncertainty.

(4) For single-overpass plumes: cross-sectional flux or IME with U_eff ≈ 1.05·U_BL (Varon et al. 2018; Sadavarte et al. 2021), excluding U < 2 m/s.

Cross-check every candidate against the SRON weekly TROPOMI plume list (Schuit et al. 2023) and against high-resolution tasking (Sentinel-2/Landsat, PRISMA/EnMAP/EMIT), using a tiered approach (Pandey et al. 2023; Jacob et al. 2022). A candidate with no plume detection in 2021–2025 is consistent with sources below TROPOMI's single-overpass threshold, not with super-emitters. Lauvaux et al. (2022) put the detection limit at about 25 t/h.

These methods need Level-2 pixels with geolocation and per-pixel time, which the GEE L3 lacks, so run them in Python on L2 files (RPRO/OFFL or the blended product) with ERA5 from the Copernicus Climate Data Store.

**Parâmetros:** Overpass timing: 13:30 local solar time is about 15:50–16:40 UTC across 35–48°W, so use ERA5 hours 16 and 17 UTC, time-interpolated. Rotation: ERA5 10-m u/v at the nearest hour; oversampling 0.01°; candidate grid 13×13 at 0.05° then 0.01°; wind 2–10 m/s. Rotated mass balance: pressure-weighted boundary-layer wind; fill ≥ 25% per region and ≥ 60% per cross-section; wind 2–10 m/s. Divergence: grid 0.05°–0.25°; background = mean of the lowest 10% within ±3–5 cells, or a surface-pressure regression at the 25th percentile; wind at about 250 m (half of a 500 m boundary layer) or the 0–500 m mean; ≥ 10 retrievals per daily background; > 7 valid days; non-divergent wind correction over the Borborema and Araripe plateaus; expected detection limit about 3 kg km⁻² h⁻¹. Cross-sectional flux: U_eff = (1.05 ± 0.17)·U_BL, discard U < 2 m/s. Plume cross-check source: SRON weekly TROPOMI plume list, 2021 onward.

**Notas:** Veefkind et al. (2023) used ERA5 at 17–19 UTC for the Permian. For Northeast Brazil the matching hours are about 16–17 UTC; this is my calculation from 13:30 local solar time and longitude. Liu et al. (2024) warn that the temporal filter used to reject artefacts can miss sources that release large amounts only occasionally, which may underestimate O&G emissions by about 40–70%. That caveat bears directly on the manuscript's 'episodic emission' hypothesis. In GEE, ERA5-Land hourly provides 10-m u/v for exploratory checks. For boundary-layer or pressure-level winds use ERA5 single and pressure levels from the Copernicus Climate Data Store.

**Referências:**

- Maasakkers, J. D., Varon, D. J., Elfarsdóttir, A., McKeever, J., Jervis, D., Mahapatra, G., et al. (2022). Using satellites to uncover large methane emissions from landfills. Sci. Adv., 8, eabn9683. https://doi.org/10.1126/sciadv.abn9683
- Sadavarte, P., et al. (2021). Methane emissions from superemitting coal mines in Australia quantified using TROPOMI satellite observations. Environ. Sci. Technol., 55, 16573–16580. https://doi.org/10.1021/acs.est.1c03976
- Fioletov, V. E., McLinden, C. A., Krotkov, N., Li, C. (2015). Lifetimes and emissions of SO2 from point sources estimated from OMI. Geophys. Res. Lett., 42, 1969–1976. https://doi.org/10.1002/2015GL063148
- Sun, K., Zhu, L., Cady-Pereira, K., Chan Miller, C., Chance, K., Clarisse, L., et al. (2018). A physics-based approach to oversample multi-satellite, multispecies observations to a common grid. Atmos. Meas. Tech., 11, 6679–6701. https://doi.org/10.5194/amt-11-6679-2018
- Schneising, O., Buchwitz, M., Reuter, M., Vanselow, S., Bovensmann, H., Burrows, J. P. (2020). Remote sensing of methane leakage from natural gas and petroleum systems revisited. Atmos. Chem. Phys., 20, 9169–9182. https://doi.org/10.5194/acp-20-9169-2020
- Beirle, S., Borger, C., Dörner, S., Li, A., Hu, Z., Liu, F., Wang, Y., Wagner, T. (2019). Pinpointing nitrogen oxide emissions from space. Sci. Adv., 5, eaax9800. https://doi.org/10.1126/sciadv.aax9800
- Liu, M., van der A, R., van Weele, M., Eskes, H., Lu, X., Veefkind, P., et al. (2021). A new divergence method to quantify methane emissions using observations of Sentinel-5P TROPOMI. Geophys. Res. Lett., 48, e2021GL094151. https://doi.org/10.1029/2021GL094151
- Liu, M., van der A, R., van Weele, M., Bryan, L., Eskes, H., Veefkind, P., et al. (2024). Current potential of CH4 emission estimates using TROPOMI in the Middle East. Atmos. Meas. Tech., 17, 5261–5277. https://doi.org/10.5194/amt-17-5261-2024
- Veefkind, J. P., Serrano-Calvo, R., de Gouw, J., Dix, B., Schneising, O., Buchwitz, M., et al. (2023). Widespread frequent methane emissions from the oil and gas industry in the Permian Basin. J. Geophys. Res. Atmos., 128, e2022JD037479. https://doi.org/10.1029/2022JD037479
- Varon, D. J., Jacob, D. J., McKeever, J., Jervis, D., Durak, B. O. A., Xia, Y., Huang, Y. (2018). Quantifying methane point sources from fine-scale satellite observations of atmospheric methane plumes. Atmos. Meas. Tech., 11, 5673–5686. https://doi.org/10.5194/amt-11-5673-2018
- Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz-Sabater, J., et al. (2020). The ERA5 global reanalysis. Q. J. R. Meteorol. Soc., 146, 1999–2049. https://doi.org/10.1002/qj.3803
- Schuit, B. J., et al. (2023). Automated detection and monitoring of methane super-emitters using satellite data. Atmos. Chem. Phys., 23, 9071–9098 (weekly plume list at earth.sron.nl/methane-emissions). https://doi.org/10.5194/acp-23-9071-2023
- Pandey, S., van Nistelrooij, M., Maasakkers, J. D., Sutar, P., Houweling, S., Varon, D. J., et al. (2023). Daily detection and quantification of methane leaks using Sentinel-3: a tiered satellite observation approach with Sentinel-2 and Sentinel-5p. Remote Sens. Environ., 296, 113716. https://doi.org/10.1016/j.rse.2023.113716
- Lauvaux, T., et al. (2022). Global assessment of oil and gas methane ultra-emitters. Science, 375, 557–561. https://doi.org/10.1126/science.abj4351

