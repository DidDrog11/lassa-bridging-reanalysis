# Project: Lassa Virus Spillover Risk Reassessment (IMSOM-SIR)

### Project Rationale

This project undertakes a rigorous re-analysis and extension of the Lassa virus (LASV) spillover model published in Bridging the gap: Using reservoir ecology and human serosurveys to estimate Lassa virus spillover in West Africa (Basinski et al. 2021).

The primary objective is to **test three key hypotheses** regarding the drivers and burden of LASV spillover:

1.  **Ecological Improvement Test:** To assess whether incorporating **interspecies competition** and co-occurrence dynamics using a **Joint Species Distribution Model (JSDM)** (Doser, Finley and Banerjee. 2023) significantly improves the accuracy of the primary reservoir's (*Mastomys natalensis*) distribution layer ($\boldsymbol{D_M}$), compared to the original single-species model.
2.  **Epidemiological Uncertainty Reduction:** To incorporate the most current literature on LASV antibody waning ($\boldsymbol{\lambda}$) into the steady-state SIR framework to reduce the current $\sim 5$-fold uncertainty in annual human infection estimates.
3.  **Case Reporting Bias Detection:** To use the model's prediction of annual incidence ($\boldsymbol{F S^*}$) as an ecological baseline against which to validate against official reported case data. This comparison will spatially identify regions where underdiagnosis or underreporting is suspected to be substantial.

The analysis utilises **R**, with packages managed via `{renv}`, and employs a two-stage hierarchical modelling approach: an Integrated Multi-Species Occupancy Model (IMSOM) for the ecological stage, followed by calibration using a Quasi-binomial GLM and conversion via an SIR model.

References

Basinski, A. J., Fichet-Calvet, E., Sjodin, A. R., et al. 2021. Bridging the Gap: Using Reservoir Ecology and Human Serosurveys to Estimate Lassa Virus Spillover in West Africa. PLOS Computational Biology, 17(3): e1008811. DOI: [10.1371/journal.pcbi.1008811](https://doi.org/10.1371/journal.pcbi.1008811)

Doser, Jeffrey W, Andrew O Finley, and Sudipto Banerjee. 2023. Joint Species Distribution Models with Imperfect Detection for High-Dimensional Spatial Data. Ecology 104(9): e4137. DOI: [10.1002/ecy.4137](https://doi.org/10.1002/ecy.4137)

***

## Project Setup and Reproducibility

This project uses a provided template (VERERNA Consortium)[https://github.com/viralemergence/r-reproducible-repo] and relies on `{renv}` to manage package dependencies.

### 🚀 Getting Started

1.  **Clone Repository:**
    ```bash
    git clone [https://github.com/DidDrog11/lassa-bridging-reanalysis.git](https://github.com/DidDrog11/lassa-bridging-reanalysis.git)
    ```
2.  **Initialize Project:** Open `scripts/project_startup.R` and follow the instructions to set up the environment and restore packages using `renv::restore()`.
3.  **Code Workflow:** All analytical code will be executed sequentially from the `main_workflow.R` script, calling modular functions stored in the `/R` folder.

***

### Environmental Data

| Predictor |	Description |	Original | Source |	Updated Source / R Acquisition Strategy |
| :---: | :--- | :--- | :--- | :--- |
| Tmu |	Mean daily Temperature (C) |	MODIS |	Use Pre-Averaged Product (Mean Annual Temperature) / WorldClim via {geodata}. |
| Pmu |	Mean daily Precipitation (mm/day)	| CHIRPS | Use Pre-Averaged Product (Annual Precipitation) / WorldClim via {geodata}. |
| Nmu	| Mean NDVI	| MODIS |	MODIS products (e.g., MOD13A2) via {MODIStools}. |
| Pcv, Ncv |	Coefficient of Variation (CV) for Precip/NDVI |	Calculated |	Calculate from monthly rasters (2001–2025). |
| Pmin, Nmin |	Average Minimum Precip/NDVI per year |	Calculated |	Calculate from monthly rasters (2001–2025). |
| Pmax, Nmax |	Average Maximum Precip/NDVI per year |	Calculated |	Calculate from monthly rasters (2001–2025). |
| Pc, Nc |	Colwell's Constancy Index (Precip/NDVI)	| Calculated |	Calculate using {dismo} or custom R function on monthly data. |
| Pm, Nm |	Colwell's Contingency Index (Precip/NDVI)	Calculated |	Calculate using {dismo} or custom R function on monthly data. |
| Pdur |	Duration of low Precip (≤1 mm/day) |	Calculated |	Calculate from monthly rasters (2001–2025). |
| Ndur |	Duration of low NDVI (≤0.5) |	Calculated |	Calculate from monthly rasters (2001–2025). |
| 11 Land Cover Layers |	Local density (fraction) of stable land cover types in surrounding 0.15∘×0.15∘ area.|	MODIS Land Cover (2001)	ESA WorldCover (2020) or recent MODIS Land Cover (e.g., MCD12Q1). Requires calculating long-term stability (2001–2025).|
| Elev |	Elevation (m) |	USGS DEM |	SRTM or ASTER GDEM via {geodata}. |
| Pop |	Human Population Density (humans/km2) |	WorldPop 2020	| WorldPop 2025 Projection via WorldPop FTP/API. |
