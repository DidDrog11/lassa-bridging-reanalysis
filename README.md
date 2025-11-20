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

## Project Checklist: Lassa Spillover Reanalysis

This checklist maps the required data and analysis steps to the project workflow.

### Stage 1: Refined Host Ecological Modeling ($\boldsymbol{D_M}$)

| Status | Component | Description & Key R Functions |
| :---: | :--- | :--- |
| $\square$ | **A. Data Aggregation & Cleaning** | Clean and format all occurrence/non-detection data for the target species, including "true absence"" points. |
| $\square$ | **B. Site Master List & Covariates** | Create master site list and extract environmental predictors ($\boldsymbol{occ.covs}$) for all unique locations. |
| $\square$ | **C. IMSOM Data Preparation** | Structure the data into the necessary input format for `intMsPGOcc()`: **Source 1 (Replicated Surveys)** and **Source 2 (Aggregated Data)** arrays, including linking indices (`sites`, `species`). |
| $\square$ | **D. Benchmark JSDM (Model A)** | Run the **Hierarchical JSDM with Detection** on **Source 1** to establish a high-confidence baseline for parameter estimates. |
| $\square$ | **E. Full Coverage JSDM (Model B)** | Run the **Integrated JSDM** (`intMsPGOcc()`) on **Sources 1 & 2 combined** to generate the full-coverage $\boldsymbol{D_M}$ layer. Use `n.omp.threads` for parallelisation. |
| $\square$ | **F. Bias Assessment** | Quantify the consistency of $\boldsymbol{\psi}$ estimates and covariance ($\boldsymbol{\rho}_{ij}$) between Model A and Model B. |
| $\square$ | **G. Final $\boldsymbol{D_M}$ Layer** | Generate the spatial raster map of predicted *M. natalensis* suitability based on the full integrated Model B. |

### Stage 2: Epidemiological Integration and Incidence

| Status | Component | Description & Data Source |
| :---: | :--- | :--- |
| $\square$ | **H. Pathogen Layer ($\boldsymbol{D_L}$)** | Reuse the original study's $\boldsymbol{D_L}$ predictions (LASV in rodents). |
| $\square$ | **I. Composite Risk Layer ($\boldsymbol{D_X}$)** | Calculate the new combined risk: $\boldsymbol{D_X} = \boldsymbol{D_M} \times \boldsymbol{D_L}$. |
| $\square$ | **J. Seroreversion Parameter ($\boldsymbol{\lambda}$)** | Search and incorporate the latest literature to define a more precise range or single value for $\boldsymbol{\lambda}_{new}$. |
| $\square$ | **K. Seroprevalence Regression** | Re-run the Quasi-binomial Regression of human seroprevalence ($\boldsymbol{F}$) against the new $\boldsymbol{D_X}$ layer. |
| $\square$ | **L. Incidence Calculation** | Apply the steady-state SIR framework (Equation 6) using $\boldsymbol{D_X}$ and $\boldsymbol{\lambda}_{new}$ to calculate the predicted annual human infection rate ($\boldsymbol{F S^*}$). |

### Stage 3: Validation and Final Output

| Status | Component | Description & Expected Finding |
| :---: | :--- | :--- |
| $\square$ | **M. Reported Case Data (H)** | **CRITICAL NEW DATA:** Gather and clean regional/national Lassa Fever case reports. |
| $\square$ | **N. External Validation** | Compare the predicted infections ($\boldsymbol{F S^*}$) to the reported case data ($\boldsymbol{H}$). |
| $\square$ | **O. Underreporting Hotspot Mapping** | Map the residuals (Predicted Infections - Reported Cases) to highlight areas of potential underreporting. |
| $\square$ | **P. Final Results** | Summarize the revised range of annual LASV human infections and report on the ecological findings (species interaction effects). |