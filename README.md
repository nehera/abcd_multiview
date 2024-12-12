# abcd_multiview

This repository supports the manuscript titled *Bayesian Integrative Mixed Modeling Framework for Analysis of the Adolescent Brain and Cognitive Development Study*. The analysis focuses on multiview data from the [ABCD Study®](https://abcdstudy.org/) to integrate heterogeneous data types, account for nested hierarchical structures, and predict behavioral outcomes using the BIPmixed framework. This framework extends the Bayesian Integrative Analysis and Prediction (BIP) model to incorporate 2-level random effects, improving predictive performance and interpretability for data with hierarchical structures.

## Authors

Aidan Neher, Apostolos Stamenos, Mark Fiecas, Sandra Safo, Thierry Chekouo @ Biostatistics and Health Data Science, University of Minnesota

## Repository Structure

The repository includes the following:

### Top-Level Files:

- **00_abcd_eda.docx**: Initial exploratory data analysis (EDA) document.
- **00_abcd_eda.qmd**: Quarto file for the EDA.
- **01_load_covars_and_outcomes.R**: Script for loading covariates and outcome data.
- **01_load_ELA.R**: Script for processing the Early Life Adversity (ELA) view.
- **01_load_imaging_data.R**: Script for processing imaging views (structural and functional MRI).
- **02_combine_data_define_sample.R**: Combines and merges the views to define the analysis sample.
- **03_generate_scree_plots.R**: Generates scree plots for determining latent component hyperparameter r.
- **data_analysis.R**: Main script for performing BIPmixed and BIP data analysis.
- **data_analysis_summary.R**: Summarizes the results of the BIPmixed and BIP analysis.
- **simulation_study.R**: Script for running simulation studies.
- **simulation_study_summary.R**: Summarizes results from simulation studies.

### Folders:

- **data/**: Contains locally stored data for analysis. Note that this directory is excluded by `.gitignore` to ensure data privacy, so you will need to populate accordingly.
- **figures/**: Stores output visualizations, e.g. scree plots and Sankey diagrams.
- **src/**: Contains the core code for the BIPmixed method and supporting functions.
- **data_analysis_results/**: Directory for output files related to data analysis results.
- **simulation_study_results/**: Directory for output files from simulation studies.

## Workflow

1. **Exploratory Data Analysis (EDA)**:
   - Scripts and documents with the prefix `00_` outline the initial exploration of data.

2. **Data Loading and Processing**:
   - Scripts with the prefix `01_` focus on loading and processing the individual views of data:
     - ELA metrics.
     - Structural MRI (cortical thickness and surface area).
     - Functional MRI (network correlations).
     - Covariates and outcomes.

3. **Sample Definition**:
   - `02_combine_data_define_sample.R` merges processed views into a unified dataset, ensuring consistent samples across views.

4. **Hyperparameter Selection**:
   - `03_generate_scree_plots.R` calculates eigenvalues of concatenated views and determines the number of latent factors using a scree plot.

5. **BIPmixed Analysis**:
   - `data_analysis.R` applies the BIPmixed framework to the processed data, performing feature selection and outcome modeling.
   - Results are summarized in `data_analysis_summary.R`.

6. **Simulation Studies**:
   - `simulation_study.R` simulates multiview data under varying random effect scenarios to evaluate BIPmixed and alternative methods.
   - Results are summarized in `simulation_study_summary.R`.

## Key Features

- **Exploratory Data Analysis**:
  - Initial data exploration guides decisions about feature inclusion, exclusion, and preprocessing.

- **BIPmixed Framework**:
  - Integrates multiview data with hierarchical structures.
  - Simultaneously performs feature selection and outcome prediction.

- **Simulation Studies**:
  - Compare BIPmixed against baseline BIP, PCA2Step, and Cooperative Learning methods across different hierarchical scenarios discussed in the manuscript.

## Outputs

- **Figures**:
  - Where scree plots for latent factor selection and Sankey plots showing feature-to-component mapping are written.

- **Summary Tables**:
  - Report feature selection performance (e.g., FPR, FNR, AUC) and prediction accuracy (e.g., MSE), and are written to data_analysis_results or simulation_study_results.

## Instructions for Reproduction

To reproduce the analyses or simulations, ensure the following:

1. **Data**:
   - Place required data files in the `data/` directory. 

2. **Execution**:
   - Run scripts sequentially according to the workflow described above, installing the required R packages as you go.

## Contact

For questions or collaborations, please contact Aidan Neher at [neher015@umn.edu](mailto:neher015@umn.edu).

---

The authors gratefully acknowledge the research support provided by the NIH T32 Interdisciplinary Biostatistics Training in Genetics and Genomics program (T32 GM132063, 2020–2025). Thierry Chekouo was supported by a National Institutes of Health (NIH)  grant: 1R35GM150537-01. Thierry Chekouo also thanks Medtronic Inc. for their support in the form of a faculty fellowship. Sandra Safo was supported by NIH NIGMS grant award number 1R35GM142695.

