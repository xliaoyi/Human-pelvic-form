# Data analysis

Takes the predicted landmarks from `../model_for_landmarks_prediction`, turns
them into quantitative pelvic and skeletal phenotypes, and prepares those
phenotypes for genome-wide association studies (GWAS) and downstream genetic
analyses.

## Notebooks (run in numbered order)

| File | Description |
| --- | --- |
| `2_combine_all_kps_info.ipynb` | Collect all DICOM/image metadata and load the HRNet (hrnet32) prediction results, converting them into a unified JSON of keypoints per individual. |
| `3_train_val_data_analysis.ipynb` | Analyze validation predictions, inspect outliers, extract phenotypes in pixels, gather image-size metadata, and convert pixels to centimeters. Compares manual annotation vs. model prediction (predictions used as training labels to demonstrate model quality). |
| `4_calculate_phenotype.ipynb` | Compute the full set of pelvic and skeletal phenotypes — lengths, angles, areas, and ratios — from the landmarks, and calibrate the pixel-to-centimeter conversion. |
| `5_gwas_data_prep.ipynb` | Assemble GWAS inputs: extract software versions and device serial numbers, build covariate files for the imaging and 400k white-British cohorts, clean data (z-score filtering), add other skeletal traits, and prepare inputs for heritability (h²) estimation. |
| `6_data_analysis.ipynb` | Phenotype-level analyses: left/right averaging, phenotypic correlations, male/female comparisons, and hierarchical clustering / grouping of phenotypes by correlation. |

## R scripts

| File | Description |
| --- | --- |
| `get_phenotypes_residual.R` | Regress covariates (e.g. standing height) out of the phenotypes to produce residualized phenotypes used as GWAS inputs, for the combined, female, and male cohorts. |
| `disease_regression_pheno.R` | Regress disease / outcome phenotypes (e.g. osteoarthritis, pain, reproductive outcomes) on the measured pelvic phenotypes. |
| `disease_regression_prs.R` | Same disease/outcome associations using polygenic risk scores (PRS) instead of measured phenotypes; merges in genetic PCs and outcome data. |
| `convert_geneid.R` | Map gene symbols to Ensembl gene IDs (ENSG) via `biomaRt` for RNA-seq gene-expression data. |

> **Note:** All scripts use absolute local paths to UK Biobank–derived data that
> is not included in this repository. Adapt the paths to your environment.
