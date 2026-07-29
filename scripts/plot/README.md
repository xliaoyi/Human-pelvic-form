# Plotting & figures

R scripts that reproduce the figures in the paper from the analysis outputs
(phenotypes, GWAS summary statistics, heritability/genetic-correlation
estimates, and post-GWAS results). Each script is self-contained and reads its
inputs from result files produced by the `../data_analysis` pipeline and the
external GWAS tools listed in the top-level README.

## Model validation

| Script | Figure content |
| --- | --- |
| `anno_vs_pred.R` | Landmark error: manual annotation vs. model prediction, per hip region. |
| `model_compare_pred_vs_anno.R` | Bar/point comparison of annotation-vs-prediction error across landmarks. |
| `model_compare_with_discrepancy.R` | Cross-visit reproducibility of the model compared against Kun et al., 2023. |
| `left_right_two_visits_correlation.R` | Correlation of left/right phenotypes between the two imaging visits (reproducibility). |

## Phenotype analyses

| Script | Figure content |
| --- | --- |
| `phenotype_compare.R` | Phenotype distributions compared between males and females (boxplots + t-tests). |
| `phenotypic_correlation.R` | Phenotypic correlation heatmap among pelvic phenotypes. |
| `female_phenotype_association.R` | Female-specific associations with hormones and reproductive traits. |
| `left_right_ratio_vs_handedness.R` | Left/right asymmetry ratio vs. handedness. |

## Heritability & genetic correlation

| Script | Figure content |
| --- | --- |
| `genetic_correlation.R` | Genetic correlation heatmap among pelvic phenotypes (both sexes). |
| `genetic_phenotypic_cor.R` | Combined genetic vs. phenotypic correlation (upper/lower triangle of one matrix). |
| `plot_genetic_phenotypic_cor.R` | Genetic vs. phenotypic correlation with Bonferroni significance. |
| `genetic_cor_between_pelvic_IDP_and_outcome_phenotypes.R` | Genetic correlation between pelvic image-derived phenotypes and outcome phenotypes. |
| `male_female_specific_h2g_rg.R` | Sex-specific heritability (h²) and genetic correlation (rg). |
| `male_vs_female_genetic_correlation.R` | Cross-sex (male vs. female) genetic correlation. |

## GWAS & post-GWAS

| Script | Figure content |
| --- | --- |
| `gwas_plot.R` | Manhattan and QQ plots, with genomic inflation factor (lambda). |
| `manhattan_and_lambda.R` | Combined multi-phenotype Manhattan plot and lambda. |
| `magma_heatmap.R` | MAGMA gene-set / tissue enrichment heatmap. |
| `SLDSC_heatmap_atacseq.R` | Stratified LD-score regression (S-LDSC) heritability enrichment for ATAC-seq annotations. |
| `SLDSC_heatmap_evo_anno.R` | S-LDSC heritability enrichment for evolutionary annotations. |
| `LCV_analysis.R` | Latent Causal Variable (LCV) analysis of genetic causality. |
| `mendelian_randomization.R` | Two-sample Mendelian randomization (`TwoSampleMR`). |

## Disease / outcome associations

| Script | Figure content |
| --- | --- |
| `disease_reg_pheno_and_prs_plot.R` | Disease/outcome associations (osteoarthritis, walking, pregnancy) for both measured phenotypes and PRS. |

> **Note:** Scripts use absolute local paths to result files that are not
> included in this repository. Adapt the paths to your environment.
