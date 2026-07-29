# The genetic architecture of and evolutionary constraints on the human pelvic form

Code accompanying:

> Xu, L., Kun, E., Pandey, D., Wang, J. Y., Brasil, M. F., Singh, T., &
> Narasimhan, V. M. (2025). **The genetic architecture of and evolutionary
> constraints on the human pelvic form.** *Science*, 388(6743), eadq1521.
> https://doi.org/10.1126/science.adq1521

This repository contains the full analysis pipeline used in the study, from
automatically landmarking full-body DXA images through phenotype extraction,
genome-wide association studies (GWAS), and the downstream genetic and
evolutionary analyses that produced the paper's figures.

---

## Overview

The project answers a simple question with a large pipeline: **which genetic
variants shape the human pelvis, and what evolutionary forces have acted on
it?** The workflow has three stages:

1. **Landmark prediction** — a deep-learning (HRNet) model places 23 anatomical
   landmarks on each imaging-cohort DXA image.
2. **Phenotyping** — pelvic and skeletal measurements (lengths, widths, angles,
   areas, and ratios) are computed from those landmarks and converted from
   pixels to centimeters.
3. **Genetics & evolution** — the phenotypes are run through GWAS, heritability
   and genetic-correlation estimation, Mendelian randomization, gene/tissue
   enrichment, and tests for selection, then visualized.

---

## Repository structure

```
Human-pelvic-form/
├── README.md
└── scripts/
    ├── model_for_landmarks_prediction/   # Stage 1: HRNet landmark model (steps 0–1)
    │   ├── 0_train_data_prepare.ipynb
    │   ├── 1_model_landmarks_pred.ipynb
    │   └── hrnet.py
    ├── data_analysis/                    # Stages 2–3: phenotyping + GWAS prep (steps 2–6)
    │   ├── 2_combine_all_kps_info.ipynb
    │   ├── 3_train_val_data_analysis.ipynb
    │   ├── 4_calculate_phenotype.ipynb
    │   ├── 5_gwas_data_prep.ipynb
    │   ├── 6_data_analysis.ipynb
    │   ├── get_phenotypes_residual.R
    │   ├── disease_regression_pheno.R
    │   ├── disease_regression_prs.R
    │   └── convert_geneid.R
    └── plot/                             # Figure-generation R scripts
        └── ... (see scripts/plot/README.md)
```

Each subdirectory has its own `README.md` describing every file.

---

## Analysis pipeline

The notebooks are numbered to reflect execution order, flowing from the model
directory into the analysis directory:

| Step | File | Purpose |
| --- | --- | --- |
| 0 | `model_for_landmarks_prediction/0_train_data_prepare.ipynb` | Build the training set; convert COCO annotations into 23 keypoints. |
| 1 | `model_for_landmarks_prediction/1_model_landmarks_pred.ipynb` | Predict landmarks on all images with the trained HRNet. |
| 2 | `data_analysis/2_combine_all_kps_info.ipynb` | Merge image metadata with model predictions into per-individual keypoints. |
| 3 | `data_analysis/3_train_val_data_analysis.ipynb` | Validate predictions, extract pixel phenotypes, convert to cm. |
| 4 | `data_analysis/4_calculate_phenotype.ipynb` | Compute all pelvic/skeletal phenotypes and calibrate pixel→cm. |
| 5 | `data_analysis/5_gwas_data_prep.ipynb` | Build covariate files, clean data, and prepare GWAS/heritability inputs. |
| 6 | `data_analysis/6_data_analysis.ipynb` | Phenotype correlations, sex comparisons, and clustering. |

Supporting R scripts in `data_analysis/` (residualization, disease/PRS
regression, gene-ID conversion) and all figure scripts in `plot/` run on the
outputs of these steps. See the per-directory READMEs for details.

---

## Requirements

The pipeline mixes Python (imaging, phenotyping) and R (statistical genetics,
figures), and relies on several external command-line genetics tools.

**Python** (notebooks, HRNet): `torch` (PyTorch), `fastai`, `numpy`, `pandas`,
`scipy`, `scikit-learn`, `matplotlib`, `seaborn`, `Pillow`, `pydicom`,
`beautifulsoup4`, `requests`.

**R** (analysis & plotting): `tidyverse` (incl. `dplyr`, `ggplot2`),
`data.table`, `biomaRt`, `pheatmap`, `ComplexHeatmap`, `circlize`, `reshape2`,
`Hmisc`, `patchwork`, `ggpubr`, `ggsci`, `ggthemes`, `RColorBrewer`, `viridis`,
`qqman`, `TwoSampleMR`, `MRInstruments`.

**External genetics tools** used in the study (run outside these scripts, with
their outputs consumed here): [BOLT-LMM](https://alkesgroup.broadinstitute.org/BOLT-LMM/),
[PLINK](https://www.cog-genomics.org/plink/),
[LDSC / S-LDSC](https://github.com/bulik/ldsc),
[GCTA](https://yanglab.westlake.edu.cn/software/gcta/),
[MAGMA](https://cncr.nl/research/magma/), and
[LCV](https://github.com/lukejoconnor/LCV).

---

## Data availability

The individual-level imaging and genetic data used here are from the
[UK Biobank](https://www.ukbiobank.ac.uk/) and are **access-controlled** — they
are not, and cannot be, distributed in this repository. Access must be obtained
through the UK Biobank application process. GWAS summary statistics and other
derived results are available as described in the paper.

Because the scripts were written for a specific compute environment, they
contain **absolute local file paths** (e.g. `/Users/.../Narasimhan_lab/...`).
These are provided for transparency and reproducibility of the analysis logic;
you will need to adapt the paths to run them on your own data.

---

## Citation

If you use this code, please cite:

```bibtex
@article{xu2025pelvic,
  title   = {The genetic architecture of and evolutionary constraints on the human pelvic form},
  author  = {Xu, Liaoyi and Kun, Eucharist and Pandey, Devansh and Wang, Jiaxue Y. and Brasil, Marianne F. and Singh, Tarjinder and Narasimhan, Vagheesh M.},
  journal = {Science},
  volume  = {388},
  number  = {6743},
  pages   = {eadq1521},
  year    = {2025},
  doi     = {10.1126/science.adq1521}
}
```

## Contact

Questions are welcome — please reach out to Liaoyi Xu by
[email](mailto:xliaoyi@utexas.edu).
