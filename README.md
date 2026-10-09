# Growth intensity in azole screening

Reproducible R code supporting the manuscript:

**Impact of growth intensity on the performance of azole agar screening in environmental *Aspergillus fumigatus* isolates**

This repository contains the analysis workflow and figure-generation code used for the manuscript. The original project database, raw laboratory files, participant-linked metadata, and geographic shapefiles are **not included** because they may contain sensitive, restricted, or third-party data.

## Repository structure

```text
growth-intensity-azole-screening/
├── README.md
├── .gitignore
├── run_all.R
├── scripts/
│   ├── 00_prepare_example_input.R
│   ├── 01_final_analysis.R
│   ├── 02_publication_figures.R
│   └── 03_supplementary_map.R
├── data/
│   ├── README.md
│   ├── input_template.csv
│   ├── input_dictionary.csv
│   └── shapefile/
│       └── README.md
└── output/
    ├── tables/
    └── figures/
```

## Analyses reproduced by this repository

The public workflow reproduces the analyses reported in the manuscript:

- derivation of the paired analysis set;
- country-level WT/NWT summary;
- exact agar growth score (0–3) by WT/NWT classification;
- NWT proportion at each exact score with exact binomial 95% confidence intervals;
- score contrasts: score 1 vs score 0 and score >=2 vs score 0;
- contribution of score 1 to screening-positive/WT discordance;
- secondary screening-performance estimates at thresholds >=1, >=2 and 3;
- conventional 2 x 2 cross-classification at the original any-growth threshold;
- cluster-bootstrap AUC for the ordinal score;
- manuscript Figure 2, Figure 3, Supplementary Figure S2 (optional map), and Supplementary Figure S3.

Analyses that are not reported in the final manuscript are intentionally not included in the public workflow.

## Required R version and packages

The scripts were designed for R 4.x and use:

- dplyr
- tidyr
- stringr
- readr
- purrr
- tibble
- ggplot2
- scales
- pROC
- writexl
- viridis
- ggalluvial
- sf
- scatterpie
- ggrepel

Missing packages are installed automatically by the scripts.

## Input data

The public repository does not contain the study dataset.

To reproduce the analyses with an authorized dataset, create:

```text
data/analysis_input.csv
```

using exactly the column structure documented in:

- `data/input_template.csv`
- `data/input_dictionary.csv`
- `data/README.md`

To reproduce the manuscript flow counts, the input should represent the full project-level analysis dataset, including records without complete paired screening and MIC results. The score-based analyses themselves use only isolates meeting the paired-analysis criteria.

The required variables are:

- de-identified isolate ID;
- de-identified citizen/source cluster ID;
- country;
- optional collection year;
- agar growth scores for itraconazole (ITC), voriconazole (VRC), and posaconazole (POS);
- MIC values for ITC, VRC, and POS.

Do not commit the real `analysis_input.csv` unless you have explicit permission to make those data public.

## Shapefile

The supplementary geographic map is optional and requires a country-level Latin America shapefile. The shapefile is not distributed in this repository.

Place all shapefile components in:

```text
data/shapefile/
```

with the main file named:

```text
LatinAmerica.shp
```

The attribute table must contain a country-name field called `CNTRY_NAME`, or `COUNTRY_FIELD` must be edited in `scripts/03_supplementary_map.R`.

See `data/shapefile/README.md` for details.

## Running the workflow

### Option 1 — run everything

From the repository root:

```r
source("run_all.R")
```

`run_all.R` performs the final analysis, generates the data-derived publication figures, and generates the supplementary map only if the required shapefile is available.

### Option 2 — run scripts separately

```r
source("scripts/01_final_analysis.R")
source("scripts/02_publication_figures.R")
source("scripts/03_supplementary_map.R") # optional
```

To test the workflow with the synthetic template:

```r
source("scripts/00_prepare_example_input.R")
source("run_all.R")
```

The synthetic example is provided only to demonstrate the expected file structure and code execution. It does **not** reproduce or approximate the manuscript results.

## Statistical definitions used in the manuscript

The final analysis uses:

- exact growth score 0–3 as the primary exposure;
- WT/NWT classification based on broth microdilution MIC results;
- ITC and VRC ECV = 1 mg/L;
- complementary EUCAST v11.0 interpretation for POS;
- exact binomial 95% confidence intervals for score-specific NWT proportions;
- citizen/source-level cluster bootstrap;
- 5,000 bootstrap replicates;
- random seed `20261005`;
- secondary screening-performance estimates at score thresholds >=1, >=2 and 3.

The secondary performance estimates describe the paired analysis set and should not be interpreted as population-level diagnostic-accuracy estimates.

## Figure numbering

The analysis-generated figures correspond to the final manuscript as follows:

- `Figure_2_NWT_proportion_by_score.*` → **Figure 2**
- `Figure_3_heatmap_score_azole_NWT.*` → **Figure 3**
- `Supplementary_Figure_S2_map_WT_NWT.*` → **Supplementary Figure S2**
- `Supplementary_Figure_S3_alluvial_score_to_WT_NWT.*` → **Supplementary Figure S3**

**Figure 1** is the study-flow diagram and **Supplementary Figure S1** contains representative laboratory photographs of the growth scores; these are not generated from the analysis dataset and therefore are not recreated by `02_publication_figures.R`.

## Privacy and data protection

The `.gitignore` file is configured to ignore:

- `data/analysis_input.csv`;
- spreadsheet/database exports;
- REDCap exports;
- shapefiles and their sidecar files;
- generated output files.

This reduces the risk of accidentally committing restricted study data.

## Data availability

The R scripts used for data processing, statistical analyses, and figure generation are publicly available in this repository. The original study data may be shared upon reasonable request to the corresponding author.

## AI-assisted review
Artificial intelligence tools, including OpenAI ChatGPT (GPT-5.6 Sol), were used to assist with review of the R scripts, including checks for code organization, consistency, documentation, and reproducibility. All analytical decisions, statistical methods, interpretation of results, and final verification of the code and outputs were performed by the authors.

## Citation

Update this section with the final article citation and DOI after publication.
