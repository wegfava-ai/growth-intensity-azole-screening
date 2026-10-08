# Growth intensity in azole screening

Reproducible R code supporting the manuscript:

**Impact of growth intensity on the performance of azole agar screening in environmental *Aspergillus fumigatus* isolates**

This repository contains the analysis workflow and figure-generation code used for the manuscript. The original project database, REDCap export, raw laboratory workbook, and geographic shapefiles are **not included** because they may contain sensitive, restricted, or third-party data.

## Repository structure

```text
growth-intensity-azole-screening/
├── README.md
├── .gitignore
├── run_all.R
├── scripts/
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

The analysis input must represent the **full screening cohort**, not only isolates with complete MIC results. This is necessary to reproduce the verification-rate analysis and STARD flow counts.

The minimum variables are:

- de-identified isolate ID
- de-identified cluster ID
- country
- optional collection year
- agar growth scores for ITC, VRC and POS
- MIC values for ITC, VRC and POS

Do not commit the real `analysis_input.csv` unless you have permission to make the data public.

## Shapefile

The supplementary map is optional and requires a Latin America country shapefile. It is not distributed here.

Place all shapefile components in:

```text
data/shapefile/
```

with the main file named:

```text
LatinAmerica.shp
```

The attribute table must contain a country-name field called `CNTRY_NAME`, or you must edit `COUNTRY_FIELD` in `scripts/03_supplementary_map.R`.

See `data/shapefile/README.md`.

## Running the analysis

### Option 1 — run everything

From the repository root:

```r
source("run_all.R")
```

`run_all.R` performs the final analysis, generates publication figures, and generates the map only if the required shapefile is present.

### Option 2 — run scripts separately

```r
source("scripts/01_final_analysis.R")
source("scripts/02_publication_figures.R")
source("scripts/03_supplementary_map.R") # optional
```

## Reproducing the published analysis

The public template contains synthetic example data and is intended only to demonstrate the required format. It does **not** reproduce the manuscript results.

To reproduce the article results, an authorized copy of the de-identified study analysis dataset must be saved as:

```text
data/analysis_input.csv
```

The final analysis uses:

- exact growth score 0–3 as the primary exposure;
- WT/NWT classification from broth microdilution MICs;
- ITC and VRC ECV = 1 mg/L;
- complementary EUCAST v11.0 interpretation for POS;
- exact binomial 95% confidence intervals;
- citizen-level cluster bootstrap;
- 5,000 bootstrap replicates;
- seed `20261005`;
- diagnostic-performance measures as secondary analyses because of partial verification.

## Important methodological note

During the later operational workflow, MIC verification was largely conditional on positive agar screening. Therefore sensitivity, specificity, predictive values, likelihood ratios, accuracy and AUC derived from the paired subset are secondary descriptive measures and should not be interpreted as unbiased population diagnostic-accuracy estimates.

## Privacy and data protection

The `.gitignore` file is configured to ignore:

- `data/analysis_input.csv`
- spreadsheet/database exports
- REDCap exports
- shapefiles and their sidecar files
- generated output files

This reduces the risk of accidentally committing restricted data.

## License

No software license is included in this repository template. Add a license only after the authors/institution have decided the intended reuse terms.

## Contact / citation

When the manuscript is accepted, update this README with the final citation and DOI.
