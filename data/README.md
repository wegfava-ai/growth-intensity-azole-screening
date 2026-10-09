# Input data

The original LATAsp databases are not included in this repository.

## File expected by the analysis

Create:

```text
data/analysis_input.csv
```

The file should contain one row per de-identified isolate record in the project-level analysis dataset.

To reproduce the manuscript flow counts, retain records even when complete screening or complete MIC results are unavailable. The score-based analyses automatically restrict to isolates with complete three-azole agar screening, complete three-azole MIC data, and interpretable WT/NWT classification.

## Required columns

| Column | Required | Description |
|---|---:|---|
| `isolate_id` | Yes | De-identified unique isolate identifier. Do not use participant names, record numbers, or other direct identifiers. |
| `cluster_id` | Yes | De-identified identifier for the citizen/source cluster used for cluster bootstrap. Repeated isolates from the same cluster must have the same value. |
| `country` | Yes for geographic summaries | Country of collection. Manuscript countries: Argentina, Brazil, Costa Rica, Guatemala, Mexico, Uruguay. |
| `collection_year` | No | Year of collection, if available. |
| `score_itc` | When ITC screening is available | Itraconazole agar growth score: 0, 1, 2, or 3. |
| `score_vrc` | When VRC screening is available | Voriconazole agar growth score: 0, 1, 2, or 3. |
| `score_pos` | When POS screening is available | Posaconazole agar growth score: 0, 1, 2, or 3. |
| `mic_itc` | When ITC MIC is available | Itraconazole MIC in mg/L. May be numeric (`1`, `2`) or censored (`>8`, `<=0.03`). |
| `mic_vrc` | When VRC MIC is available | Voriconazole MIC in mg/L. May be numeric or censored. |
| `mic_pos` | When POS MIC is available | Posaconazole MIC in mg/L. May be numeric or censored. |

Blank cells are treated as missing.

## Agar growth-score definitions

The analysis assumes the published ordinal scale used in the manuscript:

- `0`: no visible growth;
- `1`: weak/minimal growth at the inoculation site;
- `2`: growth less than the drug-free growth control;
- `3`: uninhibited growth comparable with the control.

The original operational screening rule is any visible growth (`score >= 1`).

## Reference classification implemented in the script

### Itraconazole
CLSI M57S ECV = 1 mg/L:

- WT: MIC <= 1 mg/L
- NWT: MIC > 1 mg/L

### Voriconazole
CLSI M57S ECV = 1 mg/L:

- WT: MIC <= 1 mg/L
- NWT: MIC > 1 mg/L

### Posaconazole
Complementary EUCAST v11.0 interpretation applied to CLSI-obtained MICs:

- WT/S-equivalent: MIC <= 0.125 mg/L
- NWT/R-equivalent: MIC > 0.25 mg/L
- MIC = 0.25 mg/L:
  - ITC WT -> POS WT
  - ITC NWT -> POS NWT

The labels WT and NWT are used analytically for consistency across figures and tables; the POS categorization is not a CLSI ECV-based WT/NWT designation.

## Template

`input_template.csv` contains **synthetic example records only**. It illustrates formatting and can be used to test the pipeline. It must not be used to reproduce or infer manuscript results.

`input_dictionary.csv` provides a machine-readable variable dictionary.
