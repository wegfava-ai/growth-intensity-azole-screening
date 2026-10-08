# Input data

The original LATAsp databases are not included in this repository.

## File expected by the analysis

Create:

```text
data/analysis_input.csv
```

The file should contain one row per isolate record in the **full screening dataset**.

Do not restrict the input to the paired MIC-tested isolates. The full screened cohort is required to quantify differential reference-method verification.

## Required columns

| Column | Required | Description |
|---|---:|---|
| `isolate_id` | Yes | De-identified unique isolate identifier. Do not use participant names, record numbers, or other direct identifiers. |
| `cluster_id` | Yes | De-identified identifier for the citizen/source cluster used for cluster bootstrap. Repeated isolates from the same cluster must have the same value. |
| `country` | Yes for geographic summaries | Country of collection. Expected manuscript countries: Argentina, Brazil, Costa Rica, Guatemala, Mexico, Uruguay. |
| `collection_year` | No | Year of collection, if available. |
| `score_itc` | Yes when ITC screening performed | Itraconazole agar growth score: 0, 1, 2, or 3. |
| `score_vrc` | Yes when VRC screening performed | Voriconazole agar growth score: 0, 1, 2, or 3. |
| `score_pos` | Yes when POS screening performed | Posaconazole agar growth score: 0, 1, 2, or 3. |
| `mic_itc` | Required when ITC MIC available | Itraconazole MIC in mg/L. May be numeric (`1`, `2`) or censored (`>8`, `<=0.03`). |
| `mic_vrc` | Required when VRC MIC available | Voriconazole MIC in mg/L. May be numeric or censored. |
| `mic_pos` | Required when POS MIC available | Posaconazole MIC in mg/L. May be numeric or censored. |

Blank cells are treated as missing.

## Agar growth-score definitions

The analysis assumes the published ordinal scale used in the manuscript:

- `0`: no visible growth
- `1`: weak/minimal growth at the inoculation site
- `2`: growth less than the growth-control compartment
- `3`: uninhibited growth comparable with the control

The prespecified operational screening rule is any visible growth (`score >= 1`).

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

- WT/S: MIC <= 0.125 mg/L
- NWT/R: MIC > 0.25 mg/L
- MIC = 0.25 mg/L:
  - ITC WT -> POS WT
  - ITC NWT -> POS NWT

## Template

`input_template.csv` contains **synthetic example records only**. It illustrates formatting and can be used to test the pipeline. It must not be used to reproduce or infer manuscript results.

`input_dictionary.csv` provides a machine-readable variable dictionary.
