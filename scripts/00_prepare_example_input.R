# Optional helper for testing the pipeline with synthetic example data.
# This DOES NOT reproduce manuscript results.
file.copy(
  from = file.path("data", "input_template.csv"),
  to = file.path("data", "analysis_input.csv"),
  overwrite = TRUE
)
message("Synthetic example copied to data/analysis_input.csv")
