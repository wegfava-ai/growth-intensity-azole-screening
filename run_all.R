# ============================================================
# Run complete public reproducibility workflow
# ============================================================

if (!file.exists("scripts/01_final_analysis.R")) {
  stop("Run this script from the repository root.")
}

source("scripts/01_final_analysis.R")
source("scripts/02_publication_figures.R")

shp <- file.path("data", "shapefile", "LatinAmerica.shp")

if (file.exists(shp)) {
  source("scripts/03_supplementary_map.R")
} else {
  message(
    "Supplementary map skipped: ",
    shp,
    " was not found. See data/shapefile/README.md."
  )
}
