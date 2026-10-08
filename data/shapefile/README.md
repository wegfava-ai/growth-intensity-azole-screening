# Shapefile required for the supplementary map

The original shapefile used during manuscript development is not distributed in this repository.

To reproduce the supplementary map, supply a country-level Latin America shapefile for which you have permission to use and redistribute locally.

Place all components in this folder, for example:

```text
data/shapefile/LatinAmerica.shp
data/shapefile/LatinAmerica.dbf
data/shapefile/LatinAmerica.shx
data/shapefile/LatinAmerica.prj
```

Optional files such as `.cpg` may also be present.

By default, `scripts/03_supplementary_map.R` expects the country-name attribute to be:

```text
CNTRY_NAME
```

If your shapefile uses another field such as `NAME`, `ADMIN`, or `country`, edit:

```r
COUNTRY_FIELD <- "CNTRY_NAME"
```

near the top of the map script.

The map is generated from aggregated country-level WT/NWT counts produced by `scripts/01_final_analysis.R`. No individual coordinates are required or expected.
