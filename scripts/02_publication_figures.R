# ============================================================
# Publication figures
# ============================================================
#
# Requires outputs from:
#   scripts/01_final_analysis.R
#
# Produces manuscript figures generated from analysis outputs:
#   Figure_2_NWT_proportion_by_score
#   Figure_3_heatmap_score_azole_NWT
#   Supplementary_Figure_S3_alluvial_score_to_WT_NWT
#
# Figure 1 is the study-flow diagram and Supplementary Figure S1
# contains representative laboratory photographs; neither is generated
# by this script.
# ============================================================

pkgs <- c("ggplot2", "dplyr", "tidyr", "scales", "viridis", "ggalluvial", "grid")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

tables_dir <- file.path("output", "tables")
figures_dir <- file.path("output", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

nwt_file <- file.path(tables_dir, "Primary_NWT_proportion_by_score_exact95CI.csv")
count_file <- file.path(tables_dir, "Supplement_score_x_WT_NWT_counts.csv")

if (!file.exists(nwt_file) || !file.exists(count_file)) {
  stop("Required analysis tables were not found. Run scripts/01_final_analysis.R first.")
}

nwt_score <- read.csv(nwt_file, check.names = FALSE)
score_ref <- read.csv(count_file, check.names = FALSE)

azole_levels <- c("Itraconazole", "Voriconazole", "Posaconazole")

nwt_score <- nwt_score |>
  dplyr::mutate(
    azole = factor(azole, levels = azole_levels),
    score = factor(score, levels = 0:3)
  )

score_ref <- score_ref |>
  dplyr::mutate(
    azole = factor(azole, levels = azole_levels),
    score = factor(score, levels = 0:3),
    reference = factor(reference, levels = c("WT", "NWT"))
  )

theme_latasp <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 0.5),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 1.3),
      plot.subtitle = element_text(size = base_size - 0.2),
      plot.caption = element_text(
        size = base_size - 2.2,
        colour = "grey35",
        hjust = 0
      ),
      plot.margin = margin(8, 12, 8, 10)
    )
}

save_publication <- function(plot, filename, width, height) {
  ggsave(
    file.path(figures_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    file.path(figures_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
}

# Main Figure: exact NWT proportion by score
fig1 <- ggplot(
  nwt_score,
  aes(x = score, y = proportion_nwt, colour = azole, group = azole)
) +
  geom_line(linewidth = 0.45, alpha = 0.45) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.08,
    linewidth = 0.60
  ) +
  geom_point(size = 3.0) +
  geom_text(
    aes(y = pmin(ci_upper + 0.060, 1.055), label = paste0("n=", n)),
    size = 2.9,
    colour = "grey25",
    show.legend = FALSE
  ) +
  facet_wrap(~ azole, nrow = 1) +
  scale_colour_viridis_d(option = "D", begin = 0.10, end = 0.85, guide = "none") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.08),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Agar growth score",
    y = "NWT isolates (%)",
    title = "Non-wild-type proportion across agar growth scores",
    subtitle = "Observed proportions with exact binomial 95% confidence intervals",
    caption = paste0(
      "Growth score: 0 = no visible growth; 1 = weak growth; ",
      "2 = reduced growth vs control; 3 = uninhibited growth."
    )
  ) +
  theme_latasp(base_size = 10.8)

save_publication(
  fig1,
  "Figure_2_NWT_proportion_by_score",
  width = 9.0,
  height = 4.65
)

# Main Figure: heatmap
heatmap_data <- nwt_score |>
  dplyr::mutate(
    pct_nwt = 100 * proportion_nwt,
    label = paste0(sprintf("%.1f", pct_nwt), "%\n", n_nwt, "/", n),
    azole_heat = factor(as.character(azole), levels = rev(azole_levels)),
    text_colour = ifelse(pct_nwt >= 45, "white", "black")
  )

fig2 <- ggplot(
  heatmap_data,
  aes(x = score, y = azole_heat, fill = pct_nwt)
) +
  geom_tile(colour = "white", linewidth = 1.15, width = 0.98, height = 0.92) +
  geom_text(
    aes(label = label, colour = text_colour),
    size = 3.65,
    lineheight = 0.95,
    fontface = "bold"
  ) +
  scale_colour_identity() +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%")
  ) +
  coord_fixed(ratio = 0.78) +
  labs(
    x = "Agar growth score",
    y = NULL,
    fill = "NWT",
    title = "Non-wild-type proportion by azole and growth score",
    subtitle = "Cells show NWT percentage and NWT/total isolates"
  ) +
  theme_latasp(base_size = 10.8) +
  theme(panel.border = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())

save_publication(
  fig2,
  "Figure_3_heatmap_score_azole_NWT",
  width = 7.6,
  height = 4.75
)

# Supplementary alluvial figure
alluvial_data <- score_ref |>
  dplyr::mutate(
    score_label = factor(paste0("Score ", score), levels = paste0("Score ", 0:3))
  )

score_cols <- viridis::viridis(
  n = 4, option = "D", begin = 0.08, end = 0.88
)
names(score_cols) <- paste0("Score ", 0:3)

fig3 <- ggplot(
  alluvial_data,
  aes(axis1 = score_label, axis2 = reference, y = n)
) +
  ggalluvial::geom_alluvium(
    aes(fill = score_label),
    width = 0.18,
    alpha = 0.76,
    knot.pos = 0.42,
    colour = NA
  ) +
  ggalluvial::geom_stratum(
    width = 0.18,
    fill = "grey97",
    colour = "grey45",
    linewidth = 0.40
  ) +
  ggplot2::geom_text(
    stat = "stratum",
    aes(
      label = ifelse(
        after_stat(x) == 2,
        as.character(after_stat(stratum)),
        ""
      )
    ),
    size = 3.05,
    fontface = "bold"
  ) +
  facet_wrap(~ azole, nrow = 1, scales = "free_y") +
  scale_x_discrete(
    limits = c("Growth score", "MIC classification"),
    expand = c(0.13, 0.08)
  ) +
  scale_fill_manual(values = score_cols, name = "Growth score") +
  labs(
    x = NULL,
    y = "Number of isolates",
    title = "Flow from agar growth score to WT/NWT reference classification",
    subtitle = "Ribbon width is proportional to the number of isolates",
    caption = "WT/NWT classification follows the broth microdilution MIC criteria used in the study."
  ) +
  theme_latasp(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )

save_publication(
  fig3,
  "Supplementary_Figure_S3_alluvial_score_to_WT_NWT",
  width = 10.0,
  height = 5.45
)

message("Publication figures completed.")
