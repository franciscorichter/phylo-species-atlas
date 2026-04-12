# Figure 1: Phylogenetic coverage — 3-layer bars with smooth frontier curve
# Black (described) | Purple (undated tree) | Green (dated tree)
# A smooth curve connects the right edges of all bars, marking the
# boundary between described and undescribed species.
# "Not yet represented" groups are excluded from the figure.

library(ggplot2)
library(scales)
library(ggimage)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
script_dir <- normalizePath(dirname(script_path), mustWork = FALSE)
outdir <- normalizePath(file.path(script_dir, "..", "..", "paper", "figures"),
                        mustWork = FALSE)
data_dir <- normalizePath(file.path(script_dir, "..", "..", "data",
                          "processed_manuscript_inputs"), mustWork = FALSE)

# ── Data ──
d <- read.csv(file.path(data_dir, "fig1_coverage_bars.csv"),
              stringsAsFactors = FALSE)

# Remove "Not yet represented" groups
d <- d[d$category != "Not yet represented", ]

# Percentages relative to estimated_total (= 100% = full bar)
d$desc_pct     <- d$described / d$estimated_total * 100
d$dated_pct    <- d$in_dated_tree / d$estimated_total * 100
d$undated_pct  <- d$in_undated_tree / d$estimated_total * 100
d$has_tree     <- (d$in_dated_tree + d$in_undated_tree) > 0
d$is_dated     <- d$in_dated_tree > 0
d$tree_tips    <- d$in_dated_tree + d$in_undated_tree

# Coverage = tips / described
d$coverage <- ifelse(d$has_tree & d$described > 0,
                     pmin(d$tree_tips / d$described * 100, 100), NA)

# ── Sort ──
cat_order <- c("Vertebrates", "Plants", "Arthropods",
               "Other animals", "Microbes & protists")
d$category <- factor(d$category, levels = cat_order)
# Sort: most black on top = lowest desc_pct first, then lowest tree coverage
d <- d[order(d$category, d$desc_pct, d$dated_pct + d$undated_pct), ]

# ── Assign y positions ──
gap <- 1.8
y_pos <- numeric(nrow(d))
current_y <- 1
prev_cat <- as.character(d$category[1])
for (i in seq_len(nrow(d))) {
  cat_i <- as.character(d$category[i])
  if (cat_i != prev_cat) {
    current_y <- current_y + gap
    prev_cat <- cat_i
  }
  y_pos[i] <- current_y
  current_y <- current_y + 1
}
d$y <- y_pos

# ── Category labels ──
cat_labels <- data.frame(
  category = character(), y_mid = numeric(), y_sep = numeric(),
  stringsAsFactors = FALSE)
for (cc in cat_order) {
  idx <- which(d$category == cc)
  if (length(idx) > 0) {
    cat_labels <- rbind(cat_labels, data.frame(
      category = cc,
      y_mid = mean(d$y[idx]),
      y_sep = max(d$y[idx]) + gap / 2 + 0.3,
      stringsAsFactors = FALSE))
  }
}

# ── Right-side labels ──
fmt_num <- function(x, w = 9) {
  formatC(x, format = "d", big.mark = ",", width = w)
}
d$tips_label <- ifelse(d$has_tree, fmt_num(d$tree_tips), "     --")
d$desc_label <- fmt_num(d$described)
d$cov_label <- ifelse(
  !is.na(d$coverage),
  formatC(sprintf("%.1f%%", round(d$coverage, 1)), width = 7), "")
d$cov_col <- ifelse(is.na(d$coverage), "grey60",
             ifelse(d$coverage > 80, "#1a9850",
             ifelse(d$coverage > 20, "#e6960c", "#d73027")))
d$est_label <- fmt_num(d$estimated_total)

# ── Icon paths ──
icon_dir <- normalizePath(file.path(outdir, "icons"), mustWork = FALSE)
d$icon <- file.path(icon_dir, paste0(
  gsub("[. ]+", "_", tolower(d$group)), ".png"))
d$icon_exists <- file.exists(d$icon)

# ── Bar parameters ──
bar_h <- 0.6
col_desc     <- "#1a1a1a"    # black — described species (the bar itself)
col_undated  <- "#5e3c99"    # dark purple — in undated tree
col_dated    <- "#1a9850"    # dark green — in dated tree

# ── Plot ──
p <- ggplot(d) +
  # Category separators
  geom_hline(yintercept = cat_labels$y_sep[-nrow(cat_labels)],
             color = "grey80", linewidth = 0.3) +

  # Category labels
  annotate("text", x = -30, y = cat_labels$y_mid,
           label = cat_labels$category, fontface = "bold.italic",
           size = 3.2, color = "grey35", hjust = 0.5) +

  # Layer 1: Black bar — described species (flat right edge)
  geom_rect(aes(
    xmin = 0, xmax = desc_pct,
    ymin = y - bar_h / 2, ymax = y + bar_h / 2),
    fill = col_desc, color = NA) +

  # Layer 2: Purple overlay — undated tree tips
  geom_rect(
    data = d[d$in_undated_tree > 0, ],
    aes(xmin = 0, xmax = undated_pct,
        ymin = y - bar_h / 2, ymax = y + bar_h / 2),
    fill = col_undated, color = NA) +

  # Layer 3: Green overlay — dated tree tips
  geom_rect(
    data = d[d$in_dated_tree > 0, ],
    aes(xmin = 0, xmax = dated_pct,
        ymin = y - bar_h / 2, ymax = y + bar_h / 2),
    fill = col_dated, color = NA) +

  # Group icons
  geom_image(data = d[d$icon_exists, ],
             aes(x = -3, y = y, image = icon),
             size = 0.014, asp = 1.3) +

  # Group name labels
  geom_text(aes(x = -6, y = y, label = group),
            hjust = 1, size = 2.7, color = "grey20") +

  # Right side labels
  geom_text(aes(x = 108, y = y, label = tips_label),
            hjust = 1, size = 2.1, color = "grey25", family = "mono") +
  geom_text(aes(x = 109.5, y = y, label = "/"),
            hjust = 0.5, size = 2.1, color = "grey40", family = "mono") +
  geom_text(aes(x = 119, y = y, label = desc_label),
            hjust = 1, size = 2.1, color = "grey25", family = "mono") +
  geom_text(aes(x = 130, y = y, label = cov_label),
            hjust = 1, size = 2.1, fontface = "bold",
            color = d$cov_col, family = "mono") +
  geom_text(aes(x = 145, y = y, label = est_label),
            hjust = 1, size = 2.1, color = "grey25", family = "mono") +

  # X axis
  scale_x_continuous(
    breaks = seq(0, 100, 25),
    labels = paste0(seq(0, 100, 25), "%"),
    limits = c(-42, 150),
    name = "Proportion of estimated total species"
  ) +
  scale_y_continuous(breaks = NULL, expand = expansion(add = c(1, 2))) +

  # Column headers
  annotate("text", x = 108, y = max(y_pos) + 2,
           label = "in tree", size = 2.1,
           color = "grey35", fontface = "italic", hjust = 1) +
  annotate("text", x = 119, y = max(y_pos) + 2,
           label = "described", size = 2.1,
           color = "grey35", fontface = "italic", hjust = 1) +
  annotate("text", x = 130, y = max(y_pos) + 2,
           label = "coverage", size = 2.1,
           color = "grey35", fontface = "italic", hjust = 1) +
  annotate("text", x = 145, y = max(y_pos) + 2,
           label = "estimated", size = 2.1,
           color = "grey35", fontface = "italic", hjust = 1) +

  # Reference lines
  geom_vline(xintercept = c(25, 50, 75), color = "white",
             linewidth = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 100, color = "grey50",
             linewidth = 0.3, linetype = "dashed") +

  # Legend (manual)
  annotate("rect", xmin = 0, xmax = 3, ymin = -2.5, ymax = -1.7,
           fill = col_desc) +
  annotate("text", x = 4, y = -2.1, label = "Described, not in tree",
           hjust = 0, size = 2.5, color = "grey30") +
  annotate("rect", xmin = 35, xmax = 38, ymin = -2.5, ymax = -1.7,
           fill = col_undated) +
  annotate("text", x = 39, y = -2.1, label = "In undated tree",
           hjust = 0, size = 2.5, color = "grey30") +
  annotate("rect", xmin = 62, xmax = 65, ymin = -2.5, ymax = -1.7,
           fill = col_dated) +
  annotate("text", x = 66, y = -2.1, label = "In dated tree",
           hjust = 0, size = 2.5, color = "grey30") +
  labs(
    title = "Phylogenetic coverage across the tree of life",
    y = NULL
  ) +

  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 9),
    legend.position = "none",
    plot.margin = margin(10, 20, 20, 10)
  )

ggsave(file.path(outdir, "figure_nested_circles_v2.png"), p,
       width = 12, height = 16, dpi = 300)
cat("Figure 1 (3-layer bars with frontier curve) saved.\n")
