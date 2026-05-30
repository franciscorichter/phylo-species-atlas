# Figure S1: Coverage gap bubble chart
# x = described species (log scale), y = coverage %, bubble = species gap

library(ggplot2)
library(ggrepel)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
script_dir <- normalizePath(dirname(script_path), mustWork = FALSE)
outdir <- normalizePath(file.path(script_dir, "..", "..", "figures"),
                        mustWork = FALSE)
data_dir <- normalizePath(file.path(script_dir, "..", "..", "inputs"),
                          mustWork = FALSE)

# ── Data ──
d <- read.csv(file.path(data_dir, "fig1_coverage_bars.csv"),
              stringsAsFactors = FALSE)

# Coverage metrics
d$tree_tips <- d$in_dated_tree + d$in_undated_tree
d$has_tree  <- d$tree_tips > 0
d$coverage  <- ifelse(d$has_tree & d$described > 0,
                      pmin(d$tree_tips / d$described * 100, 100), 0)
d$gap       <- d$described - d$tree_tips
d$gap       <- pmax(d$gap, 0)
d$represented <- d$is_unrepresented != "True"

# Only plot partitions with a broadly representative (partition-canonical) tree;
# the unrepresented partitions (no tree, or only sub-clade trees) sit at coverage = 0
# and are summarised separately in Table S3 + S5.
d <- d[d$represented, ]

# Category as factor
cat_order <- c("Vertebrates", "Plants", "Arthropods",
               "Other animals", "Microbes & protists")
d$category <- factor(d$category, levels = cat_order)

# Colors
cat_cols <- c("Vertebrates" = "#4575b4",
              "Plants" = "#1a9850",
              "Arthropods" = "#d73027",
              "Other animals" = "#f46d43",
              "Microbes & protists" = "#7b3294")

p <- ggplot(d, aes(x = described, y = coverage, size = gap,
                    color = category, label = group)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(size = 2.5, max.overlaps = 30, segment.size = 0.3,
                  segment.color = "grey60", show.legend = FALSE) +
  scale_x_log10(labels = scales::label_comma(),
                breaks = c(100, 1000, 10000, 100000, 1000000)) +
  scale_size_area(max_size = 25,
                  breaks = c(1000, 10000, 100000, 500000),
                  labels = scales::label_comma(),
                  name = "Species gap") +
  scale_color_manual(values = cat_cols, name = "Category") +
  labs(x = "Described species (log scale)",
       y = "Phylogenetic coverage (%)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(outdir, "figure_coverage_gaps.png"), p,
       width = 10, height = 6, dpi = 300)
cat("Figure S1 (coverage gap bubble chart) saved.\n")
