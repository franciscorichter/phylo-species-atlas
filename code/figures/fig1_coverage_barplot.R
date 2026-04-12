# Figure 1: Phylogenetic dark matter — overlaid bars
# Dark bar = total described, colored overlay = species in phylogenies

library(ggplot2)
library(scales)

script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
script_dir <- normalizePath(dirname(script_path), mustWork=FALSE)
outdir <- normalizePath(file.path(script_dir, "..", "..", "paper", "figures"),
                        mustWork = FALSE)
data_dir <- normalizePath(file.path(script_dir, "..", "..", "data",
                          "processed_manuscript_inputs"), mustWork = FALSE)

d <- read.csv(file.path(data_dir, "fig1_coverage_data.csv"),
              stringsAsFactors = FALSE)
d$tips <- pmin(d$tips, d$described)
d$coverage <- d$tips / d$described * 100

# Sort by described species count (largest at top)
d <- d[order(d$described), ]
d$group <- factor(d$group, levels = d$group)

cat_colors <- c(
  "Vertebrates"         = "#377EB8",
  "Arthropods"          = "#E41A1C",
  "Plants"              = "#4DAF4A",
  "Other animals"       = "#FF7F00",
  "Microbes & protists" = "#984EA3"
)

# Right-side labels
d$label <- sprintf("%s%%  (%s / %s)",
  ifelse(d$coverage >= 1, formatC(round(d$coverage, 1), format = "f",
         digits = 1), "<1"),
  formatC(d$tips, format = "d", big.mark = ","),
  formatC(d$described, format = "d", big.mark = ","))

p <- ggplot(d, aes(y = group)) +
  # Dark bar: total described species
  geom_col(aes(x = described), fill = "#1a1a1a", width = 0.7) +
  # Colored overlay: species in phylogenies
  geom_col(aes(x = tips, fill = category), width = 0.7) +
  # Labels
  geom_text(aes(x = described, label = label),
            hjust = -0.03, size = 2.3, color = "grey30") +
  scale_fill_manual(values = cat_colors, name = NULL) +
  scale_x_log10(
    labels = label_comma(),
    breaks = c(100, 1000, 10000, 100000, 1000000),
    expand = expansion(mult = c(0.01, 0.25))
  ) +
  labs(
    title = "The phylogenetic dark matter",
    subtitle = paste0("Colored bars: species in current phylogenies. ",
                      "Dark bars: total described species."),
    x = "Number of species",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7.5),
    plot.margin = margin(10, 80, 10, 5)
  ) +
  guides(fill = guide_legend(nrow = 1))

ggsave(file.path(outdir, "figure1_barplot_final.png"), p,
       width = 12, height = 11, dpi = 300)
ggsave(file.path(outdir, "figure1_barplot_final.pdf"), p,
       width = 12, height = 11)
cat("Figure 1 saved.\n")
