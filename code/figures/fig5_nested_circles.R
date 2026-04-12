# Figure 5: Nested circles — three levels per group
library(ggplot2)
library(scales)

script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
script_dir <- normalizePath(dirname(script_path), mustWork=FALSE)
outdir <- normalizePath(file.path(script_dir, "..", "..", "paper", "figures"), mustWork=FALSE)
data_dir <- normalizePath(file.path(script_dir, "..", "..", "data", "processed_manuscript_inputs"), mustWork=FALSE)

# ── Data ──
all_g <- read.csv(file.path(data_dir, "fig5_circle_data.csv"), stringsAsFactors=FALSE)

all_g$plot_total <- ifelse(all_g$has_denom, all_g$described, all_g$tips)
all_g$coverage <- ifelse(all_g$has_tree & all_g$has_denom,
                         round(pmin(all_g$tips / all_g$described * 100, 100), 1),
                         NA)

# ── Order by category then described ──
cat_order <- c("Vertebrates","Plants","Arthropods","Other animals","Microbes & protists","No broad tree")
all_g$category <- factor(all_g$category, levels=cat_order)
all_g <- all_g[order(all_g$category, all_g$plot_total), ]

# ── Assign y positions WITH gaps between categories ──
gap <- 2  # empty rows between categories
y_pos <- numeric(nrow(all_g))
current_y <- 1
prev_cat <- as.character(all_g$category[1])
for (i in seq_len(nrow(all_g))) {
  cat_i <- as.character(all_g$category[i])
  if (cat_i != prev_cat) {
    current_y <- current_y + gap
    prev_cat <- cat_i
  }
  y_pos[i] <- current_y
  current_y <- current_y + 1
}
all_g$y <- y_pos

# ── Category label positions (centered in each group's y range) ──
cat_labels <- data.frame(
  category = character(), y_mid = numeric(), y_top = numeric(), stringsAsFactors=FALSE)
for (cc in cat_order) {
  idx <- which(all_g$category == cc)
  if (length(idx) > 0) {
    cat_labels <- rbind(cat_labels, data.frame(
      category = cc,
      y_mid = mean(all_g$y[idx]),
      y_top = max(all_g$y[idx]) + gap/2 + 0.3,
      stringsAsFactors=FALSE))
  }
}

# ── Labels ──
all_g$count_label <- ifelse(
  all_g$has_tree & all_g$has_denom,
  paste0(format(all_g$tips, big.mark=","), " / ", format(all_g$described, big.mark=",")),
  ifelse(all_g$has_tree,
         paste0(format(all_g$tips, big.mark=","), " / --"),
         paste0("-- / ", format(all_g$described, big.mark=","))))

all_g$cov_label <- ifelse(!is.na(all_g$coverage), paste0(all_g$coverage, "%"), "")
all_g$cov_col <- ifelse(is.na(all_g$coverage), "grey60",
                 ifelse(all_g$coverage > 80, "#1a9850",
                 ifelse(all_g$coverage > 20, "#e6960c", "#d73027")))

# ── Colors ──
col_dated   <- "#1a9850"
col_undated <- "#7fcd91"
col_grey    <- "#C0C0C0"
col_black   <- "#1a1a1a"

# ── Plot ──
p <- ggplot(all_g) +
  # Category separator lines
  geom_hline(yintercept = cat_labels$y_top[-nrow(cat_labels)],
             color="grey75", linewidth=0.3, linetype="solid") +

  # Category labels — left side, rotated or bold
  annotate("text", x=log10(55), y=cat_labels$y_mid,
           label=cat_labels$category, fontface="bold.italic", size=3.5,
           color="grey35", hjust=0.5, angle=0) +

  # Outer circles
  geom_point(aes(x=log10(plot_total), y=y, size=plot_total,
                 color=ifelse(has_tree, "a_described", "d_no_tree")), alpha=0.4) +

  # Inner circles
  geom_point(data=all_g[all_g$has_tree & all_g$tips > 0, ],
             aes(x=log10(plot_total), y=y, size=tips,
                 color=ifelse(dated, "b_dated", "c_undated")), alpha=0.9) +

  # Group name labels — offset from category labels
  geom_text(aes(x=log10(150), y=y, label=group),
            hjust=0, size=3, color="grey20") +

  # Species counts
  geom_text(aes(x=log10(2500000), y=y, label=count_label),
            hjust=1, size=2.5, color="grey45", family="mono") +

  # Coverage %
  geom_text(aes(x=log10(4500000), y=y, label=cov_label),
            hjust=1, size=2.5, fontface="bold", color=all_g$cov_col) +

  # Scales
  scale_size_area(max_size=32, labels=comma, name="Species",
                  breaks=c(1000,10000,100000,500000),
                  guide=guide_legend(order=2)) +
  scale_color_manual(
    values=c("a_described"=col_grey, "b_dated"=col_dated,
             "c_undated"=col_undated, "d_no_tree"=col_black),
    labels=c("Described species","In dated tree","In undated tree","No broad tree in atlas"),
    name="",
    guide=guide_legend(order=1, override.aes=list(size=7, alpha=0.8))) +
  scale_x_continuous(
    breaks=log10(c(1000,10000,100000,1000000)),
    labels=c("1,000","10,000","100,000","1,000,000"),
    limits=c(log10(10), log10(5000000)),
    name=NULL) +
  scale_y_continuous(breaks=NULL, expand=expansion(add=c(1, 1.5))) +

  # Labels
  labs(
    title="Phylogenetic coverage across the tree of life",
    subtitle=paste0(
      "Outer circle = described species where available (otherwise tree size)  |  ",
      "Dark green = in dated tree  |  Light green = in undated tree  |  ",
      "Black = no broad tree in atlas"),
    y=NULL) +

  # Column headers
  annotate("text", x=log10(2500000), y=max(y_pos)+2.5,
           label="in tree / described", size=2.5, color="grey50", fontface="italic", hjust=1) +
  annotate("text", x=log10(4500000), y=max(y_pos)+2.5,
           label="coverage", size=2.5, color="grey50", fontface="italic", hjust=1) +

  theme_minimal(base_size=11) +
  theme(
    plot.title=element_text(face="bold", size=15, margin=margin(b=2)),
    plot.subtitle=element_text(size=9, color="grey40", margin=margin(b=12)),
    panel.grid.major.y=element_blank(),
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_line(color="grey93"),
    axis.text.y=element_blank(),
    axis.text.x=element_text(size=9.5),
    legend.position="bottom",
    legend.box="horizontal",
    legend.text=element_text(size=9.5),
    legend.margin=margin(t=8),
    plot.margin=margin(10,20,10,5))

ggsave(file.path(outdir, "figure_nested_circles.png"), p, width=12, height=16, dpi=300)
ggsave(file.path(outdir, "figure_nested_circles.pdf"), p, width=12, height=16)
cat("Figure 5 saved.\n")
