# Figure 4: Global lineage-through-time
# Requires access to tree files — set TREE_BASE below

library(ape)
library(ggplot2)
library(scales)

script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
script_dir <- normalizePath(dirname(script_path), mustWork=FALSE)
TREE_BASE <- normalizePath(file.path(script_dir, "..", "..", "data"), mustWork=FALSE)
outdir <- normalizePath(file.path(script_dir, "..", "..", "paper", "figures"), mustWork=FALSE)
data_dir <- normalizePath(file.path(script_dir, "..", "..", "data", "processed_manuscript_inputs"), mustWork=FALSE)

dated_trees <- read.csv(file.path(data_dir, "fig4_tree_manifest.csv"), stringsAsFactors=FALSE)

time_grid <- seq(0, 600, by=0.5)
ltt_list <- list()

for (i in seq_len(nrow(dated_trees))) {
  m <- dated_trees[i, ]
  full <- file.path(TREE_BASE, m$file)
  if (!file.exists(full)) next
  tr <- tryCatch({if(isTRUE(m$nexus)) read.nexus(full) else read.tree(full)}, error=function(e) NULL)
  if (is.null(tr)) next
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  bt <- tryCatch(sort(branching.times(tr), decreasing=TRUE), error=function(e) NULL)
  if (is.null(bt) || max(bt) < 1) next
  crown <- max(bt)
  ltt_list[[m$name]] <- sapply(time_grid, function(t) {
    if (t > crown) return(0); sum(bt >= t) + 1
  })
  message(sprintf("OK: %-20s %6d tips, crown=%4.0f Ma", m$name, Ntip(tr), crown))
}

present_counts <- sapply(ltt_list, function(x) x[1])

# Order by crown age: oldest at bottom of stack, youngest on top
crown_ages <- sapply(names(ltt_list), function(nm) {
  vals <- ltt_list[[nm]]
  max(which(vals > 0)) * 0.5  # time_grid step = 0.5
})
size_order <- names(sort(crown_ages, decreasing=FALSE))  # youngest first = top of stack

# Manual stacking for log-scale compatibility
ribbon_df <- data.frame()
cumsum_so_far <- rep(0, length(time_grid))
for (nm in rev(size_order)) {
  ymin <- cumsum_so_far
  ymax <- cumsum_so_far + ltt_list[[nm]]
  cumsum_so_far <- ymax
  ribbon_df <- rbind(ribbon_df, data.frame(
    time=-time_grid, ymin=ymin, ymax=ymax, name=nm, stringsAsFactors=FALSE))
}
ribbon_df$name <- factor(ribbon_df$name, levels=size_order)
ribbon_df$ymin[ribbon_df$ymin == 0] <- 0.5
ribbon_df$ymax[ribbon_df$ymax < 0.5] <- 0.5

pal <- c("Seed plants"="#2ca02c","Fish"="#1f77b4","Birds"="#17becf",
  "Squamates"="#aec7e8","Amphibians"="#98df8a","Ferns"="#006d2c",
  "Mammals"="#ff7f0e","Bees"="#ffbb78","Neotrop. fish"="#9edae5",
  "Butterflies"="#d62728","Spiders"="#c49c94","Orchids"="#bcbd22",
  "Sharks"="#7f7f7f","Solanaceae"="#dbdb8d","Cacti"="#e377c2",
  "Salamanders"="#c5b0d5","Parrots"="#f7b6d2","Primates"="#ff9896",
  "Carnivora"="#c7c7c7","Corals"="#9467bd","Sponges"="#8c564b",
  "Cephalopods"="#e7969c","Echinoderms"="#843c39","Conifers"="#393b79")

periods <- read.csv(file.path(data_dir, "fig4_geological_periods.csv"), stringsAsFactors=FALSE)
periods$mid <- (periods$start + periods$end) / 2

ymax_plot <- max(cumsum_so_far) * 1.8

p <- ggplot() +
  geom_rect(data=periods, aes(xmin=start, xmax=end, ymin=0.5, ymax=Inf),
            fill=periods$col, alpha=0.06) +
  annotate("rect", xmin=periods$start, xmax=periods$end,
           ymin=ymax_plot*0.75, ymax=ymax_plot, fill=periods$col, alpha=0.85) +
  annotate("text", x=periods$mid, y=ymax_plot*0.87,
           label=periods$name, size=2.8, fontface="bold") +
  geom_vline(xintercept=periods$start, color="grey80", linewidth=0.15) +
  geom_ribbon(data=ribbon_df, aes(x=time, ymin=ymin, ymax=ymax, fill=name), alpha=0.9) +
  scale_fill_manual(values=pal, name="Taxonomic group", breaks=size_order) +
  scale_y_log10(labels=trans_format("log10", math_format(10^.x)),
                breaks=10^(0:5), position="right", expand=c(0,0)) +
  scale_x_continuous(limits=c(-610,5), breaks=seq(-600,0,100), expand=c(0,0)) +
  coord_cartesian(ylim=c(1, ymax_plot)) +
  labs(x="Time (Ma)", y="Number of lineages") +
  theme_minimal(base_size=12) +
  theme(legend.position="right", legend.text=element_text(size=7.5),
        legend.key.size=unit(0.35,"cm"),
        legend.title=element_text(size=9, face="bold"),
        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        plot.margin=margin(5,10,10,10)) +
  guides(fill=guide_legend(ncol=1, reverse=FALSE))

ggsave(file.path(outdir, "figure_global_ltt.png"), p, width=13, height=7, dpi=300)
ggsave(file.path(outdir, "figure_global_ltt.pdf"), p, width=13, height=7)
cat("Figure 4 saved.\n")
