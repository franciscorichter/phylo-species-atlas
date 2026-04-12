# Figure 3: Dated vs undated trees by category
# Self-contained — all data embedded

library(ggplot2)

script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
script_dir <- normalizePath(dirname(script_path), mustWork=FALSE)
outdir <- normalizePath(file.path(script_dir, "..", "..", "paper", "figures"), mustWork=FALSE)
data_dir <- normalizePath(file.path(script_dir, "..", "..", "data", "processed_manuscript_inputs"), mustWork=FALSE)

dated_data <- read.csv(file.path(data_dir, "fig3_dated_undated_counts.csv"), stringsAsFactors=FALSE)
dated_data$category <- factor(dated_data$category,
  levels=c("Vertebrates","Plants","Arthropods","Microbes & protists","Other animals"))
dated_data$type <- factor(dated_data$type, levels=c("Dated","Undated"))

p <- ggplot(dated_data, aes(x=category, y=count, fill=type)) +
  geom_col() +
  scale_fill_manual(values=c("Dated"="#4DAF4A","Undated"="#E41A1C")) +
  coord_flip() +
  labs(title="Time-calibrated vs. topology-only trees",
       x=NULL, y="Number of datasets", fill="") +
  theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold"), legend.position="right")

ggsave(file.path(outdir, "figure_dated.png"), p, width=8, height=4, dpi=300)
ggsave(file.path(outdir, "figure_dated.pdf"), p, width=8, height=4)
cat("Figure 3 saved.\n")
