# Figure 2: Crown age vs species richness
# Requires access to tree files — set TREE_BASE below

library(ape)
library(ggplot2)
library(scales)
library(ggrepel)

script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
script_dir <- normalizePath(dirname(script_path), mustWork=FALSE)
TREE_BASE <- normalizePath(file.path(script_dir, "..", "..", "..", "Data"), mustWork=FALSE)
outdir <- normalizePath(file.path(script_dir, "..", "..", "figures"), mustWork=FALSE)
data_dir <- normalizePath(file.path(script_dir, "..", "..", "inputs"), mustWork=FALSE)

manifest <- read.csv(file.path(data_dir, "fig2_tree_manifest.csv"), stringsAsFactors=FALSE)
condamine_manifest <- manifest[manifest$source_kind == "directory", ]
tree_manifest <- manifest[manifest$source_kind == "file", ]
fallback_rows <- tree_manifest[!is.na(tree_manifest$crown_age_fallback), c("name", "crown_age_fallback")]
fallbacks <- setNames(fallback_rows$crown_age_fallback, fallback_rows$name)

# Read Condamine trees
cond_dir <- file.path(TREE_BASE, condamine_manifest$file[1])
cond_rows <- data.frame()
if (dir.exists(cond_dir)) {
  for (f in list.files(cond_dir, pattern=condamine_manifest$pattern[1], recursive=TRUE, full.names=TRUE)) {
    tr <- tryCatch(read.tree(f), error=function(e) NULL)
    if (is.null(tr)) next
    ca <- tryCatch(max(branching.times(tr)), error=function(e) NA)
    if (is.na(ca) || ca < 1) next
    cond_rows <- rbind(cond_rows, data.frame(
      name=basename(tools::file_path_sans_ext(f)), dataset="Condamine 2019",
      level="Family", ntips=Ntip(tr), crown_age=ca, stringsAsFactors=FALSE))
  }
}

# Read other trees
other_rows <- data.frame()
for (i in seq_len(nrow(tree_manifest))) {
  m <- tree_manifest[i, ]
  full <- file.path(TREE_BASE, m$file)
  if (!file.exists(full)) next
  tr <- tryCatch({if(isTRUE(m$nexus)) read.nexus(full) else read.tree(full)}, error=function(e) NULL)
  if (is.null(tr)) next
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  ca <- tryCatch(max(branching.times(tr)), error=function(e) NA)
  if (is.na(ca) || ca < 1) {
    if (m$name %in% names(fallbacks)) ca <- fallbacks[m$name] else next
  }
  other_rows <- rbind(other_rows, data.frame(
    name=m$name, dataset=m$dataset, level=m$level,
    ntips=Ntip(tr), crown_age=ca, stringsAsFactors=FALSE))
}

dat <- rbind(cond_rows, other_rows)
dat$is_condamine <- dat$dataset == "Condamine 2019"

# Merge crown-age 95% intervals where recoverable
hpd_path <- file.path(data_dir, "fig2_crown_age_hpds.csv")
dat$crown_age_low  <- NA_real_
dat$crown_age_high <- NA_real_
dat$hpd_method     <- "none"
if (file.exists(hpd_path)) {
  hpd <- read.csv(hpd_path, stringsAsFactors = FALSE)
  hpd <- hpd[!is.na(hpd$crown_age_low) & !is.na(hpd$crown_age_high) &
             (hpd$crown_age_high - hpd$crown_age_low) > 0.001, ]
  # Match by dataset (not name). One row per dataset in hpd table; one or more per
  # dataset in dat (most have a single tree, but some datasets have multiple
  # representative entries). Apply the same interval to all rows of the dataset.
  m <- match(dat$dataset, hpd$dataset)
  has <- !is.na(m)
  dat$crown_age_low[has]  <- hpd$crown_age_low[m[has]]
  dat$crown_age_high[has] <- hpd$crown_age_high[m[has]]
  dat$hpd_method[has]     <- hpd$hpd_method[m[has]]
}
dat$has_ci <- !is.na(dat$crown_age_low) & !is.na(dat$crown_age_high)

all_datasets <- unique(dat$dataset)
base_colors <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
  "#FF7F00","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB",
  "#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#1B9E77","#D95F02",
  "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#8DD3C7",
  "#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462")
pal <- setNames(colorRampPalette(base_colors)(length(all_datasets)), all_datasets)
# Ensure Condamine stays red
pal["Condamine 2019"] <- "#E41A1C"

shape_map <- c("Family"=16,"Order"=17,"Class"=15,"Division"=18,"Superfamily"=8,"Phylum"=4)

ci_dat <- dat[!dat$is_condamine & dat$has_ci, ]
n_with_ci <- nrow(ci_dat)
n_no_ci   <- sum(!dat$is_condamine & !dat$has_ci)

p <- ggplot() +
  geom_point(data=dat[dat$is_condamine,],
             aes(x=crown_age, y=ntips, colour=dataset),
             alpha=0.35, size=1.8, shape=16) +
  geom_errorbarh(data=ci_dat,
                 aes(y=ntips, xmin=crown_age_low, xmax=crown_age_high, colour=dataset),
                 height=0, linewidth=0.7, alpha=0.7) +
  geom_point(data=dat[!dat$is_condamine,],
             aes(x=crown_age, y=ntips, colour=dataset, shape=level),
             size=3.5, alpha=0.85) +
  geom_text_repel(data=dat[!dat$is_condamine,],
    aes(x=crown_age, y=ntips, label=name), size=2.8, max.overlaps=30,
    segment.color="grey60", segment.size=0.3, box.padding=0.4) +
  scale_y_log10(labels=comma, breaks=c(10,100,1000,10000,100000,300000)) +
  scale_colour_manual(values=pal, name="Source") +
  scale_shape_manual(values=shape_map, name="Taxonomic scope") +
  labs(title=sprintf("Crown age vs. species richness across %d dated phylogenetic trees", nrow(dat)),
       subtitle=sprintf("Horizontal bars show recoverable 95%% intervals on crown age (n = %d of %d non-Condamine source trees); remaining %d trees lack archived uncertainty.",
                        n_with_ci, n_with_ci + n_no_ci, n_no_ci),
       x="Crown age (Ma)", y="Number of species (log scale)") +
  theme_bw(base_size=12) +
  theme(legend.position="right", legend.text=element_text(size=7),
        legend.title=element_text(size=9, face="bold"),
        legend.key.size=unit(0.4,"cm"),
        plot.title=element_text(size=13, face="bold"),
        plot.subtitle=element_text(size=9, face="italic", colour="grey30"),
        panel.grid.minor=element_blank()) +
  guides(colour=guide_legend(ncol=1, override.aes=list(size=2.5)),
         shape=guide_legend(override.aes=list(size=2.5)))

ggsave(file.path(outdir, "crown_age_vs_species.png"), p, width=12, height=8, dpi=300)
ggsave(file.path(outdir, "crown_age_vs_species.pdf"), p, width=12, height=8)
cat("Figure 2 saved.\n")
