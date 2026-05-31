# Figure 3 (legacy script name fig2_crown_age.R): Crown age vs species richness
# Points coloured by ARCHIVAL UNCERTAINTY class (dating_code):
#   beast_mcc_hpd        = blue       (BEAST MCC with archived 95% HPDs)
#   posterior_sample     = green      (posterior tree distribution available)
#   bootstrap_distribution = orange   (bootstrap / treePL distribution available)
#   point_estimate       = grey       (single dated tree, no uncertainty)
#   condamine_family     = light grey (Condamine 2019 family backbone, no per-tree HPD)
# Horizontal bars retain the 7 recoverable 95% crown-age intervals from fig2_crown_age_hpds.csv.

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
rework_dir <- normalizePath(file.path(script_dir, "..", "..", "figure2-3-rework-2026-05-31"), mustWork=FALSE)

manifest <- read.csv(file.path(data_dir, "fig2_tree_manifest.csv"), stringsAsFactors=FALSE)
condamine_manifest <- manifest[manifest$source_kind == "directory", ]
tree_manifest <- manifest[manifest$source_kind == "file", ]
fallback_rows <- tree_manifest[!is.na(tree_manifest$crown_age_fallback), c("name", "crown_age_fallback")]
fallbacks <- setNames(fallback_rows$crown_age_fallback, fallback_rows$name)

# --- Per-tree dating_code mapping --------------------------------------------
# Stage 1 output: partition_slug -> dating_code. Tree manifest uses dataset
# names. Build dataset -> partition_slug bridge by mapping each non-Condamine
# source tree to the partition it canonically populates.
dating_csv <- file.path(rework_dir, "dating_code_per_tree.csv")
dating <- read.csv(dating_csv, stringsAsFactors=FALSE)

# dataset (Source author/year label in tree_manifest) -> partition_slug
dataset_to_slug <- c(
  "Smith & Brown 2018"          = "seed_plants",
  "Upham 2019"                  = "mammals",
  "Rabosky 2018"                = "fish",
  "Jetz & Pyron 2018"           = "amphibians",
  "Portik 2023"                 = "amphibians",
  "McTavish 2025"               = "birds",
  "Jetz 2012"                   = "birds",
  "Tonini 2016"                 = "squamates",
  "Kawahara 2023"               = "insects",        # butterflies subtree -> insects partition
  "Nitta 2022"                  = "ferns",
  "Vaga 2025"                   = "cnidaria",       # scleractinian corals -> cnidaria partition
  "Stein 2018"                  = "sharks",
  "Thomson 2021"                = "turtles",
  "Arnold 2010"                 = "mammals",        # primates subtree -> mammals partition
  "Garrison 2016"               = "spiders",
  "Mongiardino Koch 2022"       = "echinoderms",
  "Lavrov 2023"                 = "sponges",
  "Stewart 2025"                = "amphibians",     # salamanders subtree -> amphibians partition
  "Tagliacollo 2024"            = "fish",           # neotropical fish subtree -> fish partition
  "Burgio 2019"                 = "birds",          # parrots subtree -> birds partition
  "Nyakatura 2012"              = "mammals",        # carnivora subtree -> mammals partition
  "Henriquez-Piskulich 2024"    = "insects",        # bees subtree -> insects partition
  "Bechteler 2023"              = "bryophytes",
  "Thompson 2023"               = "seed_plants",    # orchids subtree -> seed_plants partition
  "Sarkinen 2013"               = "seed_plants",    # solanaceae subtree -> seed_plants partition
  "Thompson 2026"               = "seed_plants",    # cacti subtree -> seed_plants partition
  "Ran 2018"                    = "seed_plants",    # pinaceae subtree -> seed_plants partition
  "Basava 2024"                 = "cephalopods",
  "McGowen 2020"                = "mammals"         # cetaceans subtree -> mammals partition
)
slug_to_code <- setNames(dating$dating_code, dating$partition_slug)

# Per-dataset override: when a source tree lives in a partition whose canonical
# tree is undated/different, the SOURCE TREE itself still has its own archival
# class. These are the explicit overrides (the tree's own dating_code, not its
# partition's). Add entries here if new subtree-in-undated-partition cases
# appear (e.g., butterflies/bees living under the undated insects partition).
dataset_dating_override <- c(
  "Kawahara 2023"            = "beast_mcc_hpd",          # butterflies, BEAST MCC with HPDs
  "Henriquez-Piskulich 2024" = "bootstrap_distribution"  # bees, 1001 chronograms
)

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

# Assign dating_code per row
dat$partition_slug <- ifelse(dat$is_condamine, "condamine_family",
                             dataset_to_slug[dat$dataset])
dat$dating_code <- ifelse(dat$is_condamine, "condamine_family",
                          slug_to_code[dat$partition_slug])
# Apply per-dataset overrides (subtree-in-undated-partition cases etc.)
ovr <- dataset_dating_override[dat$dataset]
has_ovr <- !is.na(ovr)
dat$dating_code[has_ovr] <- ovr[has_ovr]
# Anything still 'undated' or NA at this point is a source tree we lack
# explicit archival info for -> classify as point_estimate.
needs_default <- is.na(dat$dating_code) | dat$dating_code == "undated"
if (any(needs_default)) {
  warning("Source trees defaulted to point_estimate (no explicit dating_code): ",
          paste(unique(dat$dataset[needs_default]), collapse=", "))
  dat$dating_code[needs_default] <- "point_estimate"
}

# Merge crown-age 95% intervals where recoverable
hpd_path <- file.path(data_dir, "fig2_crown_age_hpds.csv")
dat$crown_age_low  <- NA_real_
dat$crown_age_high <- NA_real_
dat$hpd_method     <- "none"
if (file.exists(hpd_path)) {
  hpd <- read.csv(hpd_path, stringsAsFactors = FALSE)
  hpd <- hpd[!is.na(hpd$crown_age_low) & !is.na(hpd$crown_age_high) &
             (hpd$crown_age_high - hpd$crown_age_low) > 0.001, ]
  m <- match(dat$dataset, hpd$dataset)
  has <- !is.na(m)
  dat$crown_age_low[has]  <- hpd$crown_age_low[m[has]]
  dat$crown_age_high[has] <- hpd$crown_age_high[m[has]]
  dat$hpd_method[has]     <- hpd$hpd_method[m[has]]
}
dat$has_ci <- !is.na(dat$crown_age_low) & !is.na(dat$crown_age_high)

# --- Aesthetic encoding ------------------------------------------------------
dating_levels <- c("beast_mcc_hpd", "posterior_sample", "bootstrap_distribution",
                   "point_estimate", "condamine_family")
dating_labels <- c(
  "beast_mcc_hpd"          = "BEAST MCC with HPD",
  "posterior_sample"       = "Posterior sample",
  "bootstrap_distribution" = "Bootstrap distribution",
  "point_estimate"         = "Point estimate (no uncertainty)",
  "condamine_family"       = "Condamine 2019 family backbone"
)
dating_colors <- c(
  "beast_mcc_hpd"          = "#1F77B4",  # blue
  "posterior_sample"       = "#2CA02C",  # green
  "bootstrap_distribution" = "#FF7F0E",  # orange
  "point_estimate"         = "#7F7F7F",  # grey
  "condamine_family"       = "#CCCCCC"   # light grey
)

dat$dating_code <- factor(dat$dating_code, levels=dating_levels)

shape_map <- c("Family"=16,"Order"=17,"Class"=15,"Division"=18,"Superfamily"=8,"Phylum"=4)

ci_dat <- dat[!dat$is_condamine & dat$has_ci, ]
n_with_ci <- nrow(ci_dat)
n_no_ci   <- sum(!dat$is_condamine & !dat$has_ci)
n_total_nc <- n_with_ci + n_no_ci

# Tally per-class for caption
class_counts <- table(dat$dating_code[!dat$is_condamine])

p <- ggplot() +
  # Condamine backbone first (light grey, small markers)
  geom_point(data=dat[dat$is_condamine,],
             aes(x=crown_age, y=ntips, colour=dating_code),
             alpha=0.45, size=1.6, shape=16) +
  # Horizontal 95% HPD bars (colour by dating_code so bar matches its point)
  geom_errorbarh(data=ci_dat,
                 aes(y=ntips, xmin=crown_age_low, xmax=crown_age_high, colour=dating_code),
                 height=0, linewidth=0.8, alpha=0.75) +
  # Non-Condamine source trees
  geom_point(data=dat[!dat$is_condamine,],
             aes(x=crown_age, y=ntips, colour=dating_code, shape=level),
             size=3.8, alpha=0.92, stroke=0.6) +
  geom_text_repel(data=dat[!dat$is_condamine,],
    aes(x=crown_age, y=ntips, label=name), size=2.8, max.overlaps=30,
    segment.color="grey60", segment.size=0.3, box.padding=0.4) +
  scale_y_log10(labels=comma, breaks=c(10,100,1000,10000,100000,300000)) +
  scale_colour_manual(values=dating_colors, labels=dating_labels, drop=FALSE,
                      name="Archival uncertainty") +
  scale_shape_manual(values=shape_map, name="Taxonomic scope") +
  labs(title=sprintf("Crown age vs. species richness across %d dated phylogenetic trees", nrow(dat)),
       subtitle=sprintf(
         "Colour = archival uncertainty class. Horizontal bars show recoverable 95%% crown-age intervals (n = %d of %d non-Condamine source trees).",
         n_with_ci, n_total_nc),
       x="Crown age (Ma)", y="Number of species (log scale)",
       caption=sprintf(
         "Non-Condamine source trees: BEAST MCC+HPD n=%d  |  posterior sample n=%d  |  bootstrap n=%d  |  point estimate n=%d.  Condamine 2019: %d family-level trees.",
         class_counts["beast_mcc_hpd"], class_counts["posterior_sample"],
         class_counts["bootstrap_distribution"], class_counts["point_estimate"],
         sum(dat$is_condamine))) +
  theme_bw(base_size=12) +
  theme(legend.position="right", legend.text=element_text(size=8),
        legend.title=element_text(size=9, face="bold"),
        legend.key.size=unit(0.5,"cm"),
        plot.title=element_text(size=13, face="bold"),
        plot.subtitle=element_text(size=9, face="italic", colour="grey30"),
        plot.caption=element_text(size=8, colour="grey25", hjust=0),
        panel.grid.minor=element_blank()) +
  guides(colour=guide_legend(ncol=1, order=1, override.aes=list(size=3, shape=16, alpha=1)),
         shape=guide_legend(ncol=1, order=2, override.aes=list(size=2.8, colour="grey30")))

ggsave(file.path(outdir, "crown_age_vs_species.png"), p, width=12, height=8, dpi=300)
ggsave(file.path(outdir, "crown_age_vs_species.pdf"), p, width=12, height=8)
cat("Figure 3 (crown_age_vs_species) saved with dating_code colour encoding.\n")
cat(sprintf("Class counts (non-Condamine): beast_mcc_hpd=%d, posterior=%d, bootstrap=%d, point=%d; condamine n=%d.\n",
            class_counts["beast_mcc_hpd"], class_counts["posterior_sample"],
            class_counts["bootstrap_distribution"], class_counts["point_estimate"],
            sum(dat$is_condamine)))
