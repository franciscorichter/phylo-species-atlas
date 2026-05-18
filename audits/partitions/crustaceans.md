---
partition: Crustaceans
category: Arthropods
audit_version: 1
status: in_progress
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: Appeltans 2012
  doi: 10.1016/j.cub.2012.09.036
  doi_verified: unchecked
  year_published: 2012
  cited_value: 150000

canonical_tree:
  filename: crustaceans.nwk
  provenance_group: crustaceans
  study: Wolfe JM et al. (2019)
  doi: 10.1098/rspb.2019.0079
  doi_verified: true
  crossref_match: true

methods:
  inference: bayesian
  software: [Bayesian + ML pipeline]
  data_type: UCE
  n_loci: 410
  dating_method: "410-loci anchored hybrid enrichment; 94 species across 58 decapod families."
  support_metric: posterior_probability
  divergence_ci_available: true
  note: "Decapoda-only; not representative of all crustacean diversity. 95 tips for 77k described crustaceans = 0.1% coverage."

shipped_subclades: []
candidate_subclades:
  - taxon: "Whole-tree replacement — Pancrustacea phylogenomic tree"
    rationale: |
      Wolfe 2019 ships under 'crustaceans' but covers only Decapoda
      (94 of ~9,000 decapod species; ~0.1% of all crustaceans).
      Bernot 2023 (MBE) is a Pancrustacea-wide phylogenomic synthesis:
      105 taxa across 30 of 57 crustacean orders, 576 protein-coding
      genes / 121,508 AA positions, 90 transcriptomes + 15 genomes.
      Covers Branchiopoda, Copepoda, Malacostraca, Thecostraca, Ostracoda,
      Remipedia, Cephalocarida and more — broadly representative of the
      partition.
    candidate_papers:
      - doi: 10.1093/molbev/msad175
        year: 2023
        study: Bernot JP, Owen CL, Wolfe JM, Meland K, Olesen J, Crandall KA — Major Revisions in Pancrustacean Phylogeny and Evidence of Sensitivity to Taxon Sampling
        ntips: 105
        url: https://academic.oup.com/mbe/article/40/8/msad175/7232083
        data_doi: 10.5061/dryad.dr7sqvb2h
        data_files:
          - treefiles.zip (81 KB)
          - matrices.zip (19 MB)
          - orthologs.zip (16 MB)
          - divergence_time_files.zip (38 KB)
    download_target: dryad
    priority: high
    status: proposed
    notes: |
      Dryad deposit confirmed: doi:10.5061/dryad.dr7sqvb2h. Tree files
      (81 KB ZIP) reachable via browser at
      https://datadryad.org/dataset/doi:10.5061/dryad.dr7sqvb2h but
      blocked by Cloudflare anti-bot when fetched programmatically.
      Manual download → unzip → standardize → ship. Even at
      order/family resolution (105 taxa vs Wolfe 2019's 94 decapods),
      this is a meaningful step up in breadth across crustacean classes.

resolutions:
  - issue: "Shipped tree (Wolfe 2019, 95 tips) is decapod-only — drastically under-represents the Crustacea partition (~77k described species)."
    resolved_by: "Flagged as high-priority replacement candidate. Status kept as in_progress (not verified) until a broader crustacean tree is identified or the partition is restructured."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: in_progress
---

# Crustaceans — partition audit

**Status: in_progress.** Wolfe 2019 (Proc Roy Soc B, DOI-verified) ships as
the crustacean canonical but covers only Decapoda (95 tips). This is a
**partition-tree mismatch**: a decapod-only tree shouldn't represent the full
~77k-species Crustacea partition.

Appeltans 2012 estimate (100k-250k bilateral range) appears to be paper-
published (from the World Register of Marine Species effort). High-priority
improvement candidate identified: **Bernot 2023 (MBE)** — Pancrustacean
phylogenomic synthesis with 105 taxa across 30 of 57 crustacean orders,
576 protein-coding genes. Significantly broader than Wolfe 2019's decapod
focus. Status remains in_progress until the Bernot tree is fetched,
standardized, and shipped.
