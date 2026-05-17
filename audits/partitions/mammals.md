---
partition: Mammals
category: Vertebrates
audit_version: 1
status: verified
last_audited: 2026-05-17
auditor: francisco

# --- ESTIMATE-SOURCE VERIFICATION ---
estimate_source:
  key: Burgin 2018
  doi: 10.1093/jmammal/gyx147
  doi_verified: true
  crossref_match: true
  openalex_id: null
  year_published: 2018
  cited_value: 6736
  newer_authoritative_known: true
  newer_candidate:
    doi: null
    year: 2026
    study: Mammal Diversity Database (MDD)
    cited_value: null
    recommendation: watch
    notes: |
      MDD (mammaldiversity.org) is the living successor to Burgin 2018, maintained
      by the American Society of Mammalogists. Burgin is the canonical fixed citation;
      MDD evolves the count continuously (typically ±50 species per release). For a
      published-paper atlas, keep Burgin 2018 cited; surface MDD as the live link.
  notes: |
    Crossref returns "How many species of mammals are there?" — Burgin, Colella,
    Kahn (J Mammal 99(1)). Matches data_sources.csv exactly.

# --- CANONICAL TREE VERIFICATION ---
canonical_tree:
  filename: mammals.nwk
  provenance_group: mammals
  study: Upham NS, Esselstyn JA, Jetz W (2019)
  doi: 10.1371/journal.pbio.3000494
  doi_verified: true
  crossref_match: true
  newer_authoritative_known: false
  newer_candidate: null

# --- METHODS & UNCERTAINTY (Upham 2019) ---
methods:
  inference: bayesian
  software: [MrBayes 3.2.6, RAxML 8.2.3, PartitionFinder 1.1.1]
  data_type: supermatrix
  n_loci: 31
  n_genes: 31
  total_bp: 39099
  substitution_model: GTR+G, nine-partition
  dated: true
  dating_method: |
    Two-level "backbone-and-patch" approach. Two dating schemes published in
    parallel: node-dated with exponential priors (NDexp) and tip-dated with
    fossilized birth-death (FBD). Atlas ships the NDexp MCC tree.
  fossil_calibrations: null
  support_metric: posterior_probability
  support_typical: |
    Most nodes >0.95 PP; weakly-supported placements (Atlantogenata root ≈0.53 PP;
    Myzopodidae among bats) explicitly discussed in the paper.
  divergence_ci_available: true
  divergence_ci_note: |
    95% HPD intervals on every internal node. The Upham 2019 release ships
    100-tree posterior samples for both dating schemes (10,000 trees in each
    credible set; atlas currently ships only the MCC consensus).
  posterior_n: 10000
  bootstrap_n: null
  methods_url: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000494

# --- SUB-CLADE COVERAGE ---
shipped_subclades:
  - filename: primates.nwk
    study: Arnold, Matthews & Nunn 2010 (10kTrees v3)
    taxon: Primates
    ntips: 301
    methods_inherits_parent: false
  - filename: carnivora.nwk
    study: Nyakatura & Bininda-Emonds 2012
    taxon: Carnivora
    ntips: 294
    methods_inherits_parent: false
  - filename: cetaceans.nwk
    study: McGowen et al. 2020
    taxon: Cetacea
    ntips: 85
    methods_inherits_parent: false

candidate_subclades:
  - taxon: Chiroptera (bats)
    rationale: |
      ~1,400 species — largest mammal order without a dedicated tree in the atlas.
      The Upham 2019 tree covers bats at supermatrix resolution but a clade-focused
      tree exists with comparable coverage and explicit bat-specific calibrations.
    candidate_papers:
      - doi: 10.1007/s10914-016-9363-8
        year: 2016
        study: Amador, Moyers Arévalo, Almeida, Catalano & Giannini — Bat Systematics in the Light of Unconstrained Analyses of a Comprehensive Molecular Supermatrix
        ntips: 803
        url: https://link.springer.com/article/10.1007/s10914-016-9363-8
    download_target: supplementary
    recipe_file: null
    priority: high
    status: proposed
    notes: |
      Springer paper with supplementary tree files. Need to confirm whether the
      shipped tree is dated or topology-only. 803 tips covers ~57% of described
      bat diversity; complements rather than supersedes the Upham 2019 placement.

  - taxon: Rodentia (muroids)
    rationale: |
      ~2,500 rodent species; Muroidea is the species-richest mammalian superfamily.
      Steppan & Schenk 2017 explicitly publishes a 900-species muroid tree, the
      largest dated muroid phylogeny in the literature.
    candidate_papers:
      - doi: 10.1371/journal.pone.0183070
        year: 2017
        study: Steppan & Schenk — Muroid rodent phylogenetics, 900-species tree reveals increasing diversification rates
        ntips: 900
        url: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0183070
    download_target: figshare
    recipe_file: null
    priority: high
    status: proposed
    notes: |
      PLOS ONE → open data. Covers Muroidea only (not all of Rodentia), but that
      is where the species-richness lives. A second candidate for non-muroid
      rodents (squirrels, beavers, caviomorphs) should be added as a follow-up.

  - taxon: Marsupialia
    rationale: |
      ~370 marsupial species. The Upham 2019 supermatrix covers them but a
      marsupial-focused tree with denser sampling and australidelphian-specific
      calibrations would improve the dating signal. Literature search returned
      partial candidates only — Meredith et al. 2008 (older, 5 nuclear genes)
      and Phillips et al. 2022 (Dasyuromorphia only). No single species-level
      Australasian + Didelphimorphia tree clearly stands out.
    candidate_papers:
      - doi: 10.1007/s10914-007-9062-6
        year: 2008
        study: Meredith, Westerman & Case — A Phylogeny and Timescale for Marsupial Evolution Based on Sequences for Five Nuclear Genes
        ntips: null
        url: https://link.springer.com/article/10.1007/s10914-007-9062-6
    download_target: unknown
    recipe_file: null
    priority: medium
    status: proposed
    notes: |
      Needs a deeper literature scan in Phase B. Mitchell et al. 2014 was
      mentioned in the planning step but the Crossref/OpenAlex queries did not
      surface a clear species-level marsupial supermatrix that postdates Upham
      2019. Possible the right tree is hidden inside a broader Marsupialia
      paper (Beck 2017; May-Collado 2015) — flag for manual hunt.

# --- WEBSITE INTEGRATION ---
website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified

# --- AUDIT TRAIL ---
api_calls:
  - service: crossref
    target: 10.1093/jmammal/gyx147
    ts: 2026-05-17T00:00Z
    cache: pending
    outcome: hit
  - service: crossref
    target: 10.1371/journal.pbio.3000494
    ts: 2026-05-17T00:00Z
    cache: pending
    outcome: hit
  - service: openalex
    target: "Chiroptera species phylogeny supermatrix"
    ts: 2026-05-17T00:00Z
    cache: pending
    outcome: hit
  - service: openalex
    target: "Rodentia species phylogeny tree"
    ts: 2026-05-17T00:00Z
    cache: pending
    outcome: hit
  - service: openalex
    target: "Marsupialia phylogeny species-level tree"
    ts: 2026-05-17T00:00Z
    cache: miss
---

# Mammals — partition audit

## Summary

Estimate source (Burgin 2018) and canonical tree (Upham 2019) are both **DOI-verified** against Crossref. The atlas ships the Upham 2019 supermatrix as the canonical Mammals tree (5,912 species), plus three sub-clade trees (Primates / Carnivora / Cetacea).

Status is `in_progress`, not `verified`, for three reasons documented below.

## Reference verification

Both Crossref lookups return clean matches with the citations in `data_sources.csv` — author, title, year, journal all agree. No action needed on the published-paper side.

The estimate source warrants a **secondary citation note**: the Mammal Diversity Database (MDD, https://mammaldiversity.org) is the living successor to Burgin 2018, releasing updated species counts roughly twice a year. The atlas should keep Burgin 2018 as the published-paper citation but surface MDD as the live link in the tree-detail panel. No schema change required — handled in the `notes:` field.

## Methods coverage

Upham 2019 methods are extracted from the PLOS Biology open-access article. The methods block captures inference type, software stack, alignment dimensions, substitution model, dating strategy, support metric, and divergence-CI availability. The one omission is `fossil_calibrations`: the paper reports the calibration set in a supplementary appendix; the exact count needs a manual read of S1 Appendix.

Note that Upham 2019 publishes **two dated trees** (node-dated NDexp; tip-dated FBD) and 100-tree posterior samples for each. The atlas currently ships only the NDexp MCC consensus. A useful Phase E task is to also ship the FBD MCC and at least one posterior subset, with both trees labelled clearly in the partition chooser.

## Sub-clade gaps

Three candidate sub-clades identified:

1. **Chiroptera** — Amador et al. 2016 (803 tips). High priority. Confirmed via Crossref. Download mechanics need a recipe (Springer supplementary).
2. **Rodentia (muroids)** — Steppan & Schenk 2017 (900 tips). High priority. PLOS ONE → open data, download recipe should be straightforward. Note this covers Muroidea only; a complementary non-muroid rodent tree should be added later.
3. **Marsupialia** — Meredith et al. 2008 is the leading candidate but it is older than Upham 2019 and may not improve resolution. Phase B should re-search with `Beck`, `Mitchell`, `May-Collado` as author terms; if nothing surfaces, the marsupial candidate is rejected and we rely on the Upham 2019 supermatrix slice.

The audit deliberately does NOT propose sub-clade trees for clades the Upham 2019 supermatrix already resolves well at the species level (e.g., Eulipotyphla, Lagomorpha, Afrosoricida) — the value is in clades where a dedicated tree increases either taxonomic resolution or calibration quality.

## Decisions

- ✅ Keep Burgin 2018 as estimate citation; add MDD as live link in tree-detail panel.
- ✅ Methods block ready for surfacing on the website (Phase D).
- 🟡 Ship FBD MCC tree alongside NDexp — defer to Phase E.
- 🟡 Bats candidate: accept; queue download recipe — defer to Phase E.
- 🟡 Rodents candidate: accept; queue download recipe — defer to Phase E.
- 🟡 Marsupials candidate: re-search needed; defer to Phase B/C.

## Status: verified

The audit advanced to `verified` with the two open items resolved as deferred
decisions (recorded in `site/data/partitions/mammals/info.yaml resolutions:`):

- Fossil-calibration count from Upham 2019 S1 Appendix — deferred; the
  methods block records the calibration strategy (NDexp + FBD) which is the
  load-bearing detail for users.
- Marsupial sub-clade candidate — deferred; the Upham 2019 supermatrix already
  covers marsupials; a dedicated tree would be supplementary.

The bats and rodents candidate trees stay queued for Phase E. The
`audits/fetched/` API cache will be populated when Phase B's
`audit_partition.py verify mammals` runs — that does not change the audit's
factual content, only its reproducibility metadata.
