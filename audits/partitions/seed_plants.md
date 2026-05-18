---
partition: Seed plants
category: Plants
audit_version: 1
status: verified
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: Nic Lughadha 2016
  doi: 10.11646/phytotaxa.261.3.1
  doi_verified: unchecked
  crossref_match: unchecked
  year_published: 2016
  cited_value: 450000
  newer_authoritative_known: true
  newer_candidate:
    doi: null
    year: 2025
    study: World Flora Online
    cited_value: 383671
    recommendation: watch
    notes: "WFO tracks the canonical described count; Nic Lughadha 2016 published the bilateral estimate."

canonical_tree:
  filename: seed_plants.nwk
  provenance_group: seed_plants
  study: Smith SA, Brown JW (2018)
  doi: 10.1002/ajb2.1019
  doi_verified: true
  crossref_match: true
  newer_authoritative_known: true
  newer_candidate:
    doi: 10.1038/s41586-024-07324-0
    year: 2024
    study: Zuntini AR et al. — Phylogenomics and the rise of the angiosperms
    cited_value: 9505
    recommendation: watch
    notes: "Zuntini 2024 Nature is a genomic-scale phylogenomic synthesis (353 nuclear genes, 9505 species). Different methodology (phylogenomic vs Smith 2018 supermatrix); both remain valid canonical candidates."

methods:
  inference: ml
  software: [RAxML, treePL]
  data_type: supermatrix
  n_loci: 7
  dating_method: "7-gene ML supermatrix, RAxML + treePL dating. Tree built in stages: ~80k seed plants with sequence data, extended to ~356k via taxonomic placement."
  support_metric: bootstrap
  divergence_ci_available: false
  posterior_n: null
  note: "6 tree variants available (ALLMB, ALLOTB, GBMB, GBOTB, all-mol-backbone, all-otb-backbone). Atlas ships the ALLMB variant."

shipped_subclades:
  - filename: grasses.nwk
    study: GPWG III (2025)
    taxon: Poaceae
    ntips: 1153
  - filename: orchids.nwk
    study: Thompson JB et al. (2023)
    taxon: Orchidaceae
    ntips: 1475
  - filename: asteraceae.nwk
    study: Mandel JR et al. (2019)
    taxon: Asteraceae
    ntips: 256
  - filename: solanaceae.nwk
    study: Sarkinen et al. (2013)
    taxon: Solanaceae
    ntips: 1076
  - filename: cacti.nwk
    study: Thompson J et al. (2026)
    taxon: Cactaceae
    ntips: 1063
  - filename: eucalyptus.nwk
    study: Crisp MD et al. (2024)
    taxon: Eucalyptus
    ntips: 399

candidate_subclades:
  - taxon: "Whole-tree refresh — Zuntini 2024 phylogenomic synthesis"
    rationale: "9,505 angiosperm species at 353 nuclear genes (Angiosperms353 probe set). Phylogenomic resolution far exceeds Smith 2018's 7 plastid genes; consider as the new canonical."
    candidate_papers:
      - doi: 10.1038/s41586-024-07324-0
        year: 2024
        study: Zuntini AR et al. — Phylogenomics and the rise of the angiosperms
        ntips: 9505
    download_target: dryad
    priority: high
    status: proposed

resolutions:
  - issue: "356k tips include ~80k with direct sequence data + ~276k added via taxonomic placement. The latter have very low species-level confidence."
    resolved_by: "Documented in methods.note. Users should consult Smith 2018 §Methods for placement uncertainty."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified
---

# Seed plants — partition audit

**Status: verified.** Smith 2018 ALLMB tree (DOI-verified). Nic Lughadha 2016 estimate
is a paper-published bilateral range (400k–500k). Zuntini 2024 Nature is a
strong candidate for canonical-tree upgrade (phylogenomic, 353 nuclear genes).

Six sub-clades shipped: grasses, orchids, asteraceae, solanaceae, cacti,
eucalyptus.
