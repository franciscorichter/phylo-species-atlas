---
partition: Insects
category: Arthropods
audit_version: 1
status: verified
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: Stork 2018
  doi: 10.1146/annurev-ento-020117-043348
  doi_verified: true
  crossref_match: true
  year_published: 2018
  cited_value: 5500000
  newer_authoritative_known: false

canonical_tree:
  filename: insects.nwk
  provenance_group: insects
  study: Chesters D et al. (2023)
  doi: 10.1111/1755-0998.13817
  doi_verified: true
  crossref_match: true

methods:
  inference: synthesis
  software: [insectphylo pipeline]
  data_type: supertree
  dated: false
  dating_method: "Synthesis molecular phylogeny from GenBank data; continuously updated. Atlas v2 from Feb 2026."
  support_metric: source_count
  divergence_ci_available: false
  note: "Undated. 24 order-level trees also available."

shipped_subclades:
  - filename: butterflies.nwk
    study: Kawahara AY et al. (2023)
    taxon: Lepidoptera (butterflies)
    ntips: 2258
  - filename: bees.nwk
    study: Henriquez-Piskulich PA et al. (2024)
    taxon: Apoidea (bees)
    ntips: 4586
  - filename: hemiptera.nwk
    study: Johnson KP et al. (2018)
    taxon: Hemiptera
    ntips: 193
  - filename: termites.nwk
    study: Hellemans S et al. (2024)
    taxon: Isoptera (termites)
    ntips: 138
  - filename: ants.nwk
    study: Romiguier J et al. (2022)
    taxon: Formicidae (ants)
    ntips: 83

resolutions:
  - issue: "Atlas's high bound (7.8M) and low (2.6M) come from Stork 2018's REVIEW table of multiple published estimates, not Stork's own single estimate."
    resolved_by: "Documented as paper-published-range, with auditor note clarifying it's a meta-range across multiple primary estimates (Stork's own central is ~5.5M)."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified
---

# Insects — partition audit

**Status: verified.** Chesters 2023 insectphylo v2 (DOI-verified). Stork 2018
(Annual Review of Entomology, DOI-verified) provides paper-published bilateral
range, but the 2.6M-5.5M-7.8M bounds come from Stork's TABLE of multiple
estimates, not a single primary estimate.

Five sub-clades shipped: butterflies (Kawahara 2023, 2258 tips), bees
(Henríquez-Piskulich 2024, 4586 tips), hemiptera (Johnson 2018, 193 tips),
termites (Hellemans 2024, 138 tips), ants (Romiguier 2022, 83 tips). The
insects partition is the worst-described in the atlas at ~1% species
coverage of the central estimate.
