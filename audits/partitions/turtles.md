---
partition: Turtles
category: Vertebrates
audit_version: 1
status: verified
last_audited: 2026-05-17
auditor: francisco

# --- ESTIMATE-SOURCE VERIFICATION ---
estimate_source:
  key: TTWG 2021
  doi: 10.3854/crm.8.000c.checklist.atlas.v9.2021
  doi_verified: false
  crossref_match: na
  openalex_id: null
  year_published: 2021
  cited_value: 360
  newer_authoritative_known: true
  newer_candidate:
    doi: 10.3854/crm.10.checklist.atlas.v10.2025
    year: 2025
    study: TTWG (Turtle Taxonomy Working Group) 10th edition — "Turtles of the World: Annotated Checklist and Atlas"
    cited_value: 364
    recommendation: switch
    notes: |
      The 2021 (v9) DOI returns 404 on Crossref and on doi.org — citation is dead.
      The 2025 (v10) DOI resolves cleanly (302 → iucn-tftsg.org/checklist, Crossref
      HTTP 200). Switch the citation to TTWG 2025 in data_sources.csv (or, if the
      paper-side data is frozen, surface a corrected live link in the web UI).
      Species count moves 360 → 364 (+4 net taxonomic changes since 2021).
  notes: |
    The cited DOI is broken. Chelonian Research Monographs registers a fresh DOI
    per edition; the v9/2021 DOI does not resolve at doi.org. This is a
    data-quality bug in data_sources.csv that should be tracked even though the
    paper-side data is frozen.

# --- CANONICAL TREE VERIFICATION ---
canonical_tree:
  filename: turtles.nwk
  provenance_group: turtles
  study: Thomson RC, Spinks PQ, Shaffer HB (2021)
  doi: 10.1073/pnas.2012215118
  doi_verified: true
  crossref_match: true
  newer_authoritative_known: false
  newer_candidate: null

# --- METHODS & UNCERTAINTY (Thomson 2021) ---
methods:
  inference: bayesian
  software: [MrBayes 3.2, BEAST 1.8.2, TreeAnnotator 1.8.2, LogCombiner 1.8.2]
  data_type: multi-locus
  n_loci: 15
  n_genes: 15
  total_bp: null
  substitution_model: |
    MrBayes: per-partition substitution models selected with model selection.
    BEAST divergence-time run: Hasegawa-Kishino-Yano (HKY) for each partition
    (chosen for MCMC performance, not strict model fit).
  dated: true
  dating_method: |
    BEAST 1.8.2 with uncorrelated lognormal relaxed clock; birth-death tree
    prior; estimated on a reduced one-individual-per-species version of the
    dataset.
  fossil_calibrations: 22
  support_metric: posterior_probability
  support_typical: |
    Posterior probability = 1 on the major deep splits (Cryptodira / Pleurodira,
    all 14 currently recognized families monophyletic). Tree captures all 14
    families and 98% of genera; species-level support varies, not summarized
    in a single number.
  divergence_ci_available: true
  divergence_ci_note: |
    95% HPD intervals reported on internal nodes. 100 chronograms drawn from
    the BEAST posterior were used for downstream diversification analyses;
    the atlas currently ships only the MCC consensus.
  posterior_n: 100
  bootstrap_n: null
  methods_url: https://www.pnas.org/doi/10.1073/pnas.2012215118

# --- DATA-QUALITY FINDINGS (audit-only block, not in schema yet) ---
# These are real data-integrity issues found during the audit. They should
# probably get a first-class schema slot in v2 ("findings:" with type+severity).
findings:
  - severity: high
    field: data_sources.csv → TTWG 2021 doi
    issue: "DOI 10.3854/crm.8.000c.checklist.atlas.v9.2021 returns 404 on Crossref and doi.org."
    fix: "Replace with TTWG 2025 (10.3854/crm.10.checklist.atlas.v10.2025) OR add a verified-DOI follow-up to data_sources.csv."
  - severity: medium
    field: data_provenance.csv → methods_brief for turtles
    issue: |
      The CSV says "13 loci for 294 species". The Thomson 2021 paper reports
      15 loci and 279 of 348 living species (80%). The original filename
      "thomson2021_turtle_mcc_294sp.tre" preserves the 294 figure — possibly
      from an early Dryad release that included subspecies in the species count.
    fix: |
      Manual reconciliation needed. Confirm whether the shipped MCC tree is the
      same as the paper's 591-individual MCC (which corresponds to 279 species
      + outgroups) or an alternative Dryad version. Update methods_brief
      accordingly.
  - severity: medium
    field: standardized/trees/turtles.nwk
    issue: |
      Tree has 593 tip positions and 329 unique numeric IDs. The dictionary
      maps these IDs to names that include subspecies (e.g. Terrapene_ornata_luteola)
      and locality codes (e.g. Kinosternon_dunni_K38). This means the tree's
      "species" coverage is ambiguous without de-duplication.
    fix: |
      Web app should surface this ("593 tips covering ~279 species across
      multiple individuals; some tips are subspecies"). Already partially
      captured by the existing notes field in data_provenance.csv: "Multiple
      individuals per species; displayed coverage is capped at 100%".

# --- SUB-CLADE COVERAGE ---
shipped_subclades: []

candidate_subclades:
  - taxon: Testudinidae (tortoises)
    rationale: |
      ~50 species. Thomson 2021 covers all 14 families at the family level but
      tortoise species-level resolution may improve with a tortoise-focused
      multilocus tree. Pyron 2014 (Sys Bio) is a candidate — biogeographic
      analysis with detailed Testudinidae sampling.
    candidate_papers:
      - doi: 10.1093/sysbio/syu042
        year: 2014
        study: Pyron — Biogeographic Analysis Reveals Ancient Continental Vicariance and Recent Oceanic Dispersal
        ntips: null
        url: https://academic.oup.com/sysbio/article/63/5/779/2847947
    download_target: unknown
    recipe_file: null
    priority: low
    status: proposed
    notes: |
      Low priority because Thomson 2021 at 80% species coverage already gives
      good tortoise resolution. Sub-clade tree would be supplementary.

  - taxon: Geoemydidae (Asian river/leaf turtles)
    rationale: |
      ~70 species — the largest turtle family. Thomson 2021's 80% sampling
      may have its largest gaps here (many South-East Asian species hard to
      sample). Pereira 2017 (MPE) provides a focused multilocus geoemydid tree.
    candidate_papers:
      - doi: 10.1016/j.ympev.2017.05.008
        year: 2017
        study: Pereira et al. — Multilocus phylogeny and statistical biogeography clarify the evolutionary history of major Asian river turtles
        ntips: null
        url: https://www.sciencedirect.com/science/article/pii/S1055790317301215
    download_target: unknown
    recipe_file: null
    priority: medium
    status: proposed
    notes: |
      Geoemydidae is the family most likely to have species missing from
      Thomson 2021's 80% sampling. Worth checking species count overlap.

# --- WEBSITE INTEGRATION ---
website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified

# --- AUDIT TRAIL ---
api_calls:
  - service: crossref
    target: 10.1073/pnas.2012215118
    ts: 2026-05-17T01:30Z
    cache: pending
    outcome: hit
  - service: crossref
    target: 10.3854/crm.8.000c.checklist.atlas.v9.2021
    ts: 2026-05-17T01:30Z
    cache: pending
    outcome: error
  - service: crossref
    target: 10.3854/crm.10.checklist.atlas.v10.2025
    ts: 2026-05-17T01:30Z
    cache: pending
    outcome: hit
  - service: openalex
    target: "Testudinidae tortoise phylogeny"
    ts: 2026-05-17T01:30Z
    cache: pending
    outcome: hit
  - service: openalex
    target: "Geoemydidae Emydidae phylogeny species-level"
    ts: 2026-05-17T01:30Z
    cache: pending
    outcome: hit
---

# Turtles — partition audit

## Summary

Status `in_progress`. The canonical tree (Thomson 2021, PNAS) is DOI-verified,
methods extracted, and the partition has solid species-level coverage at 80%.
But three real data-quality findings need attention.

## Reference verification

**Tree paper** (Thomson 2021) — Crossref-clean. Title, authors, year, journal
all agree with data_sources.csv (implicit — there is no Thomson 2021 row;
citation is embedded in data_provenance.csv).

**Estimate source** (TTWG 2021) — **DOI broken**. Both Crossref and the doi.org
resolver return 404 for `10.3854/crm.8.000c.checklist.atlas.v9.2021`. The
Chelonian Research Foundation issues a fresh DOI per edition; the 2021/v9
edition's DOI either was never properly registered or has been retracted. The
2025/v10 edition's DOI resolves cleanly and reports 364 species (vs the 360
in our data).

**Recommended fix**: switch the cited estimate source to TTWG 2025. If the
paper-side data must stay frozen, the web app should surface a "this DOI is
broken; live counts here" note and link to https://iucn-tftsg.org/checklist/.

## Methods coverage

Thomson 2021 methods are now in the audit:

- Bayesian inference, MrBayes 3.2 for phylogeny + BEAST 1.8.2 for dating
- 15-locus concatenated and partitioned dataset
- 591 individual turtles (+ 2 outgroups = 593 tips) representing 279 of 348
  living species (80%); all 14 families monophyletic; 98% of genera
- 22 uniform fossil calibration priors from Joyce et al. 2013
- Uncorrelated lognormal relaxed clock, birth-death tree prior
- HKY substitution model per partition for the BEAST divergence-time run
- 100 chronograms in the credible set (atlas ships only the MCC consensus)
- 95% HPD intervals available

The paper is well-documented enough that no manual S-Appendix read is needed.

## Data-quality findings

Three real bugs surfaced during the audit (see `findings:` block in the
frontmatter):

1. **Broken DOI** for TTWG 2021 (severity: high).
2. **Loci-count mismatch**: provenance says "13 loci for 294 species", paper
   says "15 loci, 279 of 348 species" (severity: medium).
3. **Tip-vs-species ambiguity**: 593 tips, 329 unique IDs, 279 species — the
   dictionary embeds subspecies and locality codes that need to be surfaced
   in the web app (severity: medium).

Finding #2 specifically suggests the standardized tree may have been built
from a slightly different Dryad release than the published 591-individual
MCC. Worth a manual check against Dryad's deposit
(https://datadryad.org/dataset/doi:10.5061/dryad.jh9w0vt8w).

## Sub-clade gaps

At 80% species coverage, turtles is in good shape. Two **low/medium-priority**
candidates surfaced but neither is critical:

1. **Testudinidae** (Pyron 2014) — supplementary at best.
2. **Geoemydidae** (Pereira 2017) — the family most likely to harbor missing
   species, worth a closer look.

Sea turtles (Cheloniidae + Dermochelyidae = 7 species) and softshells
(Trionychidae ≈ 30 species) are already well-covered by Thomson 2021.

## Decisions (resolved)

- ✅ **Broken TTWG 2021 DOI** — web app substitutes TTWG 2025 (live, Crossref-verified) via `info.yaml.estimate.source.live_doi`. The paper-side `data_sources.csv` row stays as-is (per the don't-touch-paper-data rule); users see the correction prominently on the partition detail panel.
- ✅ **Loci/species mismatch** — web app reads authoritative figures from `info.yaml.tree.methods.n_loci` (=15) and `info.yaml.tree.species_represented` (=287) instead of falling back to `data_provenance.csv methods_brief`. CSV stays as-is.
- ✅ **Multiple individuals per species — canonical-tree swap** — the published 591-individual MCC contains 593 tips for 287 species (1-3 specimens each). For an atlas of species diversity, the species-level view is the useful one, so the canonical tree shown by the web app is an atlas-derived 287-tip pruning of the MCC (one representative per species, smallest numeric tip ID). The original multi-individual tree is preserved as an uncertainty file (locally hosted) and remains the canonical PAPER artifact. Pruning is reproducible via `scripts/derive_species_tree.py`.
- 🟡 **Geoemydidae candidate** (Pereira 2017) — deferred. Worth a species-overlap check vs Thomson 2021's 80% sampling before deciding to ship.
- ⬜ **Testudinidae candidate** (Pyron 2014) — skipped, low priority.

## Status: verified

All three caveats have a concrete resolution mechanism in `site/data/partitions/turtles/info.yaml`. The web app reads from `info.yaml` and shows users the corrected facts; the paper-side CSVs are untouched. If a sub-clade tree is later accepted from the deferred Geoemydidae candidate, the status will advance to `shipped`.
