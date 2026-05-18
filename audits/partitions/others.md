---
partition: Others
category: Not yet represented
audit_version: 1
status: deferred
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: Atlas-aggregated
  doi: null
  doi_verified: na
  year_published: null
  cited_value: 470
  newer_authoritative_known: false

canonical_tree:
  filename: null
  status: not-shipped
  reason: |
    Synthetic catch-all partition aggregating ~392 described species across
    8 small lineages NOT contained in any other partition: Tuatara (1),
    Coelacanths (2), Lungfishes (6), Onychophora (~200), Placozoa (~4),
    Loricifera (~39), Cycliophora (3), Chaetognatha (~130). Each is small
    enough that no single tree exists across all of them; phylogenies are
    scattered across many papers. Bundled here so the atlas's partition
    is collectively exhaustive at the species level.

methods: null
shipped_subclades: []
candidate_subclades: []
resolutions: []

website:
  surface_methods: false
  surface_uncertainty: false
  surface_candidates: false
  badge: deferred
---

# Others — partition audit (deferred)

**Status: deferred — synthetic aggregate partition.**

This partition exists to make the atlas's chart **collectively exhaustive
at the species level**. Without it, ~392 species in 8 small lineages would
be in NO partition at all.

## Constituent clades

| Clade | Described | Estimated | Notes |
|---|---:|---:|---|
| Tuatara (Sphenodontia) | 1 | 1 | *Sphenodon punctatus* — only living rhynchocephalian, sister to Squamata |
| Coelacanths (Coelacanthiformes) | 2 | 2 | *Latimeria chalumnae* + *L. menadoensis* — "living fossils" |
| Lungfishes (Dipnoi) | 6 | 6 | Three genera in three families, sister to tetrapods |
| Onychophora (velvet worms) | 200 | 250 | Two families; sister to Tardigrada+Arthropoda |
| Placozoa | 4 | 20 | *Trichoplax* + 3 more cryptic genetic species |
| Loricifera | 39 | 60 | Meiofaunal phylum discovered 1983 |
| Cycliophora | 3 | 5 | Symbionts on lobster mouthparts |
| Chaetognatha (arrow worms) | 130 | 200 | Marine planktonic predators |

## Why these were missing

The atlas's main partitions (Mammals, Birds, Squamates, Turtles,
Crocodilians, Fish, Sharks, Hagfish, Lampreys, Amphibians for Vertebrates;
Cnidaria, Echinoderms, Bivalves, Gastropods, Sponges, Bryozoa, Cephalopods,
Annelids, Nematodes, Nemerteans, Tunicates, Ctenophores, Flatworms for
Other animals; etc.) collectively miss these 8 lineages because:

- **Tuatara** doesn't fit Squamates (sister to it, not within).
- **Coelacanths / Lungfishes** are sarcopterygian fishes; "Fish" in
  data_estimates uses Eschmeyer's Catalog of Fishes (ray-finned only).
- **Onychophora / Tardigrada** are panarthropods but not arthropods.
  Tardigrada has its own entry; Onychophora doesn't.
- **Placozoa / Loricifera / Cycliophora / Chaetognatha** are small
  phyla without category-level homes in the existing taxonomy.

## What "verified" would require

For the synthetic 'Others' partition to advance to verified, each
constituent would need its own DOI-verified estimate source. For now,
deferred — the species counts are stable monograph-level figures, just
not synthesized from a single source.
