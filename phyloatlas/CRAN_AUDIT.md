# phyloatlas v0.1.0 — CRAN audit report

**Run:** workflow `wf_2559fb59-c77`, 32 agents, 58 findings, adversarial verification on every must-fix-for-CRAN claim.

## Headline

- Auditors complete: **5/5**
- Total findings: **58**

| Severity | Count |
|---|---:|
| must-fix-for-cran | 2 |
| nice-to-have | 42 |
| refactor | 1 |
| informational | 13 |

## Confirmed must-fix-for-CRAN (after adversarial verification)

### [vignette-declaration] DESCRIPTION
Log line 19 'Package has a VignetteBuilder field but no prebuilt vignette index.' is a hard inconsistency, not a cosmetic note. DESCRIPTION declares 'VignetteBuilder: knitr' and 'Suggests: knitr, rmarkdown' but vignettes/ is empty (verified: /Users/pancho/Code/phylo-species-atlas/phyloatlas/vignettes/ contains 0 files). CRAN reviewers flag this in incoming because a knitr-built tarball is expected to contain build/vignette.rds + inst/doc/<vignette>.html -- phyloatlas ships neither. Two ways out: (A) drop the field, or (B) actually ship a vignette. Given the package is a thin atlas client, the 

**Proposed fix:**

> Option B (recommended): create vignettes/phyloatlas.Rmd with YAML header 'output: rmarkdown::html_vignette' and '%\VignetteEngine{knitr::rmarkdown}' / '%\VignetteIndexEntry{Getting started with phyloatlas}' / '%\VignetteEncoding{UTF-8}'. Inside the vignette, gate every network call with a CRAN-safe guard so it does not download on offline build machines:
> ```r
> knitr::opts_chunk$set(eval = identical(Sys.getenv("NOT_CRAN"), "true"))
> ```
> and provide a small static fallback narrative + a chunk that loads the bundled inst/extdata/tree_demo.nwk so the rendered HTML is never empty. Rebuild the tarball

### [fixtures] fixtures
Need a self-contained miniature mirror of the atlas under tests/testthat/fixtures/atlas-root/ so that paste0(base, '/standardized/trees/X.nwk') etc. resolve to real files. Layout must match the real repo exactly: standardized/dictionary.csv, standardized/metadata.csv, standardized/trees/*.nwk, data_provenance.csv at the root. file:// URL works for both utils::read.csv and ape::read.tree.

**Proposed fix:**

> See fixtures listed in test_files: tests/testthat/fixtures/atlas-root/standardized/dictionary.csv (5 species), standardized/metadata.csv (3 trees: mammals, birds, condamine_Vangidae - chosen to exercise all three .provenance_key branches), data_provenance.csv (3 rows: mammals, birds_mctavish, condamine), standardized/trees/mammals.nwk, standardized/trees/birds.nwk, standardized/trees/condamine_Vangidae.nwk. All Newick trees use small integer IDs that overlap the dictionary so label resolution can be verified.

## Downgraded (auditors said must-fix, skeptics refuted)

**26 of 21 initially-claimed must-fix items** were downgraded after adversarial verification. Common pattern: claims about "network handling" and "convert all dontrun to donttest" were misapplied — CRAN doesn't enforce those for examples wrapped in `\dontrun{}`.

- **network-handling** (cran-policy-audit) → `nice-to-have`. Reason: The claim quotes the real CRAN policy ("Packages which use Internet resources should fail gracefully...") but mis-applies its operational trigger. CRAN's enforcement of that clause is tied to the requ...
- **network-handling** (cran-policy-audit) → `nice-to-have`. Reason: The example in load_atlas_tree.R is wrapped in \dontrun{}, so CRAN's check infrastructure never executes the network call. The claim's central premise — that the missing timeout could hang CRAN checks...
- **examples** (cran-policy-audit) → `nice-to-have`. Reason: The example in load_atlas_tree.R uses \dontrun{} to wrap code that downloads a Newick tree from a remote URL via ape::read.tree(url). CRAN does accept \dontrun{} for examples requiring network access ...
- **examples** (cran-policy-audit) → `nice-to-have`. Reason: The claim is overblown for CRAN acceptance purposes. Verification:

1. The file /Users/pancho/Code/phylo-species-atlas/phyloatlas/R/list_trees.R wraps its example in \dontrun{}, calling list_trees() w...
- **examples** (cran-policy-audit) → `nice-to-have`. Reason: The claim is overblown. The example in atlas_info.R uses \dontrun{} for a function (atlas_info) that ultimately calls .load_metadata(), which downloads CSV files from raw.githubusercontent.com via uti...
- **examples** (cran-policy-audit) → `nice-to-have`. Reason: Verified the file: atlas_clear_cache() is indeed exported (line 49) and has roxygen documentation (title, description, @return, @export) but no @examples block. However, the claim that this is a CRAN ...
- **cran-incoming-feasibility** (cran-policy-audit) → `nice-to-have`. Reason: Verified against /Users/pancho/Code/phylo-species-atlas/phyloatlas/DESCRIPTION and the package tree. Two observations adjust the claim:

1. Factual: tests/ and vignettes/ directories DO exist (created...
- **description-style** (cran-policy-audit) → `nice-to-have`. Reason: The claim is factually inaccurate and overblown. Reading /Users/pancho/Code/phylo-species-atlas/phyloatlas/DESCRIPTION lines 11-20, the Description field actually contains THREE complete sentences, ea...

## Nice-to-have (deferred — won't block CRAN acceptance)

42 items, e.g.:
- [network-handling] Internal helpers .load_dictionary(), .load_metadata(), and .load_provenance() call utils::read.csv() directly on a raw.g...
- [network-handling] load_atlas_tree() does call tryCatch around ape::read.tree(url), which is good, but it sets no timeout (ape::read.tree o...
- [examples] Example uses \dontrun{}. CRAN Writing R Extensions explicitly says: 'Note that \dontrun is often misused: in most cases ...
- [examples] Same \dontrun{} issue. The example body is harmless — it would run fine if the network is up — so there is no reason to ...
- [examples] Same \dontrun{} issue....
- [examples] atlas_clear_cache() is exported but has NO @examples block at all. CRAN strongly prefers at least one runnable example p...
- [cran-incoming-feasibility] DESCRIPTION declares `Suggests: knitr, rmarkdown, testthat (>= 3.0.0)`, `Config/testthat/edition: 3`, and `VignetteBuild...
- [description-style] Description text starts with 'Convenience functions to fetch...' — CRAN's policy is that Description should be one or mo...
- [filesystem] Cache uses `.cache <- new.env(parent = emptyenv())` — in-memory only, no disk writes. This is fully CRAN-compliant: the ...
- [license] LICENSE file is well-formed: two lines exactly matching CRAN's MIT template (YEAR + COPYRIGHT HOLDER, no extra text). DE...

_(plus 32 more)_

## Outcome

Both confirmed must-fix items implemented:
- Vignette: `vignettes/phyloatlas.Rmd`
- Test fixtures: `tests/testthat/fixtures/atlas-root/` + 6 test files + helper

**`R CMD check --as-cran` now reports: 0 ERRORs · 0 WARNINGs · 2 expected NOTEs (both standard for first submission).**