# living/ — the monthly watch of the atlas

The atlas counts described eukaryotic species on the largest published species-level tree of each
group. New trees are published every month, so the counts age. `watch.py` finds the trees deposited
in a date window, measures each with the atlas's own rules, and screens it against the reference tree
of its group. It changes nothing in the atlas: every month it writes a run folder, and a curation step
decides, with evidence, whether a candidate replaces or adds a reference tree.

## What one run does

1. **Harvest.** Records of the window from Zenodo, Figshare and Dryad that match a phylogeny query,
   and the tree files of each (by file name, also inside zip, tar and gzip archives). Figshare's
   search is full-text and returns thousands of unrelated records a month, so its files (and
   Dryad's) are listed only for records whose title or abstract names a phylogeny, chronogram or
   divergence times; Zenodo returns its file lists with the search. Zenodo is queried in 7-day
   slices, and each slice's record count is checked against the total Zenodo reports. A shortfall,
   or a source that fails, puts **INCOMPLETE** at the top of the report with the numbers, and the
   run exits with status 2, so the timer shows a failed run; the report is still committed.
2. **Measure.** The first tree of each file: tips, species keys by the atlas rule
   (`scripts/canonical_trees.py`, decisions D4 and D13), the number of trees in the file, and the
   root-to-tip coefficient of variation.
3. **Assign.** Each genus is classified in the Catalogue of Life release of the paper (COL26.9 base,
   ChecklistBank 316321) through the ChecklistBank API, and each species to the group whose
   definition in `Data/col_denominators.csv` contains it. A genus name shared by homonyms is
   classified through the tree's own species of that genus. The tree's group is the group of most of
   its species (at least 60%; outgroups are expected).
4. **Screen** against the two reference trees of the group (decision D65, read from
   `Data/coverage_basis.csv`): the **primary** tree, the largest tree dated in absolute time, and the
   **secondary** tree, the largest tree of any kind. Each database is judged on its own; a larger
   undated tree never replaces a dated one.

   | Verdict | Database | Condition |
   |---|---|---|
   | LARGER | secondary | species in the group ≥ 0.9 × the secondary tree's accepted names |
   | FIRST | secondary | the group has no reference tree, and the tree holds ≥ 20 of its species |
   | DATED-LARGER | primary | ultrametric (CV ≤ 0.05) and ≥ 0.9 × the primary tree |
   | DATED-FIRST | primary | ultrametric, the group has no primary tree, and ≥ 20 species |
   | SET | primary | a set of ultrametric trees as large as the primary tree, which archives no date uncertainty |
   | REDEPOSIT? | either | added when ≥ 95% of its species are in the reference tree and < 5% are new: most likely the reference tree deposited again by a later study |

   Ultrametric stands in for dated: curation reads whether the time scale is absolute. The ledger
   also records whether a tree file carries node intervals (BEAST `height_95%_HPD`, MCMCTree
   `[&95%=…]`), which are archived date uncertainty.

   Species keys are not yet accepted names, hence the 0.9: curation counts accepted names.

Output in `runs/<since>_<until>/`: `report.md` (candidates, and Dryad records that could not be
measured), `ledger.csv.gz` (every tree file seen, with its measurements and verdict),
`records.csv.gz` (every record harvested), `run.log`. Downloads and the genus cache are in `cache/`, not tracked.

## Running it

```sh
python3 living/watch.py run --last-month              # the previous calendar month
python3 living/watch.py run --since 2026-06-01 --until 2026-10-14
python3 living/watch.py report living/runs/2026-09-01_2026-09-30   # rebuild a report from its files
```

In the ledger, trees that are not species-level (mostly gene trees) are summarised as one row per
downloaded file when a file holds more than 20 of them; the report lists one row per deposit and
group, with its largest candidate tree.

Standard library only; it needs the repository's `Data/` tables and `scripts/canonical_trees.py`,
not `Data/raw/`, so it runs on any host with network access.

## Search of the whole archive

The watch sees only the month's deposits. `watch.py search` searches the whole history of Zenodo,
Figshare and Dryad group by group, with each group's Catalogue of Life taxa (the `+` names of
`Data/col_denominators.csv` and the group's name) and the phylogeny terms, and processes every
record it has not seen before exactly as the watch does. `state/seen.csv.gz` is the ledger of every
record either mode has processed, so each record is read once; the first search covers the backlog
and later ones only what is new to the ledger. Output in `runs/search_<date>/`.

```sh
python3 living/watch.py search                       # every eukaryotic group
python3 living/watch.py search --groups Bivalves Fungi
```

## From a candidate to a change in the atlas

A candidate is a measurement, not a decision. `curate.py` measures the rest of what can be measured
for one candidate file and writes `curation/<group>/<deposit>_<file>.md` and `.json`: accepted
species (each species key matched in COL26.9 through ChecklistBank: accepted names count as
themselves, synonyms as their accepted species, genus-level matches stay unresolved), the tips
outside the group by name, the root-to-tip variation and root height, the number of trees in the
file, and the deposit's own description. A curator then reads the source study for what cannot be
measured: which tips were placed with molecular data, whether and how the tree is dated, what date
uncertainty the archive carries, and whether another file of the deposit is the one to use. The
verdict (replace the reference tree, first reference tree, or not adopted), its rule and its
evidence are added to the dossier.

```sh
python3 living/curate.py --group Scorpiones --source zenodo --record 20973872 \
    --file living/cache/downloads/<file> --member <path inside the archive>
```

An adopted tree enters the atlas as the revision's changes did: each decision recorded with its
rule, the accepted names counted again by `scripts/match_col.py`, `rebuild_atlas.py`,
`verify_atlas.py` and the release's `reproduce_paper.py` run, and the data release published only
after Francisco approves.

## Known limits

- **Dryad** serves files only to an API account. The watch takes an OAuth token (valid 10 h) from
  the account's client ID and secret: `DRYAD_CLIENT_ID`/`DRYAD_CLIENT_SECRET`, else
  `~/.config/atlas-watch/dryad.json` (chmod 600, on whitebox), else `vault get dryad`. Without
  them, Dryad records are listed in the report, not measured.
- **Trees without deposits** (figures only, or supplementary files at the publisher) are not seen.
- The query is textual; a deposit whose metadata never says phylogeny, chronogram or divergence
  times is missed.
- Only the first tree of a file is measured; a file holding several different trees counts as one.

## Schedule

On whitebox, a systemd user timer runs `monthly.sh` on the 2nd of each month at 04:30: it pulls
branch `main` into a sparse clone at `~/atlas-watch/phylo-atlas` (only `living/`, `scripts/`
and the `Data/` tables the watch reads), runs the previous calendar month and the search of the
whole archive, and commits and pushes the run folders and the seen-record ledger. Setup, once:

```sh
living/monthly.sh --setup                         # the sparse clone
mkdir -p ~/.config/systemd/user
cp ~/atlas-watch/phylo-atlas/living/systemd/atlas-watch.* ~/.config/systemd/user/
systemctl --user daemon-reload && systemctl --user enable --now atlas-watch.timer
```
