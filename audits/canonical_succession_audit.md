# Path 1 audit — per-partition canonical succession (2007–2026)

**Run:** workflow `wf_6200f856-0e9`, 122 agents, ~23 min wallclock, ~5.15 M subagent tokens.

## Headline

- Partitions surveyed: **26** (26 represented atlas partitions in S5).
- Total proposed canonicals: **96**.
- Total **accepted** (verified=true, confidence!=rejected): **95**.
- Total rejected (DOI did not resolve OR Crossref returned a different paper): **1**.

Confidence breakdown of accepted rows:

| Confidence | Count |
|---|---|
| high   | 41 |
| medium | 53 |
| low    | 1 |

## Per-partition coverage

| Partition | Proposed | Accepted | Rejected | H / M / L | Accepted citation_keys |
|---|---:|---:|---:|---|---|
| Amphibians | 4 | 4 | 0 | 1 / 3 / 0 | frost_2006, pyron_wiens_2011, pyron2014, jetz2018interplay |
| Archaea | 5 | 5 | 0 | 3 / 1 / 1 | parks2018gtdb, parks2020gtdb, rinke2021gtdb, parks2022gtdb, parks2025gtdb |
| Bacteria | 4 | 4 | 0 | 3 / 1 / 0 | hug2016, parks2018gtdb, parks2021gtdb, parks2025gtdb |
| Birds | 5 | 5 | 0 | 3 / 2 / 0 | hackett2008, jetz2012, prum_2015, stiller2024, mctavish2025 |
| Bivalves | 3 | 3 | 0 | 0 / 3 / 0 | whelan_geneva_graf_2011, lopes-lima-2017, pfeiffer2019unioverse |
| Bryophytes | 1 | 1 | 0 | 1 / 0 / 0 | bechteler2023 |
| Bryozoa | 3 | 3 | 0 | 1 / 2 / 0 | waeschenbach2012, orr2021bryozoa, orr2022bryozoa |
| Cephalopods | 2 | 2 | 0 | 1 / 1 / 0 | lindgren2012, basava2024 |
| Cnidaria (corals) | 5 | 5 | 0 | 3 / 2 / 0 | kitahara2010, huang2012, huang_roy_2015, quattrini2020, vaga2025 |
| Crocodilians | 3 | 2 | 1 | 0 / 2 / 0 | oaks2011, pan2021 |
| Crustaceans (Decapoda) | 3 | 3 | 0 | 1 / 2 / 0 | bracken2009, bracken_grissom_2014, wolfe2019 |
| Diatoms | 3 | 3 | 0 | 3 / 0 / 0 | theriot2010, nakov2018, alverson2025 |
| Echinoderms (echinoids) | 3 | 3 | 0 | 1 / 2 / 0 | kroh_smith_2010, mongiardino_koch_thompson_2021, mongiardino-koch-2022 |
| Ferns | 4 | 4 | 0 | 2 / 2 / 0 | schuettpelz2007, Lehtonen2011, testo_sundue_2016, nitta2022ftol |
| Ray-finned fish | 5 | 5 | 0 | 3 / 2 / 0 | santini2009, near2012, betancur2013fish, betancur2017, rabosky2018 |
| Fungi | 4 | 4 | 0 | 0 / 4 / 0 | james2006, hibbett2007, varga2019, li2021 |
| Gastropods | 2 | 2 | 0 | 2 / 0 / 0 | zapata2014, cunha_giribet_2019 |
| Insects | 4 | 4 | 0 | 1 / 3 / 0 | misof2014, chesters2017, Chesters2020, chesters_2023_insectphylo |
| Mammals | 7 | 7 | 0 | 1 / 6 / 0 | bininda_emonds_2007, fritz2009, kuhn2011, faurby_svenning_2015, faurby2018phylacine, upham2019, alvarez_carretero_2022 |
| Nematodes | 4 | 4 | 0 | 2 / 2 / 0 | holterman2006, vanMegen2009, smythe2019, qing2025nematoda |
| Seed plants | 8 | 8 | 0 | 4 / 4 / 0 | smith_beaulieu_donoghue_2009, zanne2014, magallon2015, qian_jin_2016_phytophylo, smith_brown_2018, janssens2020, baker2022paftol, zuntini2024angiosperms |
| Sharks | 2 | 2 | 0 | 0 / 2 / 0 | velez-zuazo_agnarsson_2011, stein2018 |
| Spiders | 3 | 3 | 0 | 1 / 2 / 0 | garrison2016, wheeler2017, kallal2021 |
| Sponges (Demospongiae) | 3 | 3 | 0 | 1 / 2 / 0 | redmond2013, plese2021, lavrov2023 |
| Squamates | 3 | 3 | 0 | 2 / 1 / 0 | pyron2013squamata, zheng_wiens_2016, Tonini 2016 |
| Turtles | 3 | 3 | 0 | 1 / 2 / 0 | guillon2012turtles, pereira2017, thomson2021 |

## Partitions needing human review (low confidence or zero accepted)

### Archaea (archaea)
- proposed=5 accepted=5 rejected=0 H/M/L=3/1/1

> Search strategy: (1) Verified both supplied anchors (Parks 2018 GTDB v1, DOI 10.1038/nbt.4229; Parks 2025 GTDB R10, DOI 10.1093/nar/gkaf1040) via Crossref. (2) Searched Google Scholar / Crossref for pre-2018 species-level archaeal reference trees — the only candidates (Wu 2009 GEBA, Hug 2016 'New view of the tree of life') do not qualify as Archaea-specific species-level canonicals: Wu 2009 sampled only 56 prokaryotes (mostly bacteria, ~13 archaea — below the 50-tip threshold), and Hug 2016 is a universal tree of 3,083 organisms with archaea as a minor (~few hundred tips) subset using only 16 ...


## Rejected entries

| Partition | Rejected citation_key |
|---|---|
| Crocodilians | condamine2019crocs |

---

Source: `path1_canonical_succession.workflow.js` (in session tmp). Discovery + verification + structured-output pipeline; all DOIs Crossref-checked before acceptance.
