#!/usr/bin/env python3
"""Rebuild figN_atlas_growth.csv with HONEST historical counting.

The original CSV used a "today's canonical only" rule: a partition appears in
year X only if its CURRENT canonical tree was published by then. Under that
rule, Birds appear in 2025 (McTavish) even though Jetz 2012 had ~10,000
species available throughout 2012-2024.

The honest rule: a partition appears in year X if ANY species-level tree
meeting atlas inclusion criteria (peer-reviewed, downloadable, species-level)
was published by then. For each partition, the species count at year X is
the best tree available at year X — typically the most-tip-rich, dated tree.

The historical record below is curated from the audit MDs (where they have
`alternative` blocks), data_provenance.csv, and verified literature
(DOI-checked via WebSearch in this session).
"""
from __future__ import annotations
import csv
import pathlib

DATA = pathlib.Path(__file__).resolve().parent.parent  # Data/
OUT_CSV = DATA.parent / "Paper" / "inputs" / "figN_atlas_growth.csv"

# Per-partition history: list of (year, n_tips, study, dated)
# "Best available species-level tree for the partition" at the year of publication.
# When a later study supersedes an earlier one, the later year's row uses the
# later tip count.
# Sources: verified via Crossref/PubMed WebSearch where uncertain;
# numbers reflect the published n_tips at the time of publication.
HISTORY = {
    "Mammals": [
        (2007, 4510, "Bininda-Emonds et al. 2007", True),
        (2019, 5911, "Upham et al. 2019", True),
    ],
    "Birds": [
        (2012, 9993, "Jetz et al. 2012", True),
        (2025, 10824, "McTavish et al. 2025", True),
    ],
    "Amphibians": [
        (2014, 3309, "Pyron 2014", True),
        (2018, 7239, "Jetz & Pyron 2018", True),
    ],
    "Fish": [
        (2013, 1410, "Betancur-R et al. 2013", True),
        (2018, 31516, "Rabosky et al. 2018", True),
    ],
    "Squamates": [
        (2013, 4161, "Pyron Burbrink & Wiens 2013", True),
        (2016, 9755, "Tonini et al. 2016", True),
    ],
    "Sharks": [
        (2012, 596, "Naylor et al. 2012", False),
        (2018, 1192, "Stein et al. 2018", True),
    ],
    "Turtles": [
        (2017, 286, "Pereira et al. 2017", True),
        (2021, 327, "Thomson et al. 2021", True),
    ],
    "Crocodilians": [
        (2011, 23, "Oaks 2011", True),
    ],
    "Hagfish": [
        (2013, 33, "Fernholm et al. 2013", True),
    ],
    "Lampreys": [
        (2025, 46, "Hughes et al. 2025", True),
    ],
    # Plants
    "Seed plants": [
        (2014, 32223, "Zanne et al. 2014", True),
        (2018, 356305, "Smith & Brown 2018", True),
    ],
    "Ferns": [
        (2016, 4000, "Testo & Sundue 2016", True),
        (2022, 6234, "Nitta et al. 2022", True),
    ],
    "Bryophytes": [
        (2023, 533, "Bechteler et al. 2023", True),
    ],
    "Conifers": [
        (2018, 578, "Leslie et al. 2018", True),
    ],
    "Lycophytes": [
        (2024, 300, "Zhou et al. 2024", True),
    ],
    "Orchids": [
        (2015, 170, "Givnish et al. 2015", True),
        (2023, 1475, "Thompson et al. 2023", True),
    ],
    "Solanaceae": [
        (2013, 1076, "Sarkinen et al. 2013", True),
    ],
    "Cacti": [
        (2026, 1063, "Thompson et al. 2026", True),
    ],
    "Grasses": [
        (2025, 1153, "GPWG III 2025", False),
    ],
    "Asteraceae": [
        (2019, 256, "Mandel et al. 2019", False),
    ],
    "Eucalyptus": [
        (2024, 399, "Crisp et al. 2024", False),
    ],
    # Arthropods
    "Insects": [
        (2017, 49338, "Chesters 2017", False),
        (2023, 53596, "Chesters et al. 2023 (insectphylo.org)", False),
    ],
    "Butterflies": [
        (2019, 2000, "Kawahara et al. 2019", True),
        (2023, 2258, "Kawahara et al. 2023", True),
    ],
    "Bees": [
        (2024, 4586, "Henriquez-Piskulich et al. 2024", True),
    ],
    "Spiders": [
        (2016, 1456, "Garrison et al. 2016", True),
    ],
    "Hemiptera": [
        (2018, 193, "Johnson et al. 2018", False),
    ],
    "Termites": [
        (2024, 138, "Hellemans et al. 2024", False),
    ],
    "Crustaceans": [
        (2019, 95, "Wolfe et al. 2019", False),
    ],
    "Ants": [
        (2022, 83, "Romiguier et al. 2022", False),
    ],
    "Myriapods": [
        (2023, 104, "Benavides et al. 2023", True),
    ],
    # Other animals
    "Corals": [
        (2025, 289, "Vaga et al. 2025", True),
    ],
    "Cephalopods": [
        (2024, 81, "Basava et al. 2024", True),
    ],
    "Gastropods": [
        (2014, 54, "Zapata et al. 2014", False),
    ],
    "Sponges": [
        (2023, 120, "Lavrov et al. 2023", True),
    ],
    "Cnidaria": [
        (2024, 1814, "DeBiasse et al. 2024", False),
    ],
    "Echinoderms": [
        (2022, 66, "Mongiardino Koch & Thompson 2022", True),
    ],
    "Bryozoa": [
        (2022, 720, "Orr et al. 2022", False),
    ],
    "Tunicates": [
        (2018, 63, "Delsuc et al. 2018", True),
    ],
    "Nemerteans": [
        (2014, 80, "Andrade et al. 2014", False),
    ],
    "Ctenophores": [
        (2017, 27, "Whelan et al. 2017", True),
    ],
    "Nematodes": [
        (2025, 166, "Qing et al. 2025", False),
    ],
    "Flatworms": [
        (2023, 83, "Brabec et al. 2023", True),
    ],
    "Bivalves": [
        (2019, 73, "Pfeiffer et al. 2019", False),
    ],
    "Annelids": [
        (2014, 60, "Weigert et al. 2014", False),
    ],
    # Microbes & protists
    "Bacteria": [
        # GTDB releases: R0 (2018) ~92K, R10 (2025) ~715K
        (2018, 92000, "Parks et al. 2018 (GTDB R0)", False),
        (2025, 136646, "Parks et al. 2025 (GTDB R10)", False),
    ],
    "Archaea": [
        (2018, 4500, "Parks et al. 2018 (GTDB R0)", False),
        (2025, 6968, "Parks et al. 2025 (GTDB R10)", False),
    ],
    "Fungi": [
        (2021, 1672, "Li et al. 2021", False),
    ],
    "Diatoms": [
        (2018, 244, "Nakov et al. 2018", True),
    ],
    "Ciliates": [
        (2024, 105, "Jiang et al. 2024", True),
    ],
    "Foraminifera": [
        (2011, 339, "Aze et al. 2011", True),
    ],
    "Dinoflagellates": [
        (2024, 80, "Cooney et al. 2024", False),
    ],
    "Brown algae": [
        (2010, 72, "Silberfeld et al. 2010", True),
    ],
    "Green algae": [
        (2021, 70, "Li et al. 2021", False),
    ],
    "Amoebozoa": [
        (2022, 113, "Tekle et al. 2022", False),
    ],
    "Choanoflagellates": [
        (2017, 47, "Carr et al. 2017", False),
    ],
    "Haptophytes": [
        (2010, 40, "Medlin et al. 2010", True),
    ],
}


def species_count_at(year: int) -> tuple[int, int, int]:
    """For a given year, return (n_partitions, n_dated_partitions, sum_best_species)."""
    n_part = 0
    n_dated = 0
    total_sp = 0
    for partition, hist in HISTORY.items():
        # Pick the best tree at or before this year (the latest one that is <= year).
        best = None
        for (y, tips, study, dated) in hist:
            if y <= year:
                if best is None or y > best[0]:
                    best = (y, tips, study, dated)
        if best is not None:
            n_part += 1
            if best[3]:
                n_dated += 1
            total_sp += best[1]
    return n_part, n_dated, total_sp


# Also compute cumulative trees (a different unit — historical accumulation
# of distinct studies, not "best per partition")
def cumulative_trees_at(year: int) -> tuple[int, int]:
    n_trees = 0
    n_dated_trees = 0
    for partition, hist in HISTORY.items():
        for (y, tips, study, dated) in hist:
            if y <= year:
                n_trees += 1
                if dated:
                    n_dated_trees += 1
    return n_trees, n_dated_trees


with OUT_CSV.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["year", "trees_cumulative", "dated_trees_cumulative",
                "partitions_covered", "species_in_canonical_trees"])
    for year in range(2010, 2027):
        n_part, _, sp = species_count_at(year)
        n_trees, n_dated = cumulative_trees_at(year)
        w.writerow([year, n_trees, n_dated, n_part, sp])
        print(f"  {year}: partitions_covered={n_part:2d}  best-tree-species={sp:>8,}  cum_trees={n_trees}  dated={n_dated}")
print(f"\nWrote {OUT_CSV}")
