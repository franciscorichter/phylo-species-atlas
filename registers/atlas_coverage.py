"""atlas_coverage -- the whole coverage layer of the phylogenetic atlas.

Every published coverage figure in the paper is ONE call of this function
against Data/coverage_basis.csv.  Numerator and denominator are column sums of
the SAME 62-row file, which is what makes a phantom denominator impossible.

    from atlas_coverage import load_basis, atlas_coverage
    b = load_basis("Data/coverage_basis.csv")
    euk = lambda r: r["is_prokaryotic"] == "no"
    atlas_coverage(b, euk)                                   # C0
    atlas_coverage(b, euk, units_col="c1_units")             # C1
    atlas_coverage(b, lambda r: euk(r) and r["is_chronogram"] == "yes",
                   units_col="c1_units", denominator_predicate=euk)   # C1d
    atlas_coverage(b, lambda r: euk(r) and r["uncertainty_available"] == "yes",
                   units_col="c1_units", denominator_predicate=euk)   # C1u
    atlas_coverage(b, lambda r: euk(r) and r["interval_recoverable"] == "yes",
                   units_col="c1_units", denominator_predicate=euk)   # C1b

c1_units is the number of distinct accepted names among the molecular tips
(c1_basis = enumerated), c0_units where every tip is molecular (all_tips), and
the expectation c0_units * f_mol only where neither holds (expectation).
"""
import csv


def load_basis(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _num(x, default=0.0):
    try:
        return float(str(x).replace(",", ""))
    except (TypeError, ValueError):
        return default


def atlas_coverage(basis, predicate=None, weight_col=None,
                   denominator_predicate=None, units_col="c0_units"):
    predicate = predicate or (lambda r: True)
    denominator_predicate = denominator_predicate or predicate
    num, n = 0.0, 0
    for r in basis:
        if not predicate(r):
            continue
        w = 1.0 if weight_col is None else _num(r.get(weight_col), 0.0)
        num += _num(r.get(units_col)) * w
        n += 1
    den = sum(_num(r.get("described")) for r in basis
              if denominator_predicate(r))
    return {"numerator": num, "denominator": den, "n_partitions": n,
            "percent": (100.0 * num / den) if den else None}
