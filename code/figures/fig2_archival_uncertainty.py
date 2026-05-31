"""Figure 2: Archival uncertainty preservation across the atlas.

For each canonical source tree (S5), we classify the dating-status
information that the published source actually preserves, and then plot
how those classes accumulate as a function of the tree's publication
year. The story this figure tells:

  * undated trees                — ``dated == no`` in S5 (phylogenomic or
    supertree topologies for which the source publication does not ship
    a time-calibrated chronogram). Time-dependent analyses cannot run on
    these without external re-dating.

  * point-dated trees            — the source publication ships node
    ages but no posterior distribution or bootstrap chronogram replicate
    set, so divergence-time uncertainty cannot be propagated downstream
    from archival material alone.

  * HPD-recoverable trees        — the source ships a posterior sample
    or a bootstrap chronogram distribution from which highest-posterior
    density intervals can be reconstructed. Only this class supports
    fully uncertainty-propagated downstream comparative analyses without
    re-running the dating pipeline.

We deliberately count trees, not species, because the question being
asked is methodological: how many of the building blocks the atlas
assembles preserve enough archival information for HPD-aware reuse.

Inputs:
  expert-panel-review-2026-05-29/table_s5_canonical_tree_provenance_DRAFT.csv

Outputs:
  figures/figure_archival_uncertainty.{png,pdf}
"""
import csv
import pathlib
import re
import matplotlib.pyplot as plt
import matplotlib as mpl

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
S5_PATH = REPO_ROOT / "audits" / "table_s5_canonical_tree_provenance.csv"
OUT_DIR = REPO_ROOT / "code" / "figures" / "outputs"

mpl.rcParams["font.family"] = "Helvetica"
mpl.rcParams["font.size"] = 9
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False

# ── Class palette (colour-blind safe) ──
COL_HPD     = "#1b7837"   # green   — HPD-recoverable
COL_POINT   = "#d9a441"   # amber   — point-dated only
COL_UNDATED = "#888888"   # grey    — undated topology

YEAR_RE = re.compile(r"(\d{4})")


def parse_year(row):
    """Best-effort publication year from tree_label, falling back to citation_key."""
    for field in ("tree_label", "citation_key", "notes"):
        text = (row.get(field) or "")
        m = YEAR_RE.search(text)
        if m:
            y = int(m.group(1))
            if 1990 <= y <= 2030:
                return y
    return None


def has_recoverable_distribution(row):
    """True iff the source archives material from which HPDs can be reconstructed.

    Strict rule: a posterior-sample count > 0 OR an explicit bootstrap /
    posterior distribution in the raw dating-method description. Mere
    reporting of a 95% HPD on a point chronogram is NOT enough — the
    distribution itself must be archived.
    """
    if (row.get("dated") or "").strip().lower() != "yes":
        return False

    posterior = (row.get("posterior_n") or "").strip()
    try:
        if posterior and int(posterior) > 0:
            return True
    except ValueError:
        pass

    raw = " ".join(
        (row.get(k) or "").lower()
        for k in ("dating_method_raw", "notes")
    )
    if re.search(r"\bbootstrap\b.*\b(distribution|chronogram|tree|sample)", raw):
        return True
    if re.search(r"\bposterior\b.*\b(distribution|sample|tree)", raw):
        return True
    if "1000 bootstrap" in raw or "1k.genus.trees" in raw:
        return True
    return False


def classify(row):
    if (row.get("dated") or "").strip().lower() != "yes":
        return "undated"
    if has_recoverable_distribution(row):
        return "hpd"
    return "point"


def main():
    rows = []
    with S5_PATH.open(newline="") as f:
        for r in csv.DictReader(f):
            year = parse_year(r)
            cls = classify(r)
            rows.append({
                "partition": r.get("partition_slug"),
                "label": r.get("tree_label"),
                "year": year,
                "class": cls,
            })

    # Exclude the Condamine bulk row (row 27): it is a single bulk row
    # representing 218 re-dated family chronograms, so counting it as
    # "one tree in 2019" is misleading. We mention it in the caption.
    rows = [r for r in rows if r["partition"] != "(condamine_bulk)"]

    classified = [r for r in rows if r["year"] is not None]
    years = sorted({r["year"] for r in classified})
    span = list(range(min(years), max(years) + 1))

    counts = {"undated": [0] * len(span), "point": [0] * len(span), "hpd": [0] * len(span)}
    for r in classified:
        idx = r["year"] - span[0]
        counts[r["class"]][idx] += 1

    totals = {k: sum(v) for k, v in counts.items()}

    fig, ax = plt.subplots(figsize=(7.2, 3.6))

    bottom = [0] * len(span)
    for key, col, label in [
        ("undated", COL_UNDATED, f"Undated topology (n={totals['undated']})"),
        ("point",   COL_POINT,   f"Point-dated only, no recoverable distribution (n={totals['point']})"),
        ("hpd",     COL_HPD,     f"HPD-recoverable from archived posterior/bootstrap (n={totals['hpd']})"),
    ]:
        ax.bar(
            span,
            counts[key],
            bottom=bottom,
            color=col,
            edgecolor="white",
            linewidth=0.6,
            label=label,
            width=0.78,
        )
        bottom = [b + c for b, c in zip(bottom, counts[key])]

    ax.set_xlabel("Publication year of canonical source tree")
    ax.set_ylabel("Canonical source trees deposited")
    ax.set_xticks(span)
    ax.set_xticklabels([str(y) for y in span], rotation=0)
    ax.set_yticks(range(0, max(bottom) + 2))
    ax.set_ylim(0, max(bottom) + 1)

    total = sum(totals.values())
    pct_hpd = 100 * totals["hpd"] / total
    ax.text(
        0.99, 0.97,
        f"{totals['hpd']} of {total} canonical trees ({pct_hpd:.0f} %) preserve a recoverable\n"
        f"divergence-time distribution. Condamine et al. (2019) re-dated family\n"
        f"chronograms (218 trees) are reported separately as point-dated.",
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=7, color="#333333",
    )

    ax.legend(loc="upper left", frameon=False, fontsize=8)
    ax.grid(axis="y", linewidth=0.4, alpha=0.4)
    ax.set_axisbelow(True)

    OUT_DIR.mkdir(exist_ok=True)
    out_png = OUT_DIR / "figure_archival_uncertainty.png"
    out_pdf = OUT_DIR / "figure_archival_uncertainty.pdf"
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")
    print()
    print(f"summary: {totals}  (total = {total})")
    for y in span:
        idx = y - span[0]
        print(f"  {y}: hpd={counts['hpd'][idx]:>2}  point={counts['point'][idx]:>2}  undated={counts['undated'][idx]:>2}")


if __name__ == "__main__":
    main()
