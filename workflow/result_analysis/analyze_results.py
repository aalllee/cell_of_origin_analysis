#!/usr/bin/env python3
"""
Standard result analysis for every sample: ranking figures, S curves, and cross-sample comparisons.

One report per (sample, gene set). Gene sets:

  baseline   every gene shared with the reference. Identical to the pipeline's correlate.py
             (verified to 1e-16), so results/<sample>/correlation.csv is reproduced exactly.
  var3k      the 3000 genes with the highest variance across the 447 (cell type, tissue) groups.

Everything is recomputed from results/<sample>/proj/fft_summaries/fft_<sample>_WPS.tsv.gz, so both
gene sets come from one source and the two modes are directly comparable.

Output layout
  <out>/<sample>/<mode>/   correlations.csv, top_bottom.png, spread.png, by_compartment.png,
                           by_tissue.png, by_celltype_origin.png, s_curves.png, s_curves.csv.gz
  <out>/comparison/<mode>/ scatter.png, rank_shift.png, tissue_order.png, spread_compare.png
  <out>/summary.csv        one row per (sample, mode)

Usage
  python analyze_results.py ../results/healthy ../results/cancer
  python analyze_results.py --all -o analysis
  python analyze_results.py --all --compare healthy cancer
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent

# ---- categories ------------------------------------------------------------------------------
ERY, PLT = "erythrocyte (anucleate)", "platelet (anucleate)"
BARRIER, LOWTURN = "behind blood barrier", "low turnover"
CONTROL_CATS = [ERY, PLT, BARRIER, LOWTURN]
CONTROL_TYPES = {
    ERY: {"erythrocyte"},
    PLT: {"platelet"},
    BARRIER: {"eye photoreceptor cell", "retinal bipolar neuron", "retina horizontal cell",
              "retinal ganglion cell", "retinal pigment epithelial cell"},
    LOWTURN: {"cardiac muscle cell", "slow muscle cell", "fast muscle cell", "cell of skeletal muscle"},
}
TYPE_TO_CAT = {t: k for k, types in CONTROL_TYPES.items() for t in types}
CAT_COLORS = {"blood/immune": "tab:red", "other": "tab:gray",
              ERY: "#2a78d6", PLT: "#e87ba4", BARRIER: "#1baf7a", LOWTURN: "#eda100"}
CAT_MARKERS = {ERY: "D", PLT: "s", BARRIER: "^", LOWTURN: "v"}
BAR_COLORS = {"immune": "tab:red", "other": "tab:gray", "negative control": "black"}
COMP_COLORS = {"immune": "tab:red", "epithelial": "tab:green",
               "endothelial": "tab:purple", "stromal": "tab:brown"}
SAMPLE_COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]

# where a cell type lives: haematopoietic / lymphoid tissues vs solid organs
LYMPHOID = {"Blood", "Bone_Marrow", "Spleen", "Thymus", "Lymph_Node"}
TISSUE_ORIGIN_COLORS = {"lymphoid": "tab:red", "solid": "tab:gray"}   # Blood/Bone_Marrow/Spleen/Thymus/Lymph_Node vs the rest
INK = "#333333"

plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 150, "savefig.facecolor": "white",
                     "font.size": 10, "axes.titlesize": 11, "legend.frameon": False,
                     "axes.spines.top": False, "axes.spines.right": False})


# ---- loading ---------------------------------------------------------------------------------
def fft_matrix(sample_dir):
    sample = sample_dir.name
    d = sample_dir / "proj" / "fft_summaries"
    for name in (f"fft_{sample}_WPS.tsv.gz", f"fft_{sample}_WPS.tsv"):
        if (d / name).exists():
            m = pd.read_csv(d / name, sep="\t", index_col=0)
            m.columns = m.columns.astype(float)
            return m[sorted(m.columns)]
    raise SystemExit(f"no FFT matrix for '{sample}' in {d}")


def split_group(g):
    cell_type, tissue = g.split("_", 1)
    return cell_type, tissue


def build_reference(ref_csv, compartment_csv):
    ref = pd.read_csv(ref_csv, index_col=0)
    ref.pop("gene_name")
    comp = pd.read_csv(compartment_csv, index_col=0)["compartment"].reindex(ref.columns)
    if comp.isna().any():
        raise SystemExit(f"groups without a compartment: {list(comp.index[comp.isna()])[:5]}")

    def category(g):
        ct = split_group(g)[0]
        return TYPE_TO_CAT.get(ct) or ("blood/immune" if comp[g] == "immune" else "other")

    cat = pd.Series({g: category(g) for g in ref.columns})
    lab = cat.map(lambda c: "negative control" if c in CONTROL_CATS
                  else "immune" if c == "blood/immune" else "other")
    tissue = pd.Series({g: split_group(g)[1] for g in ref.columns})
    celltype = pd.Series({g: split_group(g)[0] for g in ref.columns})

    # per cell type: which tissues is it seen in?
    origin = {}
    for ct, groups in celltype.groupby(celltype):
        tis = set(tissue[groups.index])
        n_lymph = len(tis & LYMPHOID)
        origin[ct] = ("blood & lymphoid only" if n_lymph == len(tis)
                      else "solid tissue only" if n_lymph == 0 else "mixed")
    return ref, comp, cat, lab, tissue, celltype, pd.Series(origin)


def pick_genes(mode, ref, index, var_n):
    "genes a mode correlates over, chosen from `index`"
    if mode == "baseline":
        return index
    return ref.loc[index].var(axis=1).nlargest(var_n).index


def gene_sets(ref, shared, var_n):
    modes = ["baseline", f"var{var_n // 1000}k"]
    return {m: pick_genes(m, ref, shared, var_n) for m in modes}


def hbox(ax, vals, pos, width=0.55):
    "horizontal boxplot; matplotlib >= 3.11 wants orientation=, older versions want vert="
    kw = dict(positions=[pos], widths=width, patch_artist=True, showfliers=False)
    try:
        return ax.boxplot(vals, orientation="horizontal", **kw)
    except TypeError:
        return ax.boxplot(vals, vert=False, **kw)


# ---- per-sample figures ----------------------------------------------------------------------
def fig_top_bottom(r, cat, comp, out, title, n=10):
    sel = pd.concat([r.iloc[:n], r.iloc[-n:]])
    shown = list(range(1, n + 1)) + list(range(len(r) - n + 1, len(r) + 1))
    y = np.arange(len(sel))[::-1].astype(float)
    y[n:] -= 0.8
    colors = [BAR_COLORS["negative control"] if cat[g] in CONTROL_CATS
              else BAR_COLORS["immune" if comp[g] == "immune" else "other"] for g in sel.index]
    fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(sel) + 1.8))
    ax.barh(y, sel.values, color=colors, alpha=0.9, height=0.75)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{k}. {ct} ({t.replace('_', ' ')})"
                        for k, (ct, t) in zip(shown, map(split_group, sel.index))], fontsize=8)
    longest = sel.abs().max()
    for yi, v in zip(y, sel.values):
        inside = abs(v) >= 0.35 * longest
        ax.annotate(f"{v:+.3f}", (v if inside else max(v, 0.0), yi),
                    xytext=((6 if v < 0 else -6) if inside else 7, 0), textcoords="offset points",
                    va="center", ha="left" if (v < 0 or not inside) else "right", fontsize=8,
                    color="white" if inside else INK, fontweight="bold" if inside else "normal",
                    clip_on=False)
    ax.axvline(0, color="k", lw=0.7)
    ax.axhline(y[n - 1] - 0.9, color="0.8", lw=1, ls="--")
    pad = 0.12 * (sel.max() - sel.min())
    ax.set_xlim(sel.min() - pad, sel.max() + pad)
    ax.set_xlabel("correlation with expression (r)")
    ax.set_title(f"{title}: strongest {n} and weakest {n} of {len(r)} groups")
    ax.legend(handles=[Patch(color=c, label=k) for k, c in BAR_COLORS.items()],
              fontsize=8, loc="lower left")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_spread(r, lab, out, title):
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for i, k in enumerate(BAR_COLORS):
        vals = r[lab[r.index] == k].values
        if not len(vals):
            continue
        bp = ax.boxplot(vals, positions=[i], widths=0.55, patch_artist=True, showfliers=False)
        bp["boxes"][0].set(facecolor=BAR_COLORS[k], alpha=0.35)
        bp["medians"][0].set(color=BAR_COLORS[k], lw=2)
        ax.scatter(rng.normal(i, 0.07, len(vals)), vals, s=7, color=BAR_COLORS[k], alpha=0.45, zorder=3)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(range(len(BAR_COLORS)))
    ax.set_xticklabels([f"{k}\n(n={int((lab == k).sum())})" for k in BAR_COLORS], fontsize=9)
    ax.set_ylabel("correlation with expression (r)")
    ax.set_title(f"{title}: correlations by label")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def _rank_boxes(ax, rank, keys, members, cat):
    rng = np.random.default_rng(0)
    for i, k in enumerate(keys):
        groups = members(k)
        vals = rank[groups]
        bp = ax.boxplot(vals, positions=[i], widths=0.55, patch_artist=True, showfliers=False)
        bp["boxes"][0].set(facecolor="tab:blue", alpha=0.35)
        bp["medians"][0].set(color="tab:blue", lw=2)
        jit = rng.normal(i, 0.07, len(vals))
        cats = np.array([cat[g] for g in groups])
        is_ctrl = np.isin(cats, list(CAT_MARKERS))
        xy = (jit, vals.values)
        ax.scatter(*(np.asarray(a)[~is_ctrl] for a in xy), s=7, color="k", alpha=0.3, zorder=3)
        for ck, mk in CAT_MARKERS.items():
            sel = cats == ck
            if sel.any():
                ax.scatter(*(np.asarray(a)[sel] for a in xy), s=60, marker=mk, color=CAT_COLORS[ck],
                           edgecolors="white", linewidths=0.6, zorder=6)


def fig_by_compartment(rank, comp, cat, out, title):
    keys = comp.groupby(comp).apply(lambda t: rank[t.index].median()).sort_values().index.tolist()
    fig, ax = plt.subplots(figsize=(8, 5))
    _rank_boxes(ax, rank, keys, lambda k: list(comp.index[comp == k]), cat)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([f"{k}\n(n={int((comp == k).sum())})" for k in keys])
    for lbl, k in zip(ax.get_xticklabels(), keys):
        lbl.set_color(COMP_COLORS.get(k, "k"))
    ax.invert_yaxis()
    ax.set_ylabel(f"rank (1 = strongest, {len(rank)} = weakest)")
    ax.set_title(f"{title}: rank by compartment")
    ax.legend(handles=[Line2D([], [], marker=m, ls="", color=CAT_COLORS[k], label=k,
                              markeredgecolor="white") for k, m in CAT_MARKERS.items()],
              fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_by_tissue(rank, tissue, cat, out, title):
    keys = tissue.groupby(tissue).apply(lambda t: rank[t.index].median()).sort_values().index
    fig, ax = plt.subplots(figsize=(16, 5.5))
    _rank_boxes(ax, rank, keys, lambda k: list(tissue.index[tissue == k]), cat)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([f"{k.replace('_', ' ')} ({int((tissue == k).sum())})" for k in keys],
                       rotation=45, ha="right")
    ax.invert_yaxis()
    ax.set_ylabel(f"rank (1 = strongest, {len(rank)} = weakest)")
    ax.set_title(f"{title}: rank by tissue, all cell types in the tissue; (n) = cell types")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_by_celltype_origin(rank, r, celltype, tissue, out, title, top_n=40):
    """one box per cell type (its tissues), each point coloured by the tissue it came from"""
    med = rank.groupby(celltype[rank.index]).median().sort_values()
    keys = list(med.index[:top_n])
    counts = celltype[rank.index].value_counts()

    fig, ax = plt.subplots(figsize=(9, 0.30 * len(keys) + 2))
    rng = np.random.default_rng(0)
    for i, ct in enumerate(keys):
        groups = list(celltype.index[(celltype == ct) & celltype.index.isin(rank.index)])
        vals = rank[groups]
        bp = hbox(ax, vals, i)
        bp["boxes"][0].set(facecolor="0.85", edgecolor="0.5")
        bp["medians"][0].set(color="0.3", lw=1.8)
        is_lymph = np.array([tissue[g] in LYMPHOID for g in groups])
        jit = rng.normal(i, 0.07, len(groups))
        for sel, key in ((is_lymph, "lymphoid"), (~is_lymph, "solid")):
            if sel.any():
                ax.scatter(np.asarray(vals)[sel], jit[sel], s=26, color=TISSUE_ORIGIN_COLORS[key],
                           alpha=0.9, edgecolors="white", linewidths=0.5, zorder=3)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([f"{ct}  ({counts[ct]})" for ct in keys], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f"rank (1 = strongest, {len(rank)} = weakest)")
    ax.set_title(f"{title}: strongest {len(keys)} of {len(med)} cell types\n"
                 f"one point per tissue the cell type was sampled in; (n) = tissues", fontsize=11)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=c, label=k, markeredgecolor="white")
                       for k, c in TISSUE_ORIGIN_COLORS.items()],
              fontsize=9, loc="lower right", title="tissue the cell comes from")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)

    agg = r.groupby(celltype[r.index]).agg(["mean", "count"]).sort_values("mean")
    lymph_frac = pd.Series({ct: np.mean([tissue[g] in LYMPHOID
                                         for g in celltype.index[(celltype == ct) & celltype.index.isin(r.index)]])
                            for ct in agg.index})
    agg["lymphoid_fraction"] = lymph_frac
    return agg


def fig_s_curves(fm, ref, genes, lab, band, out_png, out_csv, title):
    g = genes.intersection(fm.index)
    curves = pd.DataFrame({p: ref.loc[g].corrwith(fm.loc[g, p]) for p in fm.columns}).T
    curves.index.name = "period_bp"
    curves.to_csv(out_csv)
    periods = curves.index.values
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.axvspan(*band, color="gold", alpha=0.3, zorder=0)
    ax.axhline(0, color="k", lw=0.6)
    ax.plot(periods, curves.values, color="0.85", lw=0.4, alpha=0.5, zorder=1)
    med = {}
    for k, c in BAR_COLORS.items():
        sub = curves.loc[:, lab[curves.columns] == k]
        if sub.shape[1]:
            med[k] = sub.median(axis=1)
            ax.plot(periods, med[k], color=c, lw=2.6, zorder=3, label=f"{k} (n={sub.shape[1]})")
    ax.set_xlabel("period (bp)")
    ax.set_ylabel("correlation with expression (r)")
    ax.set_title(f"{title}: correlation per period\nthin grey = each of the {curves.shape[1]} groups")
    ax.legend(fontsize=9, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    in_band = (periods >= band[0]) & (periods <= band[1])
    return {f"s_curve_{k}_at_band": v[in_band].mean() for k, v in med.items()}


# ---- coverage ---------------------------------------------------------------------------------
def window_qc(sample_dir):
    """per-gene coverage from the block files: how much of each 10 kb window has reads"""
    body = sample_dir / "body"
    rows = []
    for f in sorted(body.glob("block_*.tsv.gz")):
        a = pd.read_csv(f, sep="\t", header=None, usecols=[2, 3], dtype=np.float32).to_numpy()
        cov, starts = a[:, 0], a[:, 1]
        z = (cov == 0).view(np.int8)
        if z.any():
            e = np.flatnonzero(np.diff(np.concatenate([[0], z, [0]])))
            gap = int((e[1::2] - e[::2]).max())
        else:
            gap = 0
        rows.append({"gene_id": f.name[6:-7], "zero_frac": float((cov == 0).mean()),
                     "mean_cov": float(cov.mean()), "longest_gap": gap, "starts": float(starts.sum())})
    return pd.DataFrame(rows).set_index("gene_id")


def fig_coverage_curves(sample_dir, qc, out, title, zoom=2000):
    """coverage over the whole window and WPS over the first 2 kb, for best / median / worst genes"""
    ranked = qc.sort_values("mean_cov", ascending=False)
    mid = len(ranked) // 2
    picks = [("best", ranked.index[0]), ("median", ranked.index[mid]), ("worst", ranked.index[-1])]
    fig, axes = plt.subplots(len(picks), 2, figsize=(13, 2.7 * len(picks)))
    for i, (kind, gene) in enumerate(picks):
        a = pd.read_csv(sample_dir / "body" / f"block_{gene}.tsv.gz", sep="\t", header=None,
                        usecols=[2, 4], dtype=np.float64).to_numpy()
        cov, wps = a[:, 0], a[:, 1]
        m = qc.loc[gene]
        axes[i, 0].fill_between(np.arange(len(cov)), cov, color="tab:blue", lw=0)
        axes[i, 0].set_ylabel(f"{kind}\n{gene}", fontsize=8)
        axes[i, 0].set_title(f"coverage | mean {m['mean_cov']:.2f}x, {m['zero_frac']:.0%} empty, "
                             f"gap {m['longest_gap']:.0f} bp", fontsize=9)
        axes[i, 1].plot(np.arange(zoom), wps[:zoom], color="#333333", lw=1)
        axes[i, 1].set_title(f"WPS, first {zoom // 1000} kb", fontsize=9)
    axes[-1, 0].set_xlabel("position in window (bp)")
    axes[-1, 1].set_xlabel("bp from TSS")
    fig.suptitle(f"{title}: coverage and WPS", y=1.01)
    plt.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_exclusion_scan(qc, ref, score, lab, shared, mode, var_n, out_png, out_csv, title,
                       drops=(0.0, 0.01, 0.05, 0.10, 0.25, 0.50)):
    """how the correlations move as the least-covered genes are excluded, for one gene set"""
    cov_rank = qc["mean_cov"].rank(method="first", ascending=False)   # 1 = best covered
    rows = []
    for q in drops:
        keep = qc.index[cov_rank <= len(qc) * (1 - q)]
        # re-derive the gene set from the survivors, so var3k is picked among the genes that remain
        g = pick_genes(mode, ref, shared.intersection(keep), var_n)
        r = ref.loc[g].corrwith(score.loc[g]).sort_values()
        med = {k: r[lab[r.index] == k].median() for k in BAR_COLORS}
        rows.append({"dropped": q, "genes": len(g), "strongest_r": r.min(), "top_group": r.index[0],
                     "immune_median": med["immune"], "other_median": med["other"],
                     "control_median": med["negative control"],
                     "immune_minus_other": med["immune"] - med["other"],
                     "immune_in_top20": int((lab[r.index[:20]] == "immune").sum())})
    t = pd.DataFrame(rows).set_index("dropped")
    t.to_csv(out_csv)

    x = np.arange(len(t))
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.6))
    axes[0].plot(x, t["strongest_r"], "o-", color="tab:blue")
    axes[0].set_title("strongest correlation"); axes[0].set_ylabel("r")
    axes[1].plot(x, t["immune_minus_other"], "o-", color="tab:red")
    axes[1].set_title("immune vs other separation"); axes[1].set_ylabel("difference in median r")
    axes[2].plot(x, t["genes"], "o-", color="0.5")
    axes[2].set_title("genes kept"); axes[2].set_ylabel("genes")
    for ax in axes:
        ax.set_xticks(x); ax.set_xticklabels([f"{q:.0%}" for q in t.index], fontsize=9)
        ax.set_xlabel("least-covered genes excluded")
    fig.suptitle(f"{title}: effect of excluding low-coverage genes", y=1.04)
    plt.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    return t


# ---- comparison figures ----------------------------------------------------------------------
def fig_scatter(res, cat, out, title):
    (a, ra), (b, rb) = list(res.items())[:2]
    draw = {"other": 0, "blood/immune": 1, **{k: 2 + i for i, k in enumerate(CONTROL_CATS)}}
    order = cat.sort_values(key=lambda c: c.map(draw)).index
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    for k in sorted(CAT_COLORS, key=lambda k: draw[k]):
        sel = (cat[order] == k).values
        ax.scatter(ra[order][sel], rb[order][sel], s=55 if k in CAT_MARKERS else 10,
                   marker=CAT_MARKERS.get(k, "o"), c=CAT_COLORS[k],
                   alpha=0.9 if k in CAT_MARKERS else 0.6,
                   edgecolors="white" if k in CAT_MARKERS else "none", linewidths=0.6,
                   zorder=6 if k in CAT_MARKERS else 2, label=k)
    lo = min(ra.min(), rb.min()); hi = max(ra.max(), rb.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax.set_xlabel(f"{a} r"); ax.set_ylabel(f"{b} r")
    ax.set_title(f"{title}  (Spearman ρ = {ra.corr(rb, method='spearman'):.2f})")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_rank_shift(res, celltype, comp, cat, out, title, n=15):
    (a, ra), (b, rb) = list(res.items())[:2]
    def by_ct(r):
        agg = r.groupby(celltype[r.index]).mean().sort_values()
        return agg
    ca, cb = by_ct(ra), by_ct(rb)
    common = ca.index.intersection(cb.index)
    shift = (ca[common].rank() - cb[common].rank()).sort_values()
    sel = pd.concat([shift.head(n), shift.tail(n)])
    rep = {ct: ra.index[celltype[ra.index] == ct][0] for ct in sel.index}
    colors = [BAR_COLORS["negative control"] if cat[rep[ct]] in CONTROL_CATS
              else BAR_COLORS["immune" if comp[rep[ct]] == "immune" else "other"] for ct in sel.index]
    fig, ax = plt.subplots(figsize=(8, 0.28 * len(sel) + 1.8))
    ax.barh(np.arange(len(sel)), sel.values, color=colors, alpha=0.85)
    ax.set_yticks(np.arange(len(sel)))
    ax.set_yticklabels([f"{ct}  (#{int(ca[common].rank()[ct])} → #{int(cb[common].rank()[ct])})"
                        for ct in sel.index], fontsize=8)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel(f"rank change (positive = stronger in {b}), out of {len(common)} cell types")
    ax.set_title(f"{title}: largest rank changes {a} → {b}\nexploratory: one sample per condition")
    ax.legend(handles=[Patch(color=c, label=k) for k, c in BAR_COLORS.items()],
              fontsize=8, loc="lower right")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_tissue_order(res, tissue, out, title):
    (a, ra), (b, rb) = list(res.items())[:2]
    ranks = {s: r.rank() for s, r in ((a, ra), (b, rb))}
    tis = pd.DataFrame({s: tissue.groupby(tissue).apply(lambda t: ranks[s][t.index].median())
                        for s in (a, b)})
    order = {s: tis[s].rank() for s in (a, b)}
    shift = order[a] - order[b]
    rho = tis[a].corr(tis[b], method="spearman")
    big = shift.abs().nlargest(5).index
    fig, ax = plt.subplots(figsize=(6.5, 0.32 * len(tis) + 1.8))
    for t in tis.index:
        color, lw = (("tab:red" if shift[t] > 0 else "tab:blue"), 1.8) if t in big else ("0.75", 1.0)
        ax.plot([0, 1], [order[a][t], order[b][t]], color=color, lw=lw, zorder=3 if t in big else 1)
    for x, s, ha in ((0, a, "right"), (1, b, "left")):
        for t in tis.index:
            ax.text(x + (-0.02 if ha == "right" else 0.02), order[s][t], t.replace("_", " "),
                    ha=ha, va="center", fontsize=7, color="k" if t in big else "0.45")
    ax.set_xlim(-0.5, 1.5); ax.set_xticks([0, 1]); ax.set_xticklabels([a, b])
    ax.invert_yaxis(); ax.set_ylabel("tissue rank (1 = strongest median)")
    ax.set_title(f"{title}: tissue order, {a} vs {b} (ρ = {rho:.2f})\ncoloured = 5 largest moves")
    ax.spines[["left"]].set_visible(False)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def _box(ax, vals, pos, width, color, horizontal=False):
    "one box; matplotlib >= 3.11 wants orientation=, older versions want vert="
    kw = dict(positions=[pos], widths=width, patch_artist=True, showfliers=False)
    try:
        bp = ax.boxplot(vals, orientation="horizontal" if horizontal else "vertical", **kw)
    except TypeError:
        bp = ax.boxplot(vals, vert=not horizontal, **kw)
    bp["boxes"][0].set(facecolor=color, alpha=.45, edgecolor=color)
    bp["medians"][0].set(color=color, lw=2)
    for w in bp["whiskers"] + bp["caps"]:
        w.set(color=color, alpha=.8)
    return bp


def fig_paired_boxes(res, keys, members, out, title, ylab, horizontal=False, figsize=(16, 5.5)):
    """one box per sample per key, side by side. Uses RANK, not r: correlations are not
    comparable between samples of different depth, ranks are."""
    ranks = {s: r.rank() for s, r in res.items()}
    rng = np.random.default_rng(0)
    n = len(res)
    fig, ax = plt.subplots(figsize=figsize)
    for i, k in enumerate(keys):
        groups = members(k)
        for j, s in enumerate(res):
            pos = i + (j - (n - 1) / 2) * .38
            vals = ranks[s][groups]
            _box(ax, vals, pos, .34, SAMPLE_COLORS[j], horizontal)
            jit = rng.normal(pos, .05, len(vals))
            xy = (vals.values, jit) if horizontal else (jit, vals.values)
            ax.scatter(*xy, s=7, color=SAMPLE_COLORS[j], alpha=.55, zorder=3)
    ticks = [f"{k.replace('_', ' ')} ({len(members(k))})" for k in keys]
    n_groups = len(list(res.values())[0])
    if horizontal:
        ax.set_yticks(range(len(keys))); ax.set_yticklabels(ticks, fontsize=9)
        ax.invert_yaxis(); ax.set_xlabel(ylab)
    else:
        ax.set_xticks(range(len(keys)))
        ax.set_xticklabels(ticks, rotation=45, ha="right", fontsize=9)
        ax.invert_yaxis(); ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.legend(handles=[Patch(color=SAMPLE_COLORS[j], alpha=.6, label=s) for j, s in enumerate(res)],
              fontsize=9, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_spread_compare(res, lab, out, title):
    fig, axes = plt.subplots(1, len(res), figsize=(4.6 * len(res), 4.6), sharey=True, squeeze=False)
    for ax, (s, r) in zip(axes[0], res.items()):
        rng = np.random.default_rng(0)
        for i, k in enumerate(BAR_COLORS):
            vals = r[lab[r.index] == k].values
            if not len(vals):
                continue
            bp = ax.boxplot(vals, positions=[i], widths=0.55, patch_artist=True, showfliers=False)
            bp["boxes"][0].set(facecolor=BAR_COLORS[k], alpha=0.35)
            bp["medians"][0].set(color=BAR_COLORS[k], lw=2)
            ax.scatter(rng.normal(i, 0.07, len(vals)), vals, s=6, color=BAR_COLORS[k], alpha=0.4, zorder=3)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(range(len(BAR_COLORS)))
        ax.set_xticklabels(list(BAR_COLORS), fontsize=8, rotation=20, ha="right")
        ax.set_title(s)
    axes[0][0].set_ylabel("correlation with expression (r)")
    fig.suptitle(title)
    plt.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


# ---- main ------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample_dir", type=Path, nargs="*", help="results/<sample> folders")
    ap.add_argument("--all", action="store_true", help="use every sample folder under results/")
    ap.add_argument("-o", "--output", type=Path, default=HERE / "analysis")
    ap.add_argument("--results", type=Path, default=WORKFLOW / "results")
    ap.add_argument("--ref", type=Path, default=WORKFLOW / "resources/avg_expr_by_celltype_tissue.csv")
    ap.add_argument("--compartments", type=Path, default=WORKFLOW / "resources/celltype_compartment.csv")
    ap.add_argument("--band", type=float, nargs=2, default=(193, 199), metavar=("LOW", "HIGH"))
    ap.add_argument("--var-genes", type=int, default=3000, help="size of the variable-gene set")
    ap.add_argument("--skip-coverage", action="store_true",
                    help="do not read results/<sample>/body (faster; no coverage section)")
    ap.add_argument("--max-zero-frac", type=float, default=0.95,
                    help="drop genes whose window is more than this fraction empty (default 0.95)")
    ap.add_argument("--top-celltypes", type=int, default=25,
                    help="cell types shown in the paired comparison figure (default 25)")
    ap.add_argument("--compare", nargs=2, metavar=("A", "B"), default=None,
                    help="sample pair for the comparison figures (default: first two given)")
    args = ap.parse_args()

    dirs = list(args.sample_dir)
    if args.all:
        dirs += [d for d in sorted(args.results.iterdir())
                 if d.is_dir() and (d / "proj").exists() and d not in dirs]
    if not dirs:
        raise SystemExit("no sample folders given; pass paths or --all")

    ref, comp, cat, lab, tissue, celltype, origin = build_reference(args.ref, args.compartments)
    print(f"reference: {ref.shape[0]:,} genes x {ref.shape[1]} groups")
    print(f"lymphoid tissues: {', '.join(sorted(LYMPHOID))}")

    rows, per_mode = [], {}
    for d in dirs:
        sample = d.name
        fm = fft_matrix(d)
        band_cols = [c for c in fm.columns if args.band[0] <= c <= args.band[1]]
        score = fm[band_cols].mean(axis=1)
        shared = ref.index.intersection(score.index)
        print(f"\n{sample}: {fm.shape[0]:,} genes, {len(shared):,} shared with the reference")

        qc = None
        if not args.skip_coverage and (d / "body").exists():
            cov_out = args.output / sample / "coverage"
            cov_out.mkdir(parents=True, exist_ok=True)
            qc = window_qc(d)
            qc.to_csv(cov_out / "window_qc.csv")
            fig_coverage_curves(d, qc, cov_out / "curves.png", sample)
            empty = qc.index[qc["zero_frac"] > args.max_zero_frac]
            shared = shared.difference(empty)
            print(f"  coverage: median {qc['mean_cov'].median():.2f}x, "
                  f"{qc['zero_frac'].median():.0%} of bases empty, "
                  f"{len(empty)} genes over {args.max_zero_frac:.0%} empty dropped "
                  f"-> {len(shared):,} genes")
        elif not args.skip_coverage:
            print(f"  (no {d / 'body'}, skipping coverage)")

        for mode, genes in gene_sets(ref, shared, args.var_genes).items():
            out = args.output / sample / mode
            out.mkdir(parents=True, exist_ok=True)
            r = ref.loc[genes].corrwith(score.loc[genes]).sort_values()
            title = f"{sample} [{mode}]"
            pd.DataFrame({"correlation": r, "category": cat[r.index], "label": lab[r.index],
                          "compartment": comp[r.index]}).to_csv(out / "correlations.csv")

            fig_top_bottom(r, cat, comp, out / "top_bottom.png", title)
            fig_spread(r, lab, out / "spread.png", title)
            rank = r.rank()
            fig_by_compartment(rank, comp, cat, out / "by_compartment.png", title)
            fig_by_tissue(rank, tissue, cat, out / "by_tissue.png", title)
            agg = fig_by_celltype_origin(rank, r, celltype, tissue, out / "by_celltype_origin.png", title)
            agg.to_csv(out / "by_celltype.csv")
            sc = fig_s_curves(fm, ref, genes, lab, args.band,
                              out / "s_curves.png", out / "s_curves.csv.gz", title)
            if qc is not None:
                scan = fig_exclusion_scan(qc, ref, score, lab, shared, mode, args.var_genes,
                                          out / "exclusion_scan.png", out / "exclusion_scan.csv", title)
                print("    exclusion scan: " + " | ".join(
                    f"{q:.0%}:{row.strongest_r:+.4f}" for q, row in scan.iterrows()))

            med = {k: r[lab[r.index] == k].median() for k in BAR_COLORS}
            ori = {"celltype_median_lymphoid": agg.loc[agg["lymphoid_fraction"] > 0.5, "mean"].median(),
                   "celltype_median_solid": agg.loc[agg["lymphoid_fraction"] == 0, "mean"].median()}
            rows.append({"sample": sample, "mode": mode, "genes": len(genes),
                         "median_cov": None if qc is None else qc["mean_cov"].median(),
                         "median_zero_frac": None if qc is None else qc["zero_frac"].median(),
                         "strongest_r": r.min(), "top_group": r.index[0],
                         "immune_median": med["immune"], "other_median": med["other"],
                         "control_median": med["negative control"],
                         "immune_minus_other": med["immune"] - med["other"],
                         "immune_in_top20": int((lab[r.index[:20]] == "immune").sum()),
                         **ori, **sc})
            per_mode.setdefault(mode, {})[sample] = r
            print(f"  {mode:9s} {len(genes):6,} genes | strongest {r.min():+.4f} ({r.index[0]}) | "
                  f"immune-other {med['immune'] - med['other']:+.4f} | "
                  f"immune in top 20: {int((lab[r.index[:20]] == 'immune').sum())}")

    summary = pd.DataFrame(rows).set_index(["sample", "mode"])
    args.output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output / "summary.csv", float_format="%.5f")

    # ---- comparison ---------------------------------------------------------------------
    for mode, res in per_mode.items():
        if len(res) < 2:
            continue
        pair = args.compare if args.compare else list(res)[:2]
        pair = [p for p in pair if p in res]
        if len(pair) < 2:
            print(f"skipping comparison for {mode}: {args.compare} not both present")
            continue
        sub = {p: res[p] for p in pair}
        out = args.output / "comparison" / mode
        out.mkdir(parents=True, exist_ok=True)
        title = f"[{mode}]"
        fig_scatter(sub, cat, out / f"scatter_{pair[0]}_vs_{pair[1]}.png", title)
        fig_rank_shift(sub, celltype, comp, cat, out / "rank_shift.png", title)
        fig_tissue_order(sub, tissue, out / "tissue_order.png", title)
        fig_spread_compare(res, lab, out / "spread_compare.png", f"Correlation spread by label {title}")

        nmax = len(list(sub.values())[0])
        rank_pair = {s: r.rank() for s, r in sub.items()}
        tis_order = (tissue.groupby(tissue)
                     .apply(lambda t: np.mean([rank_pair[s][t.index].median() for s in sub]))
                     .sort_values().index.tolist())
        fig_paired_boxes(sub, tis_order, lambda k: list(tissue.index[tissue == k]),
                         out / "by_tissue_compare.png",
                         f"Rank by tissue, {pair[0]} vs {pair[1]} {title}",
                         f"rank (1 = strongest, {nmax} = weakest)")

        ct_med = (celltype.groupby(celltype)
                  .apply(lambda t: np.mean([rank_pair[s][t.index].median() for s in sub]))
                  .sort_values())
        ct_order = ct_med.index[:args.top_celltypes].tolist()
        fig_paired_boxes(sub, ct_order, lambda k: list(celltype.index[celltype == k]),
                         out / "by_celltype_compare.png",
                         f"Rank by cell type, strongest {len(ct_order)}, {pair[0]} vs {pair[1]} {title}",
                         f"rank (1 = strongest, {nmax} = weakest)",
                         horizontal=True, figsize=(10, .32 * len(ct_order) + 2))
        print(f"\ncomparison [{mode}]: {pair[0]} vs {pair[1]} -> {out}")

    print(f"\n{args.output / 'summary.csv'}")
    print(summary.round(4).to_string())


if __name__ == "__main__":
    main()
