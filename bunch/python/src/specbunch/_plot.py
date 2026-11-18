"""Optional matplotlib figures for NPMLE postestimation (pip install specbunch[plot])."""
import sys
import numpy as np

KINDS = ("both", "pdf", "cdf")
COLORS = {"data": "#b3b3b3", "observed": "#1f3b73", "latent": "#d9822b", "kink": "#555555"}


def plot_kind(plot):
    """Normalise the plot option to "both", "pdf", "cdf" or None (no plot)."""
    if plot is None or plot is False:
        return None
    if plot is True:
        return "both"
    if isinstance(plot, str):
        k = plot.lower()
        if k in ("none", "false", "no", "off"):
            return None
        if k in KINDS:
            return k
    raise ValueError('plot must be "both", "pdf", "cdf", or False/None')


def _interactive():
    # Use pyplot (so the figure displays) only when the user is already in a pyplot or
    # IPython workflow; otherwise build a standalone Figure with no global state or GUI.
    if "matplotlib.pyplot" in sys.modules:
        return True
    ipython = sys.modules.get("IPython")
    return ipython is not None and getattr(ipython, "get_ipython", lambda: None)() is not None


def _view_range(grid, ecdf, kink):
    lo = grid[min(np.searchsorted(ecdf, .001, side="left"), len(grid) - 1)]
    hi = grid[min(np.searchsorted(ecdf, .999, side="left"), len(grid) - 1)]
    lo, hi = min(lo, kink), max(hi, kink)
    pad = .04 * (hi - lo if hi > lo else 1.0)
    return max(lo - pad, grid[0]), min(hi + pad, grid[-1])


def _note(r):
    lines = [r"$\hat B = %.4g$" % r.B, r"$\hat\sigma = %.4g$" % r.sigma]
    if r.raw.get("a", 0) > 0:
        lines.append(r"$\hat a = %.4g$" % r.raw["a"])
    if r.status != "converged":
        lines.append("status: " + r.status)
    return "\n".join(lines)


def draw(r, kind="both", xlabel=None, figure=None):
    """Draw the density and/or CDF view of an NpmleResult; return the Figure."""
    kind = plot_kind(kind) or "both"
    try:
        import matplotlib
    except ImportError as e:
        raise ImportError("plotting needs matplotlib: pip install 'specbunch[plot]'") from e
    raw = r.raw
    grid = raw["curves.grid"]
    kink = float(r.kink)
    xlabel = xlabel or getattr(r, "xlabel", None) or "z"
    size = (10.0, 3.9) if kind == "both" else (5.4, 3.9)
    if figure is None:
        if _interactive():
            import matplotlib.pyplot as plt
            figure = plt.figure(figsize=size)
        else:
            from matplotlib.figure import Figure
            figure = Figure(figsize=size)
    panels = ("pdf", "cdf") if kind == "both" else (kind,)
    axes = figure.subplots(1, len(panels), squeeze=False)[0]
    xlim = _view_range(grid, raw["curves.empirical_cdf"], kink)
    with matplotlib.rc_context({"axes.facecolor": "white", "figure.facecolor": "white", "font.size": 9,
                                "axes.spines.top": False, "axes.spines.right": False,
                                "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6}):
        figure.set_facecolor("white")
        for ax, panel in zip(axes, panels):
            ax.set_facecolor("white")
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            if panel == "pdf":
                lo, hi = raw["histogram.bin_lo"], raw["histogram.bin_hi"]
                width = hi - lo
                dens = raw["histogram.observed_share"] / width
                if np.array_equal(lo[1:], hi[:-1]):  # contiguous bins: one artist, fast
                    ax.stairs(dens, np.append(lo, hi[-1]), fill=True, color=COLORS["data"],
                              linewidth=0, label="Data (histogram)")
                else:
                    ax.bar(lo, dens, width=width, align="edge", color=COLORS["data"], linewidth=0,
                           label="Data (histogram)")
                ax.plot(grid, raw["curves.observed_pdf"], color=COLORS["observed"], lw=1.3,
                        label="Observed fit")
                ax.plot(grid, raw["curves.latent_pdf"], color=COLORS["latent"], lw=1.3, ls="--",
                        label="Latent fit, excl. bunchers")
                ax.set_ylabel("Density")
                ax.set_ylim(0, 1.08 * ax.get_ylim()[1])
                ax.text(.02, .97, _note(r), transform=ax.transAxes, va="top", ha="left", fontsize=8.5,
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85, "pad": 1.5})
                legend_loc = "upper right"
            else:
                ax.plot(grid, raw["curves.empirical_cdf"], color=COLORS["data"], lw=2.6,
                        drawstyle="steps-post", label="Empirical CDF")
                ax.plot(grid, raw["curves.observed_cdf"], color=COLORS["observed"], lw=1.3,
                        label="Observed fit")
                ax.plot(grid, raw["curves.latent_cdf"], color=COLORS["latent"], lw=1.3, ls="--",
                        label=r"Latent fit (jump $\hat B$ at kink)")
                ax.set_ylabel("Cumulative probability")
                ax.set_ylim(0, 1.02)
                legend_loc = "upper left"
            ax.axvline(kink, color=COLORS["kink"], lw=.7, ls=":", label="Kink", zorder=0.5)
            ax.set_xlim(*xlim)
            ax.set_xlabel(xlabel)
            handles, labels = ax.get_legend_handles_labels()
            order = sorted(range(len(labels)), key=lambda i: labels[i] not in ("Data (histogram)", "Empirical CDF"))
            ax.legend([handles[i] for i in order], [labels[i] for i in order], loc=legend_loc,
                      frameon=True, facecolor="white", edgecolor="none", framealpha=.9,
                      fontsize=7.5, handlelength=2.2)
        figure.tight_layout()
    return figure
