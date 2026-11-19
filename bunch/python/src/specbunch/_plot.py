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


# ---------------------------------------------------------------------------
# Transform-fit figure for mass results.
TCOLORS = {"empirical": "#4d4d4d", "band": "#d4d4d4", "fit": "#1f3b73", "window": "#f6d58e",
           "edge": "#b8860b", "zero": "#999999"}
XSCALES = ("scaled", "raw")


def transform_xscale(plot):
    """Normalise mass's plot option to "scaled", "raw" or None (no plot)."""
    if plot is None or plot is False:
        return None
    if plot is True:
        return "scaled"
    if isinstance(plot, str):
        k = plot.lower()
        if k in ("none", "false", "no", "off"):
            return None
        if k in XSCALES:
            return k
    raise ValueError('plot must be True, "scaled", "raw", or False/None')


def _g(x):
    return "%.4g" % x


def _transform_note(r):
    """Panel titles: method, friction and estimates (cosine panel); jump (sine panel)."""
    raw = r.raw
    method, friction = raw.get("method", "") or "cosine", raw.get("friction", "") or "gaussian"
    head = "%s fit, %s friction" % (method.capitalize(), friction)
    if r.refused:
        return head + "; refused: " + (r.reason or "status %d" % r.status), ""
    est = [r"$\hat B = %s$" % _g(raw["B"])]
    if friction == "voigt":
        est += [r"$\hat a = %s$" % _g(raw["a"]), r"$\hat s^2 = %s$" % _g(raw["s2"])]
    else:
        est.append(r"$\hat\sigma = %s$" % _g(raw["sigma"]))
    jump = r"$\widehat{\mathrm{jump}} = %s$" % _g(raw["jump"]) if "jump" in raw else ""
    return head + ":  " + ",  ".join(est), jump


def draw_transform(r, frequency="scaled", figure=None):
    """Draw the empirical-transform fit of a mass Result; return the Figure (None if no grid)."""
    xscale = transform_xscale(frequency) or "scaled"
    try:
        import matplotlib
    except ImportError as e:
        raise ImportError("plotting needs matplotlib: pip install 'specbunch[plot]'") from e
    raw = r.raw
    t = raw.get("transform.t_scaled" if xscale == "scaled" else "transform.t", np.empty(0))
    if len(t) == 0:
        return None
    joint = len(raw.get("transform.sine_scaled", ())) == len(t)
    suffix = "_scaled" if xscale == "scaled" else ""
    lo, hi = raw.get("window_lo" + suffix, np.nan), raw.get("window_hi" + suffix, np.nan)
    xlabel = (r"Scaled frequency $t\,\hat\sigma_{\mathrm{pilot}}$" if xscale == "scaled"
              else r"Frequency $t$")
    size = (10.0, 3.9) if joint else (5.6, 3.9)
    if figure is None:
        if _interactive():
            import matplotlib.pyplot as plt
            figure = plt.figure(figsize=size)
        else:
            from matplotlib.figure import Figure
            figure = Figure(figsize=size)
    panels = [("transform.cosine", "transform.cosine_se", "transform.cosine_fit",
               r"Cosine transform $\hat C(t)$", r"Empirical $\hat C(t)$", r"Fit $\hat B\,\hat\varphi(t)$")]
    if joint:
        panels.append(("transform.sine_scaled", "transform.sine_scaled_se", "transform.sine_fit",
                       r"Scaled sine transform $t\,\hat S(t)$", r"Empirical $t\,\hat S(t)$",
                       r"Fit $\widehat{\mathrm{jump}}\,\hat\varphi(t)$"))
    axes = figure.subplots(1, len(panels), squeeze=False)[0]
    band_label = r"$\pm 2$ SE"
    window_label = "Selected window [%s, %s]" % (_g(lo), _g(hi))
    notes = _transform_note(r)
    with matplotlib.rc_context({"axes.facecolor": "white", "figure.facecolor": "white", "font.size": 9,
                                "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6}):
        figure.set_facecolor("white")
        for i, (ax, (yk, sk, fk, ylabel, elabel, flabel)) in enumerate(zip(axes, panels)):
            ax.set_facecolor("white")
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            y, se, fit = raw[yk], raw[sk], raw[fk]
            if np.isfinite(lo) and np.isfinite(hi):
                ax.axvspan(lo, hi, color=TCOLORS["window"], alpha=.35, lw=0, zorder=0, label=window_label)
                for edge in (lo, hi):
                    ax.axvline(edge, color=TCOLORS["edge"], lw=.7, ls="--", zorder=.6)
            ax.axhline(0, color=TCOLORS["zero"], lw=.5, zorder=.5)
            ax.fill_between(t, y - 2 * se, y + 2 * se, color=TCOLORS["band"], lw=0, zorder=1, label=band_label)
            ax.plot(t, y, color=TCOLORS["empirical"], lw=.7, marker="o", ms=1.3, zorder=2, label=elabel)
            if np.isfinite(fit).any():
                ax.plot(t, fit, color=TCOLORS["fit"], lw=1.4, zorder=3, label=flabel)
            ax.set_xlim(t[0], t[-1])
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(notes[i], loc="left", fontsize=8.5)
            handles, labels = ax.get_legend_handles_labels()
            order = [labels.index(k) for k in (elabel, band_label, flabel, window_label) if k in labels]
            ax.legend([handles[k] for k in order], [labels[k] for k in order],
                      loc="upper right" if i == 0 else "best", frameon=True, facecolor="white",
                      edgecolor="none", framealpha=.9, fontsize=7.5, handlelength=2.2)
        figure.tight_layout()
    return figure
