"""Plots.jl recipe for `SpecBunch.NpmleResult`; loaded automatically with `using Plots`."""
module SpecBunchPlotsExt

using SpecBunch
using Plots

const DATA = colorant"#b3b3b3"
const OBSERVED = colorant"#1f3b73"
const LATENT = colorant"#d9822b"
const KINK = colorant"#555555"

_fmt(x) = string(round(x; sigdigits=4))

function _view_range(grid, ecdf, kink)
    at(p) = grid[something(findfirst(>=(p), ecdf), length(grid))]
    lo = min(at(.001), kink); hi = max(at(.999), kink)
    pad = .04 * (hi > lo ? hi - lo : 1.0)
    (max(lo - pad, grid[1]), min(hi + pad, grid[end]))
end

# One polygon for all bars: fast, and honest for non-contiguous user bins.
function _bars(lo, hi, d)
    x = Float64[]; y = Float64[]
    for i in eachindex(lo)
        append!(x, (lo[i], lo[i], hi[i], hi[i])); append!(y, (0., d[i], d[i], 0.))
    end
    x, y
end

@recipe function f(r::SpecBunch.NpmleResult; kind=:both)
    k = something(SpecBunch._plot_kind(kind), :both)
    panels = k == :both ? (:pdf, :cdf) : (k,)
    cv = r.curves; h = r.histogram; raw = r.raw
    kink = raw["kink"]; grid = cv.grid
    xl = _view_range(grid, cv.empirical_cdf, kink)
    dens = h.observed_share ./ (h.bin_hi .- h.bin_lo)
    inview = (h.bin_hi .>= xl[1]) .& (h.bin_lo .<= xl[2])
    top = 1.08 * maximum(vcat(dens[inview], cv.observed_pdf, cv.latent_pdf))
    # GR renders strings wrapped in $...$ as LaTeX: the note sits as a left-aligned panel title.
    note = "\$\\hat{B} = " * _fmt(raw["B"]) * ",\\ \\hat{\\sigma} = " * _fmt(raw["sigma"])
    get(raw, "a", 0.0) > 0 && (note *= ",\\ \\hat{a} = " * _fmt(raw["a"]))
    note *= "\$"
    raw["status"] == "converged" || (note = "status: " * raw["status"])
    attr(pdfv, cdfv) = permutedims([p == :pdf ? pdfv : cdfv for p in panels])
    layout := (1, length(panels))
    size --> (length(panels) == 2 ? (1000, 390) : (540, 390))
    background_color --> :white
    foreground_color_legend --> nothing
    background_color_legend --> RGBA(1, 1, 1, .85)
    grid --> false
    framestyle --> :axes
    legendfontsize --> 7
    guidefontsize --> 9
    tickfontsize --> 8
    titlefontsize --> 9
    titlelocation --> :left
    xguide --> "z"
    xlims --> xl
    yguide --> attr("Density", "Cumulative probability")
    ylims --> attr((0, top), (0, 1.02))
    legend --> attr(:topright, :topleft)
    title --> attr(note, "")
    left_margin --> 4Plots.mm
    bottom_margin --> 4Plots.mm
    for (i, p) in enumerate(panels)
        if p == :pdf
            bx, by = _bars(h.bin_lo, h.bin_hi, dens)
            @series begin
                subplot := i; seriestype := :shape; label := "Data (histogram)"
                fillcolor := DATA; linecolor := DATA; linewidth := 0
                bx, by
            end
            @series begin
                subplot := i; label := "Observed fit"; linecolor := OBSERVED; linewidth := 1.3
                grid, cv.observed_pdf
            end
            @series begin
                subplot := i; label := "Latent fit, excl. bunchers"; linecolor := LATENT
                linewidth := 1.3; linestyle := :dash
                grid, cv.latent_pdf
            end
        else
            @series begin
                subplot := i; seriestype := :steppost; label := "Empirical CDF"
                linecolor := DATA; linewidth := 2.6
                grid, cv.empirical_cdf
            end
            @series begin
                subplot := i; label := "Observed fit"; linecolor := OBSERVED; linewidth := 1.3
                grid, cv.observed_cdf
            end
            @series begin
                subplot := i; label := "Latent fit (jump of size B at kink)"; linecolor := LATENT
                linewidth := 1.3; linestyle := :dash
                grid, cv.latent_cdf
            end
        end
        @series begin
            subplot := i; seriestype := :vline; label := "Kink"
            linecolor := KINK; linewidth := .7; linestyle := :dot
            [kink]
        end
    end
end

function __init__()
    SpecBunch._PLOT_HOOK[] = (r, kind) -> Plots.plot(r; kind=kind)
end

end
