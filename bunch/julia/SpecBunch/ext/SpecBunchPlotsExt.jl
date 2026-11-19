"""Plots.jl recipes for `SpecBunch.NpmleResult` and mass `SpecBunch.Result` transform fits; loaded with `using Plots`."""
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

# Transform fit of a mass result: empirical transform, +/-2 SE band, fit over the whole grid,
# shaded selected window. Joint fits add the t*S(t) panel.
const EMPIRICAL = colorant"#4d4d4d"
const BAND = colorant"#d4d4d4"
const FIT = colorant"#1f3b73"
const WINDOW = colorant"#f6d58e"
const EDGE = colorant"#b8860b"
const ZERO = colorant"#999999"

function _transform_titles(raw, refused)
    method = get(raw, "method", "cosine"); friction = get(raw, "friction", "gaussian")
    head = uppercasefirst(method) * " fit, " * friction * " friction"
    refused && return (head * "; refused: " * get(raw, "reason", ""), "")
    est = "\\hat{B} = " * _fmt(raw["B"]) * (friction == "voigt" ?
        ",\\ \\hat{a} = " * _fmt(raw["a"]) * ",\\ \\hat{s}^2 = " * _fmt(raw["s2"]) :
        ",\\ \\hat{\\sigma} = " * _fmt(raw["sigma"]))
    jump = haskey(raw, "jump") ? "\$\\widehat{\\mathrm{jump}} = " * _fmt(raw["jump"]) * "\$" : ""
    ("\$\\mathrm{" * replace(head, " " => "\\ ") * ":}\\ \\ " * est * "\$", jump)
end

@recipe function f(r::SpecBunch.Result; frequency=:scaled)
    r.kind == :mass || throw(ArgumentError("transform plots are available for mass results"))
    xs = something(SpecBunch._transform_xscale(frequency), :scaled)
    tr = r.transform; raw = r.raw
    (isempty(tr) || isempty(tr.t)) && throw(ArgumentError("this mass result has no transform grid"))
    t = xs == :scaled ? tr.t_scaled : tr.t
    sfx = xs == :scaled ? "_scaled" : ""
    lo = get(raw, "window_lo" * sfx, NaN); hi = get(raw, "window_hi" * sfx, NaN)
    haswin = isfinite(lo) && isfinite(hi)
    joint = haskey(tr, :sine_scaled) && length(tr.sine_scaled) == length(t)
    titles = _transform_titles(raw, r.refused)
    panels = joint ? (1, 2) : (1,)
    attr(a, b) = permutedims([p == 1 ? a : b for p in panels])
    layout := (1, length(panels))
    size --> (joint ? (1000, 390) : (560, 390))
    background_color --> :white
    foreground_color_legend --> nothing
    background_color_legend --> RGBA(1, 1, 1, .9)
    grid --> false
    framestyle --> :axes
    legendfontsize --> 7
    guidefontsize --> 9
    tickfontsize --> 8
    titlefontsize --> 9
    titlelocation --> :left
    legend --> :topright
    xlims --> (t[1], t[end])
    xguide --> (xs == :scaled ? "\$\\mathrm{Scaled\\ frequency}\\ t\\,\\hat{\\sigma}_{\\mathrm{pilot}}\$" : "\$\\mathrm{Frequency}\\ t\$")
    yguide --> attr("\$\\mathrm{Cosine\\ transform}\\ \\hat{C}(t)\$", "\$\\mathrm{Scaled\\ sine\\ transform}\\ t\\,\\hat{S}(t)\$")
    title --> attr(titles...)
    left_margin --> 4Plots.mm
    bottom_margin --> 4Plots.mm
    for i in panels
        y, se, fit = i == 1 ? (tr.cosine, tr.cosine_se, tr.cosine_fit) : (tr.sine_scaled, tr.sine_scaled_se, tr.sine_fit)
        emp = i == 1 ? "\$\\mathrm{Empirical}\\ \\hat{C}(t)\$" : "\$\\mathrm{Empirical}\\ t\\,\\hat{S}(t)\$"
        fl = i == 1 ? "\$\\mathrm{Fit}\\ \\hat{B}\\,\\hat{\\varphi}(t)\$" : "\$\\mathrm{Fit}\\ \\widehat{\\mathrm{jump}}\\,\\hat{\\varphi}(t)\$"
        if haswin
            @series begin
                subplot := i; seriestype := :vspan; label := ""
                fillcolor := WINDOW; fillalpha := .45; linewidth := 0; linealpha := 0
                [lo, hi]
            end
            @series begin
                subplot := i; seriestype := :vline; label := ""; linecolor := EDGE; linewidth := .7; linestyle := :dash
                [lo, hi]
            end
        end
        @series begin
            subplot := i; seriestype := :hline; label := ""; linecolor := ZERO; linewidth := .5
            [0.0]
        end
        @series begin
            subplot := i; label := ""; seriestype := :path
            fillrange := y .+ 2 .* se; fillcolor := BAND; fillalpha := 1; linewidth := 0; linealpha := 0
            t, y .- 2 .* se
        end
        @series begin
            subplot := i; label := emp; linecolor := EMPIRICAL; linewidth := .8
            markershape := :circle; markersize := 1.2; markercolor := EMPIRICAL; markerstrokewidth := 0
            t, y
        end
        @series begin  # legend entry for the band, listed after the empirical curve
            subplot := i; seriestype := :shape; label := "\$\\pm 2\\ \\mathrm{SE}\$"
            fillcolor := BAND; linewidth := 0; linealpha := 0
            [NaN], [NaN]
        end
        if any(isfinite, fit)
            @series begin
                subplot := i; label := fl; linecolor := FIT; linewidth := 1.5
                t, fit
            end
        end
        if haswin  # legend entry only, so the window is listed last but drawn first
            @series begin
                subplot := i; seriestype := :shape; label := "Selected window [" * _fmt(lo) * ", " * _fmt(hi) * "]"
                fillcolor := WINDOW; fillalpha := .45; linewidth := 0; linealpha := 0
                [NaN], [NaN]
            end
        end
    end
end

function __init__()
    SpecBunch._PLOT_HOOK[] = (r, k) -> r isa SpecBunch.NpmleResult ? Plots.plot(r; kind=k) : Plots.plot(r; frequency=k)
end

end
