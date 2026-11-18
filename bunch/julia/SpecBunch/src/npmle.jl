struct NpmleOpts
 struct_size::UInt32;abi_version::UInt32
 threads::Int32;core_intervals::Int32;grid_points::Int32;exact_threshold::Int32
 max_iterations::Int32;lambda_points::Int32;refinement_points::Int32;reserved::Int32
 B::Float64;sigma::Float64;a::Float64;kink::Float64;histogram_width::Float64
 lambda_min::Float64;lambda_max::Float64;quadrature_tolerance::Float64;optimality_tolerance::Float64
end
"""Fixed-parameter spline postestimation with histogram and curve tables.

`figure` holds the Plots.jl figure drawn by `npmle` (nothing unless Plots is loaded).
"""
mutable struct NpmleResult
 status_code::Int;raw::Dict{String,Any};figure::Any
 NpmleResult(status_code,raw)=new(status_code,raw,nothing)
end
function Base.getproperty(r::NpmleResult,s::Symbol)
 s in (:status_code,:raw,:figure) && return getfield(r,s)
 raw=getfield(r,:raw)
 if s in (:histogram,:curves,:lambda_path)
  prefix=String(s)*".";names=sort!([k for k in keys(raw) if startswith(k,prefix)])
  return (;(Symbol(k[length(prefix)+1:end])=>raw[k] for k in names)...)
 end
 haskey(raw,String(s)) && return raw[String(s)]
 haskey(raw,"curves."*String(s)) && return raw["curves."*String(s)]
 throw(ArgumentError("unknown NPMLE result property $s"))
end
Base.propertynames(r::NpmleResult)=(:status_code,:raw,:figure,:histogram,:curves,:lambda_path,(Symbol(k) for k in keys(getfield(r,:raw)) if !occursin('.',k))...)
# Set by the Plots.jl package extension (ext/SpecBunchPlotsExt.jl) when Plots is loaded.
const _PLOT_HOOK=Ref{Any}(nothing)
function _plot_kind(plot)
 (plot===nothing || plot===false) && return nothing
 plot===true && return :both
 if plot isa Union{Symbol,AbstractString}
  k=Symbol(lowercase(String(plot)))
  k in (:none,:false,:no,:off) && return nothing
  k in (:both,:pdf,:cdf) && return k
 end
 throw(ArgumentError("plot must be :both, :pdf, :cdf, true or false/nothing"))
end
Base.show(io::IO,r::NpmleResult)=print(io,"NpmleResult(status=",r.status,", lambda=",r.lambda,", edf=",r.edf,")")
"""
    npmle(fit,z; weights=nothing,bins=nothing,...)

Fit nonnegative clamped cubic splines on each side with fixed atom mass and
friction. BIC selects the roughness penalty. Returns histogram shares and the
latent, observed and empirical curves; see docs/REFERENCE.md.

Main options: `weights` supplies positive observation weights; `bins` supplies
output edges without restricting the sample.

`plot=:both` (default), `:pdf` or `:cdf` draws the fitted latent and observed
distributions against the data into `result.figure` when Plots.jl is loaded
(`using Plots`), displaying it in interactive sessions; without Plots the
request is a no-op with a one-time `@info`. `plot=false` skips it. With Plots
loaded, `plot(result; kind=:cdf)` redraws via the package's recipe.

Advanced (solver settings): `threads`, `core_intervals`, `grid_points`,
`exact_threshold`, `max_iterations`, `lambda_min`, `lambda_max`,
`lambda_points`, `refinement_points`, `histogram_width`,
`quadrature_tolerance` and `optimality_tolerance`; see docs/REFERENCE.md.
"""
function npmle(fit::Result,z;weights=nothing,bins=nothing,threads=nothing,core_intervals=12,
               grid_points=1001,exact_threshold=1000,max_iterations=200,lambda_min=.01,
               lambda_max=10000.,lambda_points=15,refinement_points=9,histogram_width=nothing,
               quadrature_tolerance=1e-9,optimality_tolerance=1e-7,plot=:both)
 fit.kind==:mass && !fit.refused || throw(ArgumentError("npmle requires a coherent mass result"))
 kind=_plot_kind(plot)
 z,w=_inputs(z,weights);edges=bins===nothing ? nothing : _vector(bins)
 a=fit.raw["friction"]=="voigt" ? fit.raw["a"] : 0.
 o=NpmleOpts(UInt32(sizeof(NpmleOpts)),1,_threads(threads),Int32(core_intervals),Int32(grid_points),Int32(exact_threshold),
             Int32(max_iterations),Int32(lambda_points),Int32(refinement_points),0,fit.B,fit.sigma,a,fit.kink,
             histogram_width===nothing ? 0. : Float64(histogram_width),Float64(lambda_min),Float64(lambda_max),Float64(quadrature_tolerance),Float64(optimality_tolerance))
 out=Ref{Ptr{Cvoid}}(C_NULL)
 result=GC.@preserve z w edges begin
  d=MassData(length(z),pointer(z),_ptr(w))
  rc=ccall((:sb_npmle,libspecbunch),Cint,(Ref{MassData},Ref{NpmleOpts},Ptr{Cdouble},Csize_t,Ref{Ptr{Cvoid}}),d,o,_ptr(edges),edges===nothing ? 0 : length(edges),out)
  try
   rc==0 || throw(CoreError(rc,unsafe_string(ccall((:sb_last_error,libspecbunch),Cstring,()))))
   status,raw=_extract(out[]);NpmleResult(status,raw)
  finally
   out[]==C_NULL || ccall((:sb_result_free,libspecbunch),Cvoid,(Ptr{Cvoid},),out[])
  end
 end
 if kind!==nothing
  if _PLOT_HOOK[]===nothing
   @info "SpecBunch.npmle: load Plots.jl (`using Plots`) to draw the fitted distributions; pass plot=false to skip." maxlog=1
  else
   result.figure=_PLOT_HOOK[](result,kind)
   isinteractive() && display(result.figure)
  end
 end
 result
end
