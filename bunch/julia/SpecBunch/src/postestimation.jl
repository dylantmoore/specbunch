struct NonbuncherOpts
 struct_size::UInt32;abi_version::UInt32;threads::Int32;degree::Int32
 B::Float64;sigma::Float64;a::Float64;kink::Float64;mu_b::Float64;se_mu_b::Float64;bandwidth::Float64;bin_width::Float64
end
struct EtiOpts
 struct_size::UInt32;abi_version::UInt32;threads::Int32;integration_intervals::Int32
 B::Float64;sigma::Float64;a::Float64;kink::Float64;window_lo::Float64;window_hi::Float64;net_below::Float64;net_above::Float64
end
"""Core postestimation outputs, with a local_fit NamedTuple of vectors when available."""
struct PostestimationResult
 status_code::Int;raw::Dict{String,Any}
end
function Base.getproperty(r::PostestimationResult,s::Symbol)
 s in (:status_code,:raw) && return getfield(r,s)
 raw=getfield(r,:raw)
 s==:refused && return getfield(r,:status_code)!=0
 s==:flags && return split(raw["flags"])
 if s in (:local_fit,:saez,:continuous)
  prefix=String(s)*".";names=sort!([k for k in keys(raw) if startswith(k,prefix)])
  return (;(Symbol(k[length(prefix)+1:end])=>raw[k] for k in names)...)
 end
 haskey(raw,String(s)) && return raw[String(s)]
 throw(ArgumentError("unknown postestimation property $s"))
end
Base.propertynames(r::PostestimationResult)=(:status_code,:raw,:refused,:flags,:local_fit,:saez,:continuous,(Symbol(k) for k in keys(getfield(r,:raw)) if !occursin('.',k))...)
Base.show(io::IO,r::PostestimationResult)=print(io,"PostestimationResult(status=",r.status,", flags=",r.flags,")")
function _post_finish(rc,handle)
 try
  rc==0 || throw(CoreError(rc,unsafe_string(ccall((:sb_last_error,libspecbunch),Cstring,()))))
  status,raw=_extract(handle);PostestimationResult(status,raw)
 finally
  handle==C_NULL || ccall((:sb_result_free,libspecbunch),Cvoid,(Ptr{Cvoid},),handle)
 end
end
function _post_mass(fit,kink)
 fit isa Result && fit.kind==:mass && !fit.refused || throw(ArgumentError("a coherent mass result is required"))
 k=kink===nothing ? fit.kink : Float64(kink)
 k==fit.kink || throw(ArgumentError("kink must match the mass fit"))
 (fit.B,fit.sigma,fit.raw["friction"]=="voigt" ? fit.raw["a"] : 0.,k)
end
"""
    nonbuncher_mean(z,x,mass_fit,buncher_mean_fit=nothing; mu_b=nothing,se_mu_b=0,...)

Triangular local-linear nonbuncher mean at the kink (bandwidth=.05 log units).
Degree=2 is a local quadratic check; bin_width=.0025. Binning and observed-density
smoothing occur in C++. Supply a buncher_mean result or mu_b. Returns m_NB, se,
pi_kink, flags and local_fit. SE is not a confidence interval: mass/friction and
smoothed density are fixed, with buncher-mean variance propagated assuming
independence. Report the Gaussian/Voigt pair. No plot is drawn.
"""
function nonbuncher_mean(z,x,mass_fit,buncher_mean_fit=nothing;mu_b=nothing,se_mu_b=0.,
                        kink=nothing,weights=nothing,bandwidth=.05,degree=1,bin_width=.0025,threads=nothing)
 B,s,a,k=_post_mass(mass_fit,kink)
 if buncher_mean_fit!==nothing
  mu_b===nothing || throw(ArgumentError("supply buncher_mean_fit or mu_b, not both"))
  buncher_mean_fit isa Result && buncher_mean_fit.kind==:bmean && !buncher_mean_fit.refused || throw(ArgumentError("a coherent buncher_mean result is required"))
  buncher_mean_fit.kink==k || throw(ArgumentError("buncher mean and mass fit must have the same kink"))
  mu_b=buncher_mean_fit.mu;se_mu_b=buncher_mean_fit.se
 end
 mu_b===nothing && throw(ArgumentError("buncher_mean_fit or mu_b is required"))
 (z===nothing || x===nothing) && throw(ArgumentError("z and x are required"))
 z,x,w=_inputs(z,x,weights)
 o=NonbuncherOpts(UInt32(sizeof(NonbuncherOpts)),1,_threads(threads),Int32(degree),B,s,a,k,Float64(mu_b),Float64(se_mu_b),Float64(bandwidth),Float64(bin_width));out=Ref{Ptr{Cvoid}}(C_NULL)
 GC.@preserve z x w begin
  d=MeanData(length(z),pointer(z),pointer(x),_ptr(w))
  rc=ccall((:sb_nonbuncher_mean,libspecbunch),Cint,(Ref{MeanData},Ref{NonbuncherOpts},Ref{Ptr{Cvoid}}),d,o,out)
  _post_finish(rc,out[])
 end
end
function _eti_estimate(fit,z;tax_below=nothing,tax_above=nothing,net_below=nothing,net_above=nothing,
                       kink=nothing,weights=nothing,integration_intervals=200,threads=nothing)
 B,s,a,k=_post_mass(fit,kink)
 if net_below===nothing && net_above===nothing
  (tax_below===nothing || tax_above===nothing) && throw(ArgumentError("supply both tax rates or both net-of-tax rates"))
  net_below=1-Float64(tax_below);net_above=1-Float64(tax_above)
 elseif tax_below!==nothing || tax_above!==nothing || net_below===nothing || net_above===nothing
  throw(ArgumentError("supply exactly one complete rate pair"))
 end
 z===nothing && throw(ArgumentError("z is required"));z,w=_inputs(z,weights)
 o=EtiOpts(UInt32(sizeof(EtiOpts)),1,_threads(threads),Int32(integration_intervals),B,s,a,k,fit.raw["window_lo"],fit.raw["window_hi"],Float64(net_below),Float64(net_above));out=Ref{Ptr{Cvoid}}(C_NULL)
 GC.@preserve z w begin
  d=MassData(length(z),pointer(z),_ptr(w))
  rc=ccall((:sb_eti,libspecbunch),Cint,(Ref{MassData},Ref{EtiOpts},Ref{Ptr{Cvoid}}),d,o,out)
  _post_finish(rc,out[])
 end
end
"""Elasticity of taxable income by shift-and-complete; no forward transform completion or derivative rules."""
module ETI
import ..SpecBunch: _eti_estimate
export estimate
"""estimate(fit,z; tax_below,tax_above) or net_below/net_above returns saez and continuous e/d, shoulders and flags. Input z is in log units; shoulders are finite-resolution CF functionals. No SE or confidence interval is computed. See REFERENCE.md."""
estimate(args...;kwargs...)=_eti_estimate(args...;kwargs...)
end
