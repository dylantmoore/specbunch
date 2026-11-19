module SpecBunch
using SpecBunch_jll: libspecbunch
export mass, buncher_mean, npmle, NpmleResult, POLICIES, Result, CoreError, version, default_threads
include("transport.jl")

"""Shared core default: OMP_NUM_THREADS, else available performance/physical cores."""
default_threads() = ccall((:sb_default_threads, libspecbunch), Cint, ())

const POLICIES = (:tiered,:first)
struct MassData
    n::Csize_t; z::Ptr{Cdouble}; w::Ptr{Cdouble}
end
struct MeanData
    n::Csize_t; z::Ptr{Cdouble}; x::Ptr{Cdouble}; w::Ptr{Cdouble}
end
struct MassOpts
    struct_size::UInt32; abi_version::UInt32
    threads::Int32; check_all_windows::Int32; seed::UInt64
    distribution_checks::Int32; policy::Int32; lp_iterations::Int64; lp_retry_iterations::Int64
end
struct MeanOpts
    struct_size::UInt32; abi_version::UInt32
    threads::Int32; reserved_band::Int32; reserved_visibility::Int32; R::Int32; seed::UInt64
end
function _massopts(;threads=nothing,seed=1,policy=:tiered,check_all_windows=false,lp_iterations=5000,lp_retry_iterations=31250)
    p=Symbol(policy);p in POLICIES || throw(ArgumentError("policy must be tiered or first"))
    check_all_windows isa Bool || throw(ArgumentError("check_all_windows must be boolean"))
    min(lp_iterations,lp_retry_iterations)>0 || throw(ArgumentError("LP budgets must be positive"))
    MassOpts(UInt32(sizeof(MassOpts)),1,_threads(threads),check_all_windows ? 2 : 0,UInt64(seed),1,p==:first ? 2 : 0,Int64(lp_iterations),Int64(lp_retry_iterations))
end
function _masscall(method,friction,data,o,out)
    ccall((:sb_mass,libspecbunch),Cint,(Ref{MassData},Ref{MassOpts},Cint,Cint,Ref{Ptr{Cvoid}}),data,o,method==:joint,friction==:voigt,out)
end
"""Weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits sine; the companions are log-linear OLS/LAD.
The transform NamedTuple contains t, t_scaled, cosine, cosine_se, cosine_fit; joint fits
add sine_scaled, sine_scaled_se and sine_fit. n_points describes the selected window.

`plot=true` (default, same as `:scaled`) draws the empirical cosine transform with a ±2 SE
band, the fitted curve over the whole grid and the shaded selected window into
`result.figure` when Plots.jl is loaded (joint fits add a t·S(t) panel); `:raw` uses raw t
on the x-axis instead of t·σ_pilot. Without Plots the request is a no-op with a one-time
`@info`; `plot=false` or `nothing` skips it. `plot(result; frequency=:raw)` redraws. With
`npmle=true` the NPMLE figure is `result.npmle.figure`."""
function mass(z;kink=0.0,weights=nothing,method=:cosine,friction=:gaussian,policy=:tiered,threads=nothing,seed=1,
              check_all_windows=false,lp_iterations=5000,lp_retry_iterations=31250,npmle=false,plot=true)
    z===nothing && throw(ArgumentError("z is required"))
    Symbol(method) in (:cosine,:joint) || throw(ArgumentError("method must be cosine or joint"))
    Symbol(friction) in (:gaussian,:voigt) || throw(ArgumentError("friction must be gaussian or voigt"))
    xscale=_transform_xscale(plot);npmle isa Bool || throw(ArgumentError("npmle must be boolean"))
    z,w=_inputs(z,weights);original_z=z;kink==0 || (z=z .- Float64(kink))
    o=_massopts(;threads,seed,policy,check_all_windows,lp_iterations,lp_retry_iterations)
    out=Ref{Ptr{Cvoid}}(C_NULL)
    GC.@preserve z w begin
        d=MassData(length(z),pointer(z),_ptr(w));rc=_masscall(Symbol(method),Symbol(friction),d,o,out)
        ans=_finish(rc,out[],:mass);ans.metadata[:kink]=Float64(kink);ans.metadata[:figure]=nothing
        xscale===nothing || _attach_figure!(ans,xscale,"mass")
        if npmle ans.metadata[:npmle]=ans.refused ? nothing : SpecBunch.npmle(ans,original_z;weights=w,threads) end
        ans
    end
end
"""Cosine-ratio buncher mean; coherent bands failing the fit test are flagged poor_fit."""
function buncher_mean(z,x;kink=0.0,weights=nothing,threads=nothing,seed=1,draws=499)
    (z===nothing || x===nothing) && throw(ArgumentError("z and x are required"))
    z,x,w=_inputs(z,x,weights);kink==0 || (z=z .- Float64(kink))
    o=MeanOpts(UInt32(sizeof(MeanOpts)),1,_threads(threads),2,1,Int32(draws),UInt64(seed));out=Ref{Ptr{Cvoid}}(C_NULL)
    GC.@preserve z x w begin
        d=MeanData(length(z),pointer(z),pointer(x),_ptr(w));rc=ccall((:sb_bmean,libspecbunch),Cint,(Ref{MeanData},Ref{MeanOpts},Ref{Ptr{Cvoid}}),d,o,out)
        ans=_finish(rc,out[],:bmean);ans.metadata[:kink]=Float64(kink);ans
    end
end
include("npmle.jl")
include("postestimation.jl")
export nonbuncher_mean, PostestimationResult, ETI
end
