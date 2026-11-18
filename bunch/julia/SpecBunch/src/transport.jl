
"""A C ABI error, retaining the integer code and calling thread's message."""
struct CoreError <: Exception
    code::Int
    message::String
end
Base.showerror(io::IO, e::CoreError) = print(io, "CoreError(", e.code, "): ", e.message)

"""Owned output: every C field in raw, plus summary properties and candidate columns."""
struct Result
    status::Int
    raw::Dict{String,Any}
    kind::Symbol
    metadata::Dict{Symbol,Any}
end
Result(status,raw,kind)=Result(status,raw,kind,Dict{Symbol,Any}())
function Base.getproperty(r::Result, s::Symbol)
    s in (:status, :raw, :kind, :metadata) && return getfield(r, s)
    haskey(getfield(r,:metadata),s) && return getfield(r,:metadata)[s]
    raw = getfield(r, :raw); kind = getfield(r, :kind)
    s == :refused && return r.status != 0
    s == :reason && return get(raw, "reason", "")
    s == :flags && return split(get(raw, "flags_str", ""),','; keepempty=false)
    s == :flag_bits && return Int(get(raw, "flags", 0))
    s == :tier && return haskey(raw,"tier") ? Int(raw["tier"]) : nothing
    s == :policy && return get(raw,"policy","")
    s == :se_scale && return kind == :mass ? "log(B)" : "levels"
    s == :se && return get(raw, kind == :mass ? "se_log_B" : "se", NaN)
    s == :window && return get(raw, "win", r.band)
    s == :band && return kind == :mass ? (get(raw,"window_lo",NaN),get(raw,"window_hi",NaN)) : (get(raw,"band_lo",NaN),get(raw,"band_hi",NaN))
    if s == :windows
        keys_=sort!([k for k in keys(raw) if startswith(k,"windows.")])
        checks=("cdf_inconsistent","moments_inconsistent","no_valid_distribution","distribution_check_inconclusive","residual_cdf_inconsistent","poor_fit")
        columns=Dict{Symbol,Any}(Symbol(k[9:end]) => (kind==:mass && k[9:end] in checks ? [x<0 ? "not_evaluated" : x==0 ? "fail" : "pass" for x in raw[k]] : raw[k]) for k in keys_)
        if kind==:bmean && haskey(columns,:poor_fit) columns[:flags]=[x>0 ? "poor_fit" : "" for x in columns[:poor_fit]] end
        return (;(k=>columns[k] for k in sort!(collect(keys(columns))))...)
    end
    if s in (:simple_ols,:simple_lad)
        prefix=String(s)*"."
        return Dict(k[length(prefix)+1:end]=>v for (k,v) in raw if startswith(k,prefix))
    end
    s == :diagnostics && return raw
    haskey(raw,String(s)) && return raw[String(s)]
    throw(ArgumentError("unknown result property $s; all named fields are in raw"))
end
Base.propertynames(r::Result, private::Bool=false) = Tuple(unique(vcat(
    [:status,:raw,:kind,:refused,:reason,:flags,:tier,:se,:se_scale,:window,:band,:windows,:diagnostics,:simple,:flag_bits,:policy],
    r.kind == :mass ? [:simple_ols,:simple_lad] : r.kind == :bmean ? Symbol[] : Symbol[],
    Symbol.(collect(keys(r.raw))))))
function Base.show(io::IO, r::Result)
    key = r.kind == :mass ? :B : r.kind == :bmean ? :mu : :tau
    print(io, "Result(", r.kind, ", ", key, "=", get(r.raw,String(key),NaN),
          ", SE[", r.se_scale, "]=", r.se, ", status=", r.status,
          ", tier=",r.tier,", flags=",r.flags, isempty(r.reason) ? "" : ", reason="*r.reason, ")")
end
function Base.show(io::IO, ::MIME"text/plain", r::Result)
    println(io, r.kind, r.refused ? " | refused" : " | estimate", " (status ",r.status,")")
    println(io, rpad("field",16), "value")
    fields = r.kind == :mass ? (:B,:sigma,:se) : r.kind == :bmean ? (:mu,:se,:pointwise_median) : (:tau,:se)
    for field in fields
        println(io,rpad(string(field),16),getproperty(r,field))
    end
    println(io,rpad("band/window",16),r.window)
    println(io,rpad("tier / flags",16),r.tier," / ",join(r.flags,", "))
    isempty(r.reason) || println(io,rpad("reason",16),r.reason)
    print(io,"Model-based SE [",r.se_scale,"]; not a confidence interval.")
    c=r.windows
    if !isempty(c)
        print(io,"\nCandidate columns (",length(first(c))," rows): ",join(string.(keys(c)),", "))
    end
end

_vector(x::AbstractVector) = x isa Vector{Float64} ? x : Float64[ismissing(v) ? NaN : Float64(v) for v in x]
_vector(x) = throw(ArgumentError("inputs must be one-dimensional vectors"))
_ptr(x) = x === nothing ? Ptr{Cdouble}(C_NULL) : pointer(x)
function _inputs(xs...)
    a = map(x->x === nothing ? nothing : _vector(x), xs)
    lengths = [length(x) for x in a if x !== nothing]
    isempty(lengths) || all(==(first(lengths)),lengths) || throw(DimensionMismatch("column lengths differ"))
    a
end
_threads(x) = x === nothing ? Int32(0) : Int32(x)

function _extract(h)
    raw = Dict{String,Any}()
    status = Int(ccall((:sb_result_status, libspecbunch), Cint, (Ptr{Cvoid},), h))
    for kind in 0:2
        p = ccall((:sb_result_names, libspecbunch), Cstring, (Ptr{Cvoid},Cint),h,kind)
        p == C_NULL && continue
        for name in split(unsafe_string(p),'\n'; keepempty=false)
            if kind == 0
                v=Ref{Cdouble}()
                rc=ccall((:sb_result_scalar, libspecbunch),Cint,(Ptr{Cvoid},Cstring,Ref{Cdouble}),h,name,v)
                rc == 0 || error("scalar getter failed: $name ($rc)")
                raw[name]=v[]
            elseif kind == 1
                n=Ref{Csize_t}(0)
                v=ccall((:sb_result_vector, libspecbunch),Ptr{Cdouble},(Ptr{Cvoid},Cstring,Ref{Csize_t}),h,name,n)
                raw[name]=n[] == 0 ? Float64[] : copy(unsafe_wrap(Vector{Float64},v,Int(n[])))
            else
                v=ccall((:sb_result_string, libspecbunch),Cstring,(Ptr{Cvoid},Cstring),h,name)
                raw[name]=v == C_NULL ? "" : unsafe_string(v)
            end
        end
    end
    status, raw
end
function _finish(rc, h, kind)
    try
        rc == 0 || throw(CoreError(rc,unsafe_string(ccall((:sb_last_error, libspecbunch),Cstring,()))))
        status, raw = _extract(h)
        Result(status,raw,kind)
    finally
        h == C_NULL || ccall((:sb_result_free, libspecbunch),Cvoid,(Ptr{Cvoid},),h)
    end
end
"""Return package/core versions and the loaded core's ABI version."""
version() = (package="0.1.0", core=unsafe_string(ccall((:sb_version,libspecbunch),Cstring,())),
             abi=ccall((:sb_abi_version,libspecbunch),UInt32,()))
function __init__()
    version().abi == 1 || error("unsupported core ABI")
end
