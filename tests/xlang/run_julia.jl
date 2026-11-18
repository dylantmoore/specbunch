# spectral-xlang-v1: unchanged baseline, all public APIs and all thread settings.
using SpecBunch, JSON3, SHA
const B = get(ENV,"SB_JULIA_BUILD",joinpath(get(ENV,"SPECPKG_BUILD",joinpath(tempdir(),"spectral-build")),"julia"))
const FIXTURES = joinpath(dirname(B),"xlang","fixtures")

# Exact port of fixtures.py. No reassociation/FMA; only dyadic arithmetic.
function simulate(n; seed=1, shape="base", rd=false)
    function uniform(stream)
        x=UInt64.(0:n-1) .+ (UInt64(seed)+UInt64(stream)*0x9e3779b97f4a7c15)
        x=(x .⊻ (x .>> 30)) .* 0xbf58476d1ce4e5b9
        x=(x .⊻ (x .>> 27)) .* 0x94d049bb133111eb
        x=x .⊻ (x .>> 31)
        Float64.(x .>> 11) .* 2.0^-53
    end
    function normal(stream)
        x=zeros(n)
        for j in 0:11
            x .+= uniform(stream+j)
        end
        x .- 6.0
    end
    latent=normal(1); noise=normal(20)
    if rd
        d=Float64.(uniform(40) .< ifelse.(latent .>= 0,0.8,0.2))
        return Dict("w"=>latent .+ 0.12.*noise, "d"=>d,
                    "y"=>(1.5.*d .+ 0.25.*latent) .+ normal(60).*0.2)
    end
    member=uniform(40) .< 0.2
    z=ifelse.(member,0.0,latent) .+ 0.04.*noise
    if shape == "skew"
        z=ifelse.(member,0.0,(latent.*latent).*0.4 .- 0.35) .+ 0.08.*noise
    elseif shape == "refusal"
        z=latent
    end
    Dict("z"=>z,"x"=>ifelse.(member,2.0,0.5 .+ latent.*0.2) .+ normal(60).*0.1)
end
function arrays(case)
    if haskey(case,:recipe)
        r=case.recipe
        return simulate(Int(r.n); seed=Int(r.seed), shape=get(r,:shape,"base"),rd=get(r,:rd,false))
    end
    out=Dict{String,Vector{Float64}}()
    for (key,col) in pairs(case.columns)
        bytes=read(joinpath(FIXTURES,col.file))
        bytes2hex(sha256(bytes)) == col.sha256 || error("input SHA mismatch")
        u=reinterpret(UInt64,bytes)
        a=reinterpret(Float64,ltoh.(u))
        length(a) == col.n || error("input length mismatch")
        out[String(key)]=a
    end
    out
end
function public_call(case,a,threads)
    opts=Dict(Symbol(k)=>v for (k,v) in pairs(case.options))
    opts[:threads]=threads
    weights=get(a,"weights",nothing)
    if case.entry == "npmle"
        mo=Dict(Symbol(k)=>(k==:method ? Symbol(v) : v) for (k,v) in pairs(case.mass_options))
        fit=SpecBunch.mass(a["z"];weights,threads,mo...)
        return SpecBunch.npmle(fit,a["z"];weights,opts...)
    end
    case.entry == "rd" && return SpecRD.rd(a["w"],a["d"],a["y"];weights,opts...)
    case.entry == "bmean" && return SpecBunch.buncher_mean(a["z"],a["x"];weights,opts...)
    SpecBunch.mass(a["z"];weights,method=Symbol(case.entry[6:end]),opts...)
end
bits(x)=string(reinterpret(UInt64,Float64(x));base=16,pad=16)
function compare(r,e)
    n=0
    (r isa SpecBunch.NpmleResult ? r.status_code : r.status) == e.status || error("handle status mismatch")
    n+=1
    expected=Set(String(k) for kind in (:scalars,:vectors,:strings) for k in keys(e[kind]))
    Set(keys(r.raw)) == expected || error("named output key mismatch")
    for (k,v) in pairs(e.scalars)
        bits(r.raw[String(k)]) == v.bits || error("scalar $k: $(bits(r.raw[String(k)])) != $(v.bits)")
        n+=1
    end
    for (k,v) in pairs(e.vectors)
        got=r.raw[String(k)]
        length(got) == length(v) || error("vector $k length")
        for i in eachindex(got)
            bits(got[i]) == v[i].bits || error("vector $k index $i")
            n+=1
        end
    end
    for (k,v) in pairs(e.strings)
        r.raw[String(k)] == v || error("string $k")
        n+=1
    end
    n
end
function main()
    manifest_path=joinpath(FIXTURES,"manifest.json.gz")
    manifest=JSON3.read(read(`gzip -dc $manifest_path`,String))
    manifest.format == "spectral-xlang-v1" || error("unexpected format")
    comparisons=0; calls=0
    for case in manifest.cases
        a=arrays(case)
        for threads in (1,2,4,8,0)
            r=public_call(case,a,threads)
            expected_threads = threads > 0 ? threads : (case.entry == "rd" ? SpecRD.default_threads() : SpecBunch.default_threads())
            @assert r.raw["threads_used"] == expected_threads
            delete!(r.raw,"threads_used")
            try
                comparisons+=compare(r,case.expected)
            catch err
                error("$(case.id), threads=$threads: $(sprint(showerror,err))")
            end
            calls+=1
        end
        println("PASS ",case.id); flush(stdout)
        GC.gc()
    end
    receipt=(fixtures=length(manifest.cases),calls=calls,comparisons=comparisons,
             threads=[1,2,4,8,0],manifest_sha256=bytes2hex(sha256(read(manifest_path))),
             julia=string(VERSION),package=SpecBunch.version())
    write(joinpath(B,"xlang_result.json"),JSON3.write(receipt))
    println(JSON3.write(receipt))
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
