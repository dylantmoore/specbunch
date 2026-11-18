# Local artifact QA, isolated from the user's Julia depot and checkout.
using Pkg, UUIDs, Libdl, TOML
using Pkg.Artifacts
B = get(ENV,"SB_JULIA_BUILD",joinpath(get(ENV,"SPECPKG_BUILD",joinpath(tempdir(),"spectral-build")),"julia"))
R = normpath(joinpath(@__DIR__,"../../.."))
specs = Pkg.PackageSpec[]
names = String[]
overrides = Dict{String,Any}()
for (family,name,uuid) in [
    ("bunch","SpecBunch","7035fe15-a909-5e98-9b00-51961afcec2a"),
    ("rd","SpecRD","32c5df02-eca5-59f7-98f5-23b3f2277a3a")]
    family != "bunch" && continue
    dest=joinpath(B,"packages",name)
    isdir(dest) && rm(dest;recursive=true)
    mkpath(dirname(dest)); cp(joinpath(R,family,"julia",name),dest)
    libname=(Sys.iswindows() ? "" : "lib")*lowercase(name)*"."*Libdl.dlext
    lib=joinpath(B,"build",family,libname)
    isfile(lib) || error("build core first: $lib")
    hash=create_artifact() do dir
        mkpath(joinpath(dir,"lib")); cp(lib,joinpath(dir,"lib",libname))
    end
    overrides[string(hash)]=artifact_path(hash)
    stub=joinpath(B,"packages",name*"_jll"); mkpath(joinpath(stub,"src"))
    write(joinpath(stub,"Project.toml"),"name = \"$(name)_jll\"\nuuid = \"$uuid\"\nversion = \"0.1.0+0\"\n[deps]\nArtifacts = \"56f22d72-fd6d-98f1-02f0-08ddc0907c33\"\n")
    bind_artifact!(joinpath(stub,"Artifacts.toml"),"core",hash;force=true)
    write(joinpath(stub,"src",name*"_jll.jl"),"module $(name)_jll\nusing Artifacts\nexport lib$(lowercase(name))\nconst lib$(lowercase(name)) = joinpath(artifact\"core\",\"lib\",$(repr(libname)))\nend\n")
    push!(specs,PackageSpec(path=stub),PackageSpec(path=dest));push!(names,name)
end
mkpath(joinpath(DEPOT_PATH[1],"artifacts"))
open(joinpath(DEPOT_PATH[1],"artifacts","Overrides.toml"),"w") do io
    TOML.print(io,overrides)
end
Pkg.activate(joinpath(B,"env"));Pkg.develop(specs);Pkg.test(names)
