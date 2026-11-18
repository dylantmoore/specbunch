using Libdl
B=get(ENV,"SB_JULIA_BUILD",joinpath(get(ENV,"SPECPKG_BUILD",joinpath(tempdir(),"spectral-build")),"julia"))
h=Ref{Ptr{Cvoid}}(C_NULL)
z=[sin(i*0.17)+cos(i*0.31)*0.04 for i in 1:300]
x=cos.(z)
function foreign()
    h[] == C_NULL && (h[]=Libdl.dlopen(joinpath(B,"foreign_omp.dylib"),Libdl.RTLD_LOCAL))
    ccall(Libdl.dlsym(h[],:foreign_parallel),Cint,(Cint,),3) == 3 || error("foreign team mismatch")
end
function bunch(t=1)
    Base.invokelatest(SpecBunch.mass,z;threads=t,distribution_checks=0)
end
function rdcall(t=1)
    Base.invokelatest(SpecRD.rd,z,Float64.(z .> 0),x;threads=t)
end
haskey(ENV,"KMP_DUPLICATE_LIB_OK") && error("duplicate-runtime override must be absent")
for node in split(only(ARGS),",")
    node == "bunch" && (@eval using SpecBunch)
    node == "rd" && (@eval using SpecRD)
    node == "foreign" ? foreign() : node == "bunch" ? Base.invokelatest(bunch) : Base.invokelatest(rdcall)
end
for t in (1,2,4,8)
    b=bunch(t); r=rdcall(t); foreign()
    jobs=[Threads.@spawn(bunch(t)),Threads.@spawn(rdcall(t))]
    bb,rr=fetch.(jobs)
    isequal(bb.raw,b.raw) || error("concurrent bunch mismatch")
    isequal(rr.raw,r.raw) || error("concurrent RD mismatch")
    foreign()
end
haskey(ENV,"KMP_DUPLICATE_LIB_OK") && error("package changed override")
println("PASS ",only(ARGS)," interleaved/concurrent 1/2/4/8; foreign team=3")
