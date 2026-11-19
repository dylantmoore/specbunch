# Plot smoke test; needs Plots.jl in the active environment (not a test dependency).
# Usage: julia --project=<env with SpecBunch and Plots> test/plots.jl [output.png]
using SpecBunch, Plots, Test
gr()
u(i,s)=(x=UInt64(i)+(UInt64(23)+UInt64(s)*0x9e3779b97f4a7c15);x=(x ⊻ (x>>30))*0xbf58476d1ce4e5b9;
        x=(x ⊻ (x>>27))*0x94d049bb133111eb;Float64((x ⊻ (x>>31))>>11)*2.0^-53)
nrm(i,s)=sum(u(i,s+j) for j in 0:11)-6.0
z=[(u(i,40)<.15 ? 0. : .2+nrm(i,1))+.08*nrm(i,20) for i in 0:3999]
fit=mass(z;threads=1)
@testset "npmle plots" begin
    r=npmle(fit,z;threads=1)
    @test r.figure isa Plots.Plot && length(r.figure.subplots)==2
    @test length(npmle(fit,z;threads=1,grid_points=101,lambda_points=3,refinement_points=3,plot=:cdf).figure.subplots)==1
    @test npmle(fit,z;threads=1,grid_points=101,lambda_points=3,refinement_points=3,plot=false).figure===nothing
    @test plot(r;kind=:pdf) isa Plots.Plot
    isempty(ARGS) || (savefig(r.figure,ARGS[1]);println("saved ",ARGS[1]," B=",fit.B," sigma=",fit.sigma))
end
@testset "mass transform plots" begin
    @test fit.figure isa Plots.Plot && length(fit.figure.subplots)==1
    @test mass(z;threads=1,plot=false).figure===nothing && mass(z;threads=1,plot=nothing).figure===nothing
    j=mass(z;method=:joint,friction=:voigt,threads=1)
    @test j.figure isa Plots.Plot && length(j.figure.subplots)==2
    @test plot(fit;frequency=:raw) isa Plots.Plot
    both=mass(z;threads=1,npmle=true)
    @test both.figure isa Plots.Plot && both.npmle.figure isa Plots.Plot && both.figure!==both.npmle.figure
    @test_throws ArgumentError plot(buncher_mean(z,z;threads=1,draws=19))
    length(ARGS)>=3 && (savefig(fit.figure,ARGS[2]);savefig(plot(fit;frequency=:raw),replace(ARGS[2],".png"=>"_raw.png")))
end
