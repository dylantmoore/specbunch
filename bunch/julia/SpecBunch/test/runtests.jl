using Test, SpecBunch
@testset "SpecBunch ABI and inputs" begin
    @test SpecBunch.version().abi == 1
    @test SpecBunch.version().package == "0.1.0"
    @test sizeof(SpecBunch.MassOpts) == 48
    @test sizeof(SpecBunch.MeanOpts) == 32
    @test_throws ArgumentError mass([1.0]; method=:bad)
    @test_throws ArgumentError mass([1.0]; policy=:bad)
    @test_throws DimensionMismatch mass(zeros(100); weights=ones(99))
    @test_throws ArgumentError mass(zeros(10,10))
    @test_throws InexactError mass(zeros(100); seed=-1)
    @test_throws InexactError mass(zeros(100); threads=1.5)
    @test mass(zeros(10)).reason == "insufficient_data"
    @test_throws CoreError mass(zeros(100); method=:joint,weights=ones(100))
    z=[sin(i*1.71)*0.4 + sin(i*0.31)*0.03 for i in 1:300]
    x=[1.0+cos(i*0.7) for i in 1:300]
    a=buncher_mean(z,x;threads=1,draws=19)
    b=buncher_mean(z,x;threads=2,draws=19)
    @test a.threads_used == 1 && b.threads_used == 2
    @test isequal(filter(p->p.first!="threads_used",a.raw),filter(p->p.first!="threads_used",b.raw))
    @test a.mu === a.raw["mu"] || isequal(a.mu,a.raw["mu"])
    @test a.refused == (a.status != 0)
    @test occursin("not a confidence interval",sprint(show,MIME"text/plain"(),a))
    @test_throws DimensionMismatch buncher_mean(z,x[1:end-1])
    missingz=Union{Missing,Float64}[z;missing]
    @test isequal(buncher_mean(missingz,[x;NaN];threads=1,draws=19).mu,a.mu)
    for (method,friction) in ((:cosine,:gaussian),(:cosine,:voigt),(:joint,:gaussian),(:joint,:voigt))
        r=mass(zeros(200);method,friction,threads=1)
        @test r.windows isa NamedTuple
        @test r.transform isa NamedTuple
        @test haskey(r.transform,:cosine)
        @test r.se_scale == "log(B)"
        @test isequal(r.band,r.window)
        @test r.flags isa AbstractVector
        @test r.flag_bits == Int(get(r.raw,"flags",0))
        GC.gc()
        @test r.raw isa Dict{String,Any}
    end
    withenv("OMP_NUM_THREADS"=>nothing) do
        @test buncher_mean(z,x;).threads_used == SpecBunch.default_threads()
        for t in (1,3,16)
            @test buncher_mean(z,x;threads=t).threads_used == t
        end
    end
    withenv("OMP_NUM_THREADS"=>"3") do
        @test SpecBunch.default_threads() == 3
        @test buncher_mean(z,x;).threads_used == 3
        @test buncher_mean(z,x;threads=16).threads_used == 16
    end
end

@testset "Spline postestimation" begin
    uniform(i,stream)=Float64(((UInt64(i)+UInt64(stream)*0x9e3779b97f4a7c15)*0xbf58476d1ce4e5b9) >> 11)*2.0^-53
    normal(i,stream)=sum(uniform(i,stream+j) for j in 0:11)-6.0
    # The package fixture generator uses independent mixed integer streams.
    function u(i,s)
        x=UInt64(i)+(UInt64(23)+UInt64(s)*0x9e3779b97f4a7c15)
        x=(x ⊻ (x>>30))*0xbf58476d1ce4e5b9;x=(x ⊻ (x>>27))*0x94d049bb133111eb
        Float64((x ⊻ (x>>31))>>11)*2.0^-53
    end
    norm(i,s)=sum(u(i,s+j) for j in 0:11)-6.0
    z=[(u(i,40)<.2 ? 0. : norm(i,1))+.04*norm(i,20) for i in 0:1499]
    fit=mass(z;threads=1,npmle=true)
    @test !fit.refused && fit.npmle isa NpmleResult
    result=npmle(fit,z;threads=1,grid_points=101,lambda_points=5,refinement_points=3)
    if result.status != "converged"
        @info "NPMLE convergence diagnostics" B=fit.B sigma=fit.sigma diagnostics=Dict(k=>result.raw[k] for k in ("lambda","kkt","iterations","lambda_path.lambda","lambda_path.kkt"))
    end
    @test result.status=="converged"
    @test result.curves isa NamedTuple && result.histogram isa NamedTuple
    @test all(result.coefficients .>=0)
    @test abs(result.empirical_cdf[end]-1)<1e-12
    @test all(diff(result.latent_cdf).>=-1e-12)
    @test result.histogram.fitted_share == result.histogram.nonbuncher_share+result.histogram.buncher_share
    # Without Plots.jl loaded the default plot request is a no-op; plot=false is always silent.
    if !isdefined(Main,:Plots)
        @test result.figure===nothing && fit.npmle.figure===nothing && fit.figure===nothing
    end
    @test mass(z;threads=1,plot=false).figure===nothing
    @test_throws ArgumentError mass(z;plot=:log)
    @test npmle(fit,z;threads=1,grid_points=51,lambda_points=3,refinement_points=3,plot=false).figure===nothing
    @test_throws ArgumentError npmle(fit,z;plot=:hist)
end

@testset "Local nonbuncher mean and ETI" begin
    z=collect(range(-.5,.5,length=5000));x=fill(.4,length(z))
    f=SpecBunch.Result(0,Dict{String,Any}("B"=>.08,"sigma"=>.03,"a"=>0.,"friction"=>"gaussian","window_lo"=>12.,"window_hi"=>45.),:mass)
    f.metadata[:kink]=0.
    r=nonbuncher_mean(z,x,f;mu_b=.4)
    @test !r.refused
    @test abs(r.m_NB-.4)<1e-10
    @test length(r.local_fit.center)==40
    @test r.pi_kink>0
    @test abs(nonbuncher_mean(z,x,f;mu_b=.4,degree=2,threads=4).m_NB-.4)<1e-10
    e=SpecBunch.ETI.estimate(f,z;net_below=1.34,net_above=1.)
    @test !e.refused
    @test e.saez.e≈e.saez.d/log(1.34)
    @test_throws SpecBunch.CoreError SpecBunch.ETI.estimate(f,z;net_below=1.,net_above=1.)
end

@testset "NPMLE ETI shoulders" begin
 z=Float64.((0:399).-200)./256
 fixture=Dict(split(line,',')[1]=>parse.(Float64,split(line,',')[2:end]) for line in readlines(joinpath(@__DIR__,"eti_shoulders.csv"))[2:end])
 f=SpecBunch.Result(0,Dict{String,Any}("B"=>.08,"sigma"=>.03,"a"=>0.,"kink"=>0.,"friction"=>"gaussian","window_lo"=>12.,"window_hi"=>45.),:mass)
 for basis in (:linear,:cubic)
  n=npmle(f,z;basis,plot=false,core_intervals=4,grid_points=25,lambda_points=3,refinement_points=3,threads=1,net_below=1.34,net_above=1.)
  @test n.latent_left≈fixture[String(basis)][1] rtol=1e-8 atol=1e-12
  @test n.eti.npmle.saez.e≈fixture[String(basis)][3] rtol=1e-8 atol=1e-12
  @test n.eti.npmle.basis==String(basis)
  @test n.eti.npmle.shoulder_left==n.latent_left==n.coefficients[1]
  @test n.eti.npmle.saez.e≈2*f.B/(n.latent_left+n.latent_right)/log(1.34)
  f.metadata[:npmle]=n
  @test SpecBunch.ETI.estimate(f,z;net_below=1.34,net_above=1.).npmle.saez.e==n.eti.npmle.saez.e
 end
end
