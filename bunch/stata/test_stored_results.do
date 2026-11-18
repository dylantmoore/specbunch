* Compare the complete e() interface before/after optimization, not just ABI fields.
* Arguments: isolated net/PLUS root and output filename, both outside the repo.
args root output
clear all
set more off
sysdir set PLUS "`root'/plus"
sysdir set PERSONAL "`root'/personal"
net install specbunch, from("`root'/net/specbunch") replace
net install specrd, from("`root'/net/specrd") replace
capture erase "`output'"
mata:
void snapshot(string scalar path, string scalar label)
{
    real scalar f,i,j,k
    string colvector names
    string scalar name
    real matrix values
    string matrix stripes
    f=fopen(path,"a"); fput(f,"CASE "+label)
    names=sort(st_dir("e()","numscalar","*"),1)
    for(i=1;i<=rows(names);i++) fput(f,"S "+names[i]+" "+strofreal(st_numscalar("e("+names[i]+")"),"%21x"))
    names=sort(st_dir("e()","macro","*"),1)
    for(i=1;i<=rows(names);i++) fput(f,"T "+names[i]+" "+st_global("e("+names[i]+")"))
    names=sort(st_dir("e()","matrix","*"),1)
    for(i=1;i<=rows(names);i++) {
        name="e("+names[i]+")"; values=st_matrix(name)
        fput(f,"M "+names[i]+" "+strofreal(rows(values))+" "+strofreal(cols(values)))
        stripes=st_matrixrowstripe(name)
        for(j=1;j<=rows(stripes);j++) fput(f,"R "+stripes[j,1]+" "+stripes[j,2])
        stripes=st_matrixcolstripe(name)
        for(j=1;j<=rows(stripes);j++) fput(f,"C "+stripes[j,1]+" "+stripes[j,2])
        for(j=1;j<=rows(values);j++) for(k=1;k<=cols(values);k++) fput(f,strofreal(values[j,k],"%21x"))
    }
    fclose(f)
}
end
set obs 1200
set seed 520
quietly generate double z=cond(runiform()<.2,0,rnormal())+.04*rnormal()
quietly generate double x=2+.2*z+.1*rnormal()
quietly generate double d=z>0
quietly generate double wt=1+mod(_n,3)
foreach entry in gaussian voigt joint {
    local method cosine
    local friction `entry'
    if "`entry'"=="joint" {
        local method joint
        local friction gaussian
    }
    quietly specbunch mass z, method(`method') friction(`friction') threads(1) seed(1)
    quietly specbunch mass z, method(`method') friction(`friction') threads(1) seed(1)
    mata: snapshot(st_local("output"),st_local("entry"))
}
quietly specbunch mass z, threads(1) seed(1)
quietly specbunch mass z, threads(1) seed(1)
mata: snapshot(st_local("output"),"no_candidates")
foreach entry in gaussian voigt {
    local method cosine
    local friction `entry'
    quietly specbunch mass z [pw=wt], method(`method') friction(`friction') threads(1) seed(1)
    quietly specbunch mass z [pw=wt], method(`method') friction(`friction') threads(1) seed(1)
    mata: snapshot(st_local("output"),"weighted_"+st_local("entry"))
}
quietly specbunch buncher z x, threads(1)
quietly specbunch buncher z x, threads(1)
mata: snapshot(st_local("output"),"bmean")
quietly specbunch buncher z x [pw=wt], threads(1)
quietly specbunch buncher z x [pw=wt], threads(1)
mata: snapshot(st_local("output"),"weighted_bmean")
quietly specrd x d z, threads(1)
quietly specrd x d z, threads(1)
mata: snapshot(st_local("output"),"rd")
quietly specrd x d z [fw=wt], threads(1)
quietly specrd x d z [fw=wt], threads(1)
mata: snapshot(st_local("output"),"weighted_rd")
quietly replace z=0
quietly specbunch mass z, threads(1) policy(strict)
quietly specbunch mass z, threads(1) policy(strict)
mata: snapshot(st_local("output"),"refusal")
quietly specbunch buncher z x, threads(1)
quietly specbunch buncher z x, threads(1)
mata: snapshot(st_local("output"),"bmean_refusal")
display "STORED_RESULTS_PASS"
exit, clear
