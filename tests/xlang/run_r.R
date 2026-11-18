#!/usr/bin/env Rscript
# All inputs are simulated. Reads the unchanged spectral-xlang-v1 oracle.
args <- commandArgs(trailingOnly = TRUE)
root <- path.expand(if(length(args)) args[[1L]] else Sys.getenv("SPECPKG_BUILD", file.path(tempdir(),"spectral-build")))
B <- file.path(root, "r")
# Never leave a stale PASS receipt behind a failed run.
unlink(file.path(B, "xlang_result.json"))
options(error = function() {
  jsonlite::write_json(list(result = "FAIL", error = geterrmessage()),
                      file.path(B, "xlang_result.json"), auto_unbox = TRUE, pretty = TRUE)
  quit(save = "no", status = 1L)
})
.libPaths(c(file.path(B, "lib"), .libPaths()))
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[[1L]])
xlang <- dirname(normalizePath(gsub("~+~", " ", script, fixed=TRUE)))
fixture_root <- file.path(root, "xlang/fixtures")
manifest <- jsonlite::fromJSON(gzfile(file.path(fixture_root, "manifest.json.gz")), simplifyVector=FALSE)
stopifnot(manifest$format == "spectral-xlang-v1")
# Only the normative generator runs in Python, never an estimator.
load_recipe <- function(recipe) {
  tmp <- tempfile("recipe-", tmpdir=file.path(B,"tmp")); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE))
  jsonlite::write_json(recipe, file.path(tmp,"recipe.json"), auto_unbox=TRUE)
  code <- paste0("import sys,json;sys.path.insert(0,sys.argv[1]);from fixtures import simulate;from pathlib import Path;p=Path(sys.argv[2]);a=simulate(**json.loads((p/'recipe.json').read_text()));[(v.astype('<f8').tofile(p/(k+'.bin'))) for k,v in a.items()]")
  python <- Sys.getenv("SPECTRAL_PYTHON", file.path(root,"py/xlang-venv/bin/python"))
  status <- system2(python,c("-B", "-c", shQuote(code),shQuote(xlang),shQuote(tmp)))
  stopifnot(status==0L)
  paths <- list.files(tmp,"\\.bin$",full.names=TRUE)
  setNames(lapply(paths,function(p) readBin(p,"double",n=recipe$n,size=8,endian="little")),sub("\\.bin$","",basename(paths)))
}
load_case <- function(case) {
  if(!is.null(case$recipe)) return(load_recipe(case$recipe))
  lapply(case$columns,function(col) {
    p <- file.path(fixture_root,col$file)
    stopifnot(digest::digest(file=p,algo="sha256")==col$sha256, file.info(p)$size==8*col$n)
    readBin(p,"double",n=col$n,size=8,endian="little")
  })
}
public <- function(case,a,threads) {
  opts <- c(case$options,list(threads=threads))
  if(case$entry=="npmle") {
    fit <- do.call(specbunch::mass,c(list(z=a$z,weights=a$weights,threads=threads),case$mass_options))
    return(do.call(specbunch::npmle,c(list(fit=fit,z=a$z,weights=a$weights),opts)))
  }
  if(case$entry=="rd") return(do.call(specrd::rd,c(list(w=a$w,d=a$d,y=a$y,weights=a$weights),opts)))
  if(case$entry=="bmean") return(do.call(specbunch::buncher_mean,c(list(z=a$z,x=a$x,weights=a$weights),opts)))
  do.call(specbunch::mass,c(list(z=a$z,weights=a$weights,method=sub("mass_","",case$entry)),opts))
}
bits <- function(x) {
  b <- writeBin(as.double(x),raw(),size=8,endian="big")
  if(!length(b)) return(character())
  apply(matrix(sprintf("%02x",as.integer(b)),nrow=8),2,paste0,collapse="")
}
compare <- function(ans,e,id) {
  if((if(!is.null(ans$status_code)) ans$status_code else ans$status)!=e$status) stop(id,": status")
  for(k in 0:2) {
    expected <- e[[c("scalars","vectors","strings")[[k+1L]]]]
    actual <- ans$raw[ans$kinds==k]
    if(!setequal(names(actual),names(expected))) stop(id,": field names, kind ",k)
    for(name in names(expected)) {
      want <- expected[[name]]
      good <- if(k==2L) identical(actual[[name]],want) else {
        want <- if(k==0L) want$bits else vapply(want,`[[`,"", "bits")
        identical(unname(bits(actual[[name]])),unname(want))
      }
      if(!good) stop(id,": bit mismatch: ",name)
    }
  }
  1+length(e$scalars)+length(e$strings)+sum(lengths(e$vectors))
}
# Both libraries are initialized in this process, with no KMP override.
stopifnot(!nzchar(Sys.getenv("KMP_DUPLICATE_LIB_OK")))
print(specbunch::version())
values <- 0L; calls <- 0L
cases <- manifest$cases
if(nzchar(Sys.getenv("SPECTRAL_CASE_FILTER"))) cases <- Filter(function(x) grepl(Sys.getenv("SPECTRAL_CASE_FILTER"),x$id),cases)
for(case in cases) {
  a <- load_case(case)
  first <- NULL
  for(threads in c(1L,2L,4L,8L,0L)) {
    ans <- public(case,a,threads)
    expected_threads <- if(threads > 0) threads else if(case$entry=="rd") specrd::default_threads() else specbunch::default_threads()
    build <- if(case$entry=="rd") specrd::version()$build else specbunch::version()$build
    if(grepl("openmp=off",build,fixed=TRUE)) expected_threads <- 1L
    stopifnot(ans$raw$threads_used == expected_threads)
    keep <- names(ans$raw) != "threads_used"
    ans$raw <- ans$raw[keep]; ans$kinds <- ans$kinds[keep]
    values <- values + compare(ans,case$expected,paste(case$id,threads))
    if(is.null(first)) first <- ans$raw else stopifnot(identical(first,ans$raw))
    calls <- calls + 1L
  }
  cat(case$id,"bitwise 1/2/4/8/default PASS\n"); flush.console()
  rm(a,ans,first); invisible(gc())
}
receipt <- list(cases=length(cases),calls=calls,compared_values=values,threads=c(1,2,4,8,0),result="PASS",builds=list(specbunch::version()))
jsonlite::write_json(receipt,file.path(B,"xlang_result.json"),pretty=TRUE,auto_unbox=TRUE)
print(receipt)
