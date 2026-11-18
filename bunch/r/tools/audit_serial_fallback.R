#!/usr/bin/env Rscript
# Characterize the unsupported-SIMD fallback without changing the reference expected outputs.
# Exact suffixes distinguish rank1.B / c2_summaries numerical diagnostics from decisions.
args <- commandArgs(trailingOnly=TRUE)
root <- path.expand(if(length(args)) args[[1L]] else Sys.getenv("SPECPKG_BUILD", file.path(tempdir(),"spectral-build")))
B <- file.path(root,"r")
.libPaths(c(file.path(B,"fallback-lib"),.libPaths()))
repo <- Sys.getenv("SPECTRAL_REPO", getwd())
xlang <- file.path(repo,"tests/xlang")
fixture_root <- file.path(root,"xlang/fixtures")
lines <- readLines(file.path(xlang,"run_r.R"))
eval(parse(text=lines[grep('^load_recipe <-',lines):(grep('^# Both libraries',lines)-1L)]))
manifest <- jsonlite::fromJSON(gzfile(file.path(fixture_root,"manifest.json.gz")),simplifyVector=FALSE)
versions <- list(specbunch=specbunch::version(),specrd=specrd::version())
stopifnot(all(vapply(versions,function(x) !grepl("fopenmp-simd",x$build,fixed=TRUE),TRUE)))
decode <- function(h) {
  if(!length(h)) return(numeric())
  b <- substring(rep(h,each=8L),rep(seq(1L,15L,2L),length(h)),rep(seq(2L,16L,2L),length(h)))
  readBin(as.raw(strtoi(b,16L)),"double",n=length(h),size=8L,endian="big")
}
# A decision flip is a changed result status, any string, or a changed discrete
# status/check/flag/tier/mask/window-index/rank field (including candidate fields).
decision_field <- function(name) grepl("(^|[._])(status|flags|tier|fallback|mask|rank|range_ok|a|e|kind|answered)$",name,ignore.case=TRUE)
cache <- file.path(B,"fallback-cache"); dir.create(cache,showWarnings=FALSE)
identity <- list(manifest=digest::digest(file=file.path(fixture_root,"manifest.json.gz"),algo="sha256"),
                 binaries=vapply(c("specbunch","specrd"),function(pkg) digest::digest(file=getLoadedDLLs()[[pkg]][["path"]],algo="sha256"),""))
identity_path <- file.path(cache,"identity.rds")
if(file.exists(identity_path)) stopifnot(identical(readRDS(identity_path),identity)) else saveRDS(identity,identity_path)
rows <- list(); flips <- list(); shapes <- list(); comparisons <- 0; numeric_mismatches <- 0
unmatched_numeric_values <- 0
max_relative <- 0; max_absolute <- 0; max_zero_absolute <- 0
relative_at <- NULL; absolute_at <- NULL; nonfinite_mismatches <- 0
for(case in manifest$cases) {
  cachefile <- file.path(cache,paste0(case$id,".rds"))
  if(file.exists(cachefile)) ans <- readRDS(cachefile) else {
    a <- load_case(case); ans <- public(case,a,1L); saveRDS(ans,cachefile,compress=TRUE); rm(a)
  }
  e <- case$expected
  local_bits <- 0; local_flips <- list(); local_shapes <- list(); local_max <- 0
  comparisons <- comparisons+1
  if(ans$status!=e$status) local_flips$handle_status <- list(actual=ans$status,expected=e$status)
  for(k in 0:2) {
    expected <- e[[c("scalars","vectors","strings")[[k+1L]]]]
    actual <- ans$raw[ans$kinds==k]
    if(!setequal(names(actual),names(expected))) local_shapes[[paste0("field_names_kind_",k)]] <- list(missing=setdiff(names(expected),names(actual)),added=setdiff(names(actual),names(expected)))
    for(name in names(expected)) {
      if(k==2L) {
        comparisons <- comparisons+1
        if(!identical(actual[[name]],expected[[name]])) local_flips[[name]] <- list(actual=actual[[name]],expected=expected[[name]])
        next
      }
      want_bits <- if(k==0L) expected[[name]]$bits else vapply(expected[[name]],`[[`,"","bits")
      got <- actual[[name]]; want <- decode(want_bits)
      comparisons <- comparisons+length(want)
      if(length(got)!=length(want)) {
        local_shapes[[name]] <- list(actual_length=length(got),expected_length=length(want))
        unmatched_numeric_values <- unmatched_numeric_values+abs(length(got)-length(want))
        common <- seq_len(min(length(got),length(want)))
        got <- got[common]; want <- want[common]; want_bits <- want_bits[common]
      }
      bad_bits <- which(bits(got)!=want_bits)
      local_bits <- local_bits+length(bad_bits)
      finite <- is.finite(got)&is.finite(want)
      semantically_equal <- (is.nan(got)&is.nan(want)) | (!is.na(got)&!is.na(want)&got==want)
      changed <- which(!semantically_equal)
      nonfinite_mismatches <- nonfinite_mismatches+sum(!finite[changed])
      if(decision_field(name) && length(changed)) local_flips[[name]] <- list(count=length(changed),indices=head(changed,10L),actual=head(got[changed],10L),expected=head(want[changed],10L))
      indices <- which(finite)
      if(length(indices)) {
        delta <- abs(got[indices]-want[indices]); i <- which.max(delta)
        if(delta[i]>max_absolute) { max_absolute <- delta[i]; absolute_at <- list(case=case$id,field=name,index=indices[i],actual=got[indices[i]],expected=want[indices[i]]) }
        zero <- indices[want[indices]==0]
        if(length(zero)) max_zero_absolute <- max(max_zero_absolute,abs(got[zero]))
        nz <- indices[want[indices]!=0]
        if(length(nz)) {
          rel <- abs(got[nz]-want[nz])/abs(want[nz]); i <- which.max(rel)
          local_max <- max(local_max,rel[i])
          if(rel[i]>max_relative) { max_relative <- rel[i]; relative_at <- list(case=case$id,field=name,index=nz[i],actual=got[nz[i]],expected=want[nz[i]]) }
        }
      }
    }
  }
  numeric_mismatches <- numeric_mismatches+local_bits
  if(length(local_flips)) flips[[case$id]] <- local_flips
  if(length(local_shapes)) shapes[[case$id]] <- local_shapes
  rows[[case$id]] <- list(entry=case$entry,numeric_bit_mismatches=local_bits,max_relative=local_max,decision_field_flips=length(local_flips),shape_changes=length(local_shapes))
  cat(case$id,":",local_bits,"numeric bit mismatches;",length(local_flips),"decision/string fields changed;",length(local_shapes),"shape changes; max relative",format(local_max,digits=12),"\n");flush.console()
  rm(ans);invisible(gc())
}
receipt <- list(result="CHARACTERIZED",cases=length(rows),calls=length(rows),threads=1,compared_values=comparisons,numeric_bit_mismatches=numeric_mismatches,max_relative=max_relative,max_relative_definition="abs(actual-expected)/abs(expected), finite nonzero expected values",relative_at=relative_at,max_absolute=max_absolute,absolute_at=absolute_at,max_absolute_at_expected_zero=max_zero_absolute,nonfinite_mismatches=nonfinite_mismatches,decision_flip_cases=length(flips),decision_flips=flips,shape_change_cases=length(shapes),shape_changes=shapes,unmatched_numeric_values=unmatched_numeric_values,by_case=rows,builds=versions)
# Separately quantify the primary point estimate; ranked candidate positions may
# refer to different windows when the fallback changes candidate ordering.
headline <- list()
for(case in manifest$cases) {
  field <- if(case$entry=="rd") "tau" else if(case$entry=="bmean") "mu" else "B"
  got <- readRDS(file.path(cache,paste0(case$id,".rds")))$raw[[field]]
  want <- decode(case$expected$scalars[[field]]$bits)
  row <- headline[[case$entry]]
  if(is.null(row)) row <- list(finite_pairs=0L,max_relative=0,max_absolute=0,relative_at=NULL,nonfinite_changes=0L)
  if(is.finite(got)&&is.finite(want)) {
    row$finite_pairs <- row$finite_pairs+1L
    row$max_absolute <- max(row$max_absolute,abs(got-want))
    if(want!=0 && abs(got-want)/abs(want)>row$max_relative) {
      row$max_relative <- abs(got-want)/abs(want)
      row$relative_at <- list(case=case$id,field=field,actual=got,expected=want)
    }
  } else if(!identical(bits(got),case$expected$scalars[[field]]$bits)) row$nonfinite_changes <- row$nonfinite_changes+1L
  headline[[case$entry]] <- row
}
receipt$headline_point_estimates <- headline
receipt$decision_definition <- "Changed handle status, any string, or exact-suffix status/flags/tier/fallback/mask/rank/range_ok/a/e/kind/answered numeric fields; excludes continuous rank1.* estimates and c2_summaries."
receipt$comparison_alignment <- "Names and vector positions; unequal lengths and unmatched values reported separately. Changed candidate order can compare different windows."
jsonlite::write_json(receipt,file.path(B,"fallback_result.json"),pretty=TRUE,auto_unbox=TRUE,digits=17,na="string")
print(receipt[c("cases","compared_values","numeric_bit_mismatches","max_relative","relative_at","max_absolute","max_absolute_at_expected_zero","nonfinite_mismatches","decision_flip_cases","shape_change_cases","unmatched_numeric_values")])
