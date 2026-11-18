





#include <cstring>
#include <memory>
#include <string>
#include <vector>
#include "mean_engine.hpp"
#include "public_result.hpp"
#include "sbc_internal.hpp"
#include "specbunch.h"

namespace {

int bmean_impl(const sb_data_xy* data, const sb_opts_bmean* opts_in, sb_result** out) {
  if (!data || !out) return sbc::fail(SBC_E_NULL, "data or out is NULL");
  *out = nullptr;
  sb_opts_bmean o; sb_opts_bmean_init(&o);
  if (opts_in) {
    if (int rc = sbc::check_opts_header(opts_in, 8, SB_ABI_VERSION)) return rc;
    sb_opts_bmean tmp; sb_opts_bmean_init(&tmp);
    std::memcpy(&tmp, opts_in, std::min<size_t>(opts_in->struct_size, sizeof tmp));
    o = tmp; o.struct_size = sizeof o;
  }
  if (o.reserved_band < 0 || o.reserved_band > 2 || o.reserved_visibility < 0 || o.reserved_visibility > 2 || o.R < 19)
    return sbc::fail(SBC_E_INPUT, "reserved_band and reserved_visibility must be 0..2, R >= 19");
  if (data->n == 0) return sbpublic::insufficient(out,0,"buncher_mean",o.threads);
  if (!data->z || !data->y) return sbc::fail(SBC_E_INPUT, "sample pointers are NULL");
  const bool weighted = data->w != nullptr;
  std::vector<double> V, Y, W, T;
  V.reserve(data->n); Y.reserve(data->n); T.reserve(data->n); if (weighted) W.reserve(data->n);
  for (size_t i = 0; i < data->n; ++i) {
    const double wi = weighted ? data->w[i] : 1.0;
    if (!std::isfinite(data->z[i]) || !std::isfinite(data->y[i]) || !(std::isfinite(wi) && wi > 0)) continue;
    V.push_back(data->z[i]); Y.push_back(data->y[i]); T.push_back(1.0); if (weighted) W.push_back(wi);
  }
  if(V.size()<100)return sbpublic::insufficient(out,V.size(),"buncher_mean",o.threads);
  auto* R = new sbc::Result();
  std::unique_ptr<sbc::Result> guard(R);
  p2::Result r;
  {
    sbc::ThreadGuard tg(o.threads);
    R->set("threads_used", tg.count());
    std::vector<double> om, cw;
    p2::cell_weights(V.size(), weighted ? W.data() : nullptr, T.data(), weighted ? 1 : 0, om, cw);
    p2::Config cf; cf.logu = 1; cf.seed = o.seed; cf.R = o.R; cf.reserved_band = o.reserved_band; cf.reserved_visibility = o.reserved_visibility;
    r = p2::analyze(V, Y, om, cw, cf);
  }
  const double sig = r.sigma_hat, D = r.D, NaN = std::numeric_limits<double>::quiet_NaN();
  auto put = [&](const std::string& p, const p2::Est& e, bool headline) {
    const double lo = e.a >= 0 ? (e.a + 1) * D : NaN, hi = e.e >= 0 ? (e.e + 1) * D : NaN;
    R->set(p + "mu", e.est); R->set(p + "status", e.status); R->sets(p + "reason", e.status == 0 ? "" : p2::STATUS_NAME[e.status]);
    R->set(p + "mask", e.mask); R->set(p + "band_lo", lo); R->set(p + "band_hi", hi); R->set(p + "band_lo_u", lo * sig);
    R->set(p + "band_hi_u", hi * sig); R->set(p + "m", e.a >= 0 ? e.e - e.a + 1 : 0); R->set(p + "p", e.p); R->set(p + "se", e.se);
    R->set(p + "den", e.den); R->set(p + "rank", e.rank); R->set(p + "pw_mean", e.pw_mean); R->set(p + "pw_median", e.pw_med);
    R->set(p + "imag", e.imag);
    (void)headline;
  };
  put("", r.cos_chk, true);

  R->status = r.cos_chk.status;
  R->sets("se_note", "se is the standard error of the uncorrected least-squares ratio on the selected band; it is NOT a confidence interval");
  R->set("n", (double)V.size()); R->set("n_input", (double)data->n); R->set("weighted", weighted ? 1.0 : 0.0);
  R->set("n_pos", (double)r.n_pos); R->set("W", r.W); R->set("W2", r.W2);
  R->set("sigma_hat", r.sigma_hat); R->set("sigma_fallback", r.sigma_fallback); R->set("sig_filt", r.sig_filt); R->set("sig_iter", r.sig_iter);
  R->set("tc", r.tc); R->set("tc_u", r.tc * sig); R->set("nvis", r.nvis); R->set("frequency_band_visible", r.frequency_band_visible); R->set("refine", r.scan_refine);
  R->set("nband", r.nband); R->set("chol", r.chol_jitter); R->set("q1", r.q1); R->set("q2", r.q2); R->set("ng", r.ng);
  R->set("ysd", r.ysd); R->set("ymin", r.ymin); R->set("ymax", r.ymax); R->set("D", r.D);
  R->set("pass_cos", r.pass_cos); R->set("pass_cpx", r.pass_cpx); R->set("pfull_cos", r.pfull_cos);
  R->set("seed", (double)o.seed); R->set("R", o.R); R->set("reserved_band", o.reserved_band); R->set("reserved_visibility", o.reserved_visibility);
  const char* columns[]={"rank","band_lo","band_hi","band_lo_scaled","band_hi_scaled","mu","pointwise_median","se","fit_p","coherent","poor_fit","selected"};
  for(int j=0;j<12;++j) { std::vector<double> v;for(const auto &e:r.windows) {
    double lo=(e.a+1)*D,hi=(e.e+1)*D;
    double q[]={double(e.rank),lo,hi,lo*sig,hi*sig,e.est,e.pw_med,e.se,e.p,double(e.status==0),double(e.p<0.1),double(e.rank==r.cos_chk.rank && r.cos_chk.status==0)};
    v.push_back(q[j]);}R->setv(std::string("windows.")+columns[j],std::move(v));}
  R->set("flags",r.cos_chk.status==0 && r.cos_chk.p<0.1?32:0);
  R->sets("flags_str",r.cos_chk.status==0 && r.cos_chk.p<0.1?"poor_fit":"");
  sbpublic::finish(*R,"buncher_mean");
  R->finish_names();
  *out = reinterpret_cast<sb_result*>(guard.release());
  return SBC_OK;
}

}  

extern "C" {
SB_API void sb_opts_bmean_init(sb_opts_bmean* o) {
  if (!o) return;
  std::memset(o, 0, sizeof *o);
  o->struct_size = sizeof(sb_opts_bmean); o->abi_version = SB_ABI_VERSION;
  o->threads = 0; o->reserved_band = 2; o->reserved_visibility = 1; o->R = 499; o->seed = 1;
}
SB_API int sb_bmean(const sb_data_xy* data, const sb_opts_bmean* opts, sb_result** out) {
  return sbc::guarded([&] { return bmean_impl(data, opts, out); });
}
}
