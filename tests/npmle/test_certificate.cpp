// A final certificate margin must not change the iteration path or accept a poor fit.
#include "../../bunch/core/src/npmle_fit.hpp"
#include <cassert>
int main() {
  sbnpmle::Likelihood model;
  model.m=2; model.n=2; model.effective_n=1;
  model.design={1,0,0,1}; model.atom={.1,.1}; model.penalty={0,0,0,0};
  const std::vector<double> initial={.5,.5};
  model.counts={.5+8e-8,.5-8e-8};
  auto near=sbnpmle::optimize(model,0,initial,0,1e-7);
  assert(near.kkt>1e-7 && near.kkt<1e-6 && near.converged);
  assert(near.p==initial);
  model.counts={.5+8e-6,.5-8e-6};
  auto far=sbnpmle::optimize(model,0,initial,0,1e-7);
  assert(far.kkt>1e-6 && !far.converged);
  assert(far.p==initial);
}
