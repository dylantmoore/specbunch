#include <stdexcept>
extern "C" void spectral_r_stop(int) { throw std::runtime_error("Disabled standalone solver exit"); }
