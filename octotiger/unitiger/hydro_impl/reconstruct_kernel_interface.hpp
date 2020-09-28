#pragma once

//#define TVD_TEST

#include <array>
#include <vector>

#include "octotiger/cuda_util/cuda_global_def.hpp"
#include "octotiger/hydro_defs.hpp"
#include "octotiger/unitiger/hydro.hpp"
#include "octotiger/unitiger/safe_real.hpp"

#include "octotiger/unitiger/physics.hpp"
#include "octotiger/unitiger/physics_impl.hpp"

#include "octotiger/unitiger/hydro_impl/flux_kernel_interface.hpp"

template <typename T>
CUDA_CALLABLE_METHOD inline T copysign_wrapper(const T& tmp1, const T& tmp2) {
    return std::copysign(tmp1, tmp2);
}
template <typename T>
CUDA_CALLABLE_METHOD inline T abs_wrapper(const T& tmp1) {
    return std::abs(tmp1);
}
template <typename T>
CUDA_CALLABLE_METHOD inline T minmod_wrapper(const T& a, const T& b) {
    return (copysign_wrapper<T>(0.5, a) + copysign_wrapper<T>(0.5, b)) *
        min_wrapper<T>(abs_wrapper<T>(a), abs_wrapper<T>(b));
}
template <typename T>
CUDA_CALLABLE_METHOD inline T minmod_theta_wrapper(const T& a, const T& b, const T& c) {
    return minmod_wrapper<T>(c * minmod_wrapper<T>(a, b), 0.5 * (a + b));
}

void reconstruct_ppm_experimental(std::vector<std::vector<safe_real>>& q,
    const std::vector<safe_real>& u, bool smooth, bool disc_detect,
    const std::vector<std::vector<double>>& disc);

const hydro::recon_type<NDIM>& reconstruct_experimental(const hydro::state_type& U_,
    const hydro::x_type& X, safe_real omega, const size_t nf_, const int angmom_index_,
    const std::vector<bool>& smooth_field_, const std::vector<bool>& disc_detect_);
