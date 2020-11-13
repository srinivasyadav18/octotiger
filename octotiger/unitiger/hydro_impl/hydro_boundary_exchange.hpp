#include "octotiger/grid.hpp"

#include <aligned_buffer_util.hpp>
#include <buffer_manager.hpp>
#ifdef OCTOTIGER_HAVE_CUDA
#include <cuda_buffer_util.hpp>
#include <cuda_runtime.h>
#include <stream_manager.hpp>
#include "octotiger/cuda_util/cuda_helper.hpp"
#endif
#include "octotiger/util/vec_scalar_host_wrapper.hpp"

void complete_hydro_amr_boundary_cpu(const double dx, const bool energy_only,
    const std::vector<std::vector<real>>& ushad, const std::vector<std::atomic<int>>& is_coarse,
    const std::array<double, NDIM>& xmin, std::vector<std::vector<real>>& u);
void complete_hydro_amr_boundary_vc(const double dx, const bool energy_only,
    const std::vector<std::vector<real>>& Ushad, const std::vector<std::atomic<int>>& is_coarse,
    const std::array<double, NDIM>& xmin, std::vector<std::vector<double>>& U);
void launch_complete_hydro_amr_boundary_cuda(
    stream_interface<hpx::cuda::experimental::cuda_executor, pool_strategy>& executor, double dx,
    bool energy_only, const std::vector<std::vector<real>>& ushad,
    const std::vector<std::atomic<int>>& is_coarse, const std::array<double, NDIM>& xmin,
    std::vector<std::vector<real>>& u);

CUDA_GLOBAL_METHOD inline double minmod_cuda(double a, double b) {
    return (copysign(0.5, a) + copysign(0.5, b)) * std::min(std::abs(a), abs(b));
}

CUDA_GLOBAL_METHOD inline double minmod_cuda_theta(double a, double b, double c) {
    return minmod_cuda(c * minmod_cuda(a, b), 0.5 * (a + b));
}

CUDA_GLOBAL_METHOD inline double limiter(const double a, const double b) {
    return minmod_cuda_theta(a, b, 64. / 37.);
}

template <typename T, typename mask_t, typename index_t>
CUDA_GLOBAL_METHOD inline void complete_hydro_amr_boundary_inner_loop(const double dx, const bool energy_only,
    const double* __restrict__ unified_ushad, const int* __restrict__ coarse,
    const double* __restrict__ xmin, double* __restrict__ unified_uf, const int i0, const int j0,
    const int k0, const int nfields, const mask_t mask, const index_t k, const int iii0) {
    const int field_offset = HS_N3 * 8;
    for (int ir = 0; ir < 2; ir++) {
        for (int jr = 0; jr < 2; jr++) {
            for (int kr = 0; kr < 2; kr++) {
                #pragma unroll
                for (int f = 0; f < nfields; f++) {
                    if (!energy_only || f == egas_i) {
                        const int oct_index = ir * 4 + jr * 2 + kr;

                        const auto is = ir % 2 ? +1 : -1;
                        const auto js = jr % 2 ? +1 : -1;
                        const auto ks = kr % 2 ? +1 : -1;
                        // const auto& u0 = Ushad[f][iii0];
                        // const auto& uc = Ushad[f];
                        const double* uc = unified_ushad + f * HS_N3;

                        const T u0 = load_value<T>(unified_ushad, f * HS_N3 + iii0);
                        const T uc_x = load_value<T>(uc, iii0 + is * HS_DNX);
                        const T uc_y = load_value<T>(uc, iii0 + js * HS_DNY);
                        const T uc_z = load_value<T>(uc, iii0 + ks * HS_DNZ);
                        const T uc_x_neg = load_value<T>(uc, iii0 - is * HS_DNX);
                        const T uc_y_neg = load_value<T>(uc, iii0 - js * HS_DNY);
                        const T uc_z_neg = load_value<T>(uc, iii0 - ks * HS_DNZ);

                        const auto s_x = limiter_wrapper(uc_x - u0, u0 - uc_x_neg);
                        const auto s_y = limiter_wrapper(uc_y - u0, u0 - uc_y_neg);
                        const auto s_z = limiter_wrapper(uc_z - u0, u0 - uc_z_neg);

                        const T uc_xy = load_value<T>(uc, iii0 + is * HS_DNX + js * HS_DNY);
                        const T uc_xz = load_value<T>(uc, iii0 + is * HS_DNX + ks * HS_DNZ);
                        const T uc_yz = load_value<T>(uc, iii0 + js * HS_DNY + ks * HS_DNZ);
                        const T uc_xy_neg = load_value<T>(uc, iii0 - is * HS_DNX - js * HS_DNY);
                        const T uc_xz_neg = load_value<T>(uc, iii0 - is * HS_DNX - ks * HS_DNZ);
                        const T uc_yz_neg = load_value<T>(uc, iii0 - js * HS_DNY - ks * HS_DNZ);

                        const auto s_xy = limiter_wrapper(uc_xy - u0, u0 - uc_xy_neg);
                        const auto s_xz = limiter_wrapper(uc_xz - u0, u0 - uc_xz_neg);
                        const auto s_yz = limiter_wrapper(uc_yz - u0, u0 - uc_yz_neg);

                        const T uc_xyz =
                            load_value<T>(uc, iii0 + is * HS_DNX + js * HS_DNY + ks * HS_DNZ);
                        const T uc_xyz_neg =
                            load_value<T>(uc, iii0 - is * HS_DNX - js * HS_DNY - ks * HS_DNZ);
                        const auto s_xyz = limiter_wrapper(uc_xyz - u0, u0 - uc_xyz_neg);
                        // auto &uf = Uf[f][iii0][ir][jr][kr];
                        // T uf = ;
                        // load_value<T>(unified_uf, f * field_offset + 8 * iii0 + oct_index);
                        T uf = u0;
                        uf += (9.0 / 64.0) * (s_x + s_y + s_z);
                        uf += (3.0 / 64.0) * (s_xy + s_yz + s_xz);
                        uf += (1.0 / 64.0) * s_xyz;
                        store_value<T>(unified_uf, f * field_offset + iii0 + oct_index * HS_N3, uf);
                    }
                }
                if (!energy_only) {
                    const int oct_index = ir * 4 + jr * 2 + kr;
                    const auto i1 = 2 * i0 - H_BW + ir;
                    const auto j1 = 2 * j0 - H_BW + jr;
                    const auto k1 = 2 * (k0 + k) - H_BW + kr;
                    const auto x = (i1) *dx + xmin[XDIM];
                    const auto y = (j1) *dx + xmin[YDIM];
                    const auto z = (k1) *dx + xmin[ZDIM];

                    const T sx_val =
                        load_value<T>(unified_uf, sx_i * field_offset + iii0 + oct_index * HS_N3);
                    const T sy_val =
                        load_value<T>(unified_uf, sy_i * field_offset + iii0 + oct_index * HS_N3);
                    const T sz_val =
                        load_value<T>(unified_uf, sz_i * field_offset + iii0 + oct_index * HS_N3);

                    T lx_val =
                        load_value<T>(unified_uf, lx_i * field_offset + iii0 + oct_index * HS_N3);
                    lx_val -= y * sz_val - z * sy_val;
                    store_value<T>(
                        unified_uf, lx_i * field_offset + iii0 + oct_index * HS_N3, lx_val);

                    T ly_val =
                        load_value<T>(unified_uf, ly_i * field_offset + iii0 + oct_index * HS_N3);
                    ly_val += x * sz_val - z * sx_val;
                    store_value<T>(
                        unified_uf, ly_i * field_offset + iii0 + oct_index * HS_N3, ly_val);

                    T lz_val =
                        load_value<T>(unified_uf, lz_i * field_offset + iii0 + oct_index * HS_N3);
                    lz_val -= x * sy_val - y * sx_val;
                    store_value<T>(
                        unified_uf, lz_i * field_offset + iii0 + oct_index * HS_N3, lz_val);
                }
            }
        }
    }
    if (!energy_only) {
        T zx = 0, zy = 0, zz = 0, rho = 0;
        for (int ir = 0; ir < 2; ir++) {
            for (int jr = 0; jr < 2; jr++) {
                #pragma unroll
                for (int kr = 0; kr < 2; kr++) {
                    const int oct_index = ir * 4 + jr * 2 + kr;
                    T lx_val =
                        load_value<T>(unified_uf, lx_i * field_offset + iii0 + oct_index * HS_N3);
                    T ly_val =
                        load_value<T>(unified_uf, ly_i * field_offset + iii0 + oct_index * HS_N3);
                    T lz_val =
                        load_value<T>(unified_uf, lz_i * field_offset + iii0 + oct_index * HS_N3);
                    zx += lx_val / 8.0;
                    zy += ly_val / 8.0;
                    zz += lz_val / 8.0;
                    //			rho += Uf[rho_i][iii0][ir][jr][kr] / 8.0;
                }
            }
        }
        for (int ir = 0; ir < 2; ir++) {
            for (int jr = 0; jr < 2; jr++) {
                #pragma unroll
                for (int kr = 0; kr < 2; kr++) {
                    //					const auto factor =
                    // Uf[rho_i][iii0][ir][jr][kr]
                    /// rho;
                    const auto factor = 1.0;
                    const int oct_index = ir * 4 + jr * 2 + kr;
                    zx *= factor;
                    zy *= factor;
                    zz *= factor;
                    store_value<T>(unified_uf, lx_i * field_offset + iii0 + oct_index * HS_N3, zx);
                    store_value<T>(unified_uf, ly_i * field_offset + iii0 + oct_index * HS_N3, zy);
                    store_value<T>(unified_uf, lz_i * field_offset + iii0 + oct_index * HS_N3, zz);
                }
            }
        }
        for (int ir = 0; ir < 2; ir++) {
            for (int jr = 0; jr < 2; jr++) {
                #pragma unroll
                for (int kr = 0; kr < 2; kr++) {
                    const int oct_index = ir * 4 + jr * 2 + kr;
                    const auto i1 = 2 * i0 - H_BW + ir;
                    const auto j1 = 2 * j0 - H_BW + jr;
                    // TODO(daissgr) Replace with vectorization
                    const auto k1 = 2 * (k0 + k) - H_BW + kr;
                    // std::cout << k1 << std::endl;
                    const auto x = (i1) *dx + xmin[XDIM];
                    const auto y = (j1) *dx + xmin[YDIM];
                    const auto z = (k1) *dx + xmin[ZDIM];
                    // std::cout << z << std::endl;

                    const T sx_val =
                        load_value<T>(unified_uf, sx_i * field_offset + iii0 + oct_index * HS_N3);
                    const T sy_val =
                        load_value<T>(unified_uf, sy_i * field_offset + iii0 + oct_index * HS_N3);
                    const T sz_val =
                        load_value<T>(unified_uf, sz_i * field_offset + iii0 + oct_index * HS_N3);

                    T lx_val =
                        load_value<T>(unified_uf, lx_i * field_offset + iii0 + oct_index * HS_N3);
                    lx_val += y * sz_val - z * sy_val;
                    select_wrapper<T>(lx_val, mask, lx_val, T(0));
                    store_value<T>(
                        unified_uf, lx_i * field_offset + iii0 + oct_index * HS_N3, lx_val);

                    T ly_val =
                        load_value<T>(unified_uf, ly_i * field_offset + iii0 + oct_index * HS_N3);
                    ly_val -= x * sz_val - z * sx_val;
                    select_wrapper<T>(ly_val, mask, ly_val, T(0));
                    store_value<T>(
                        unified_uf, ly_i * field_offset + iii0 + oct_index * HS_N3, ly_val);

                    T lz_val =
                        load_value<T>(unified_uf, lz_i * field_offset + iii0 + oct_index * HS_N3);
                    lz_val += x * sy_val - y * sx_val;
                    select_wrapper<T>(lz_val, mask, lz_val, T(0));
                    store_value<T>(
                        unified_uf, lz_i * field_offset + iii0 + oct_index * HS_N3, lz_val);
                }
            }
        }
    }
}
