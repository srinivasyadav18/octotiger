//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <Kokkos_Core.hpp>
#include "octotiger/geometry.hpp"          // todo reasonably include this from physics
#include "octotiger/unitiger/hydro.hpp"    // todo reasonably include this from physics
#include "octotiger/unitiger/physics.hpp"

#include <hpx/timing/high_resolution_timer.hpp>

class [[nodiscard]] scoped_timer
{
public:
    scoped_timer(std::string const& label)
      : label(label)
      , timer() {}
    ~scoped_timer() {
        std::ostringstream s;
        s << label << ": " << timer.elapsed() << " seconds" << std::endl;
        std::cerr << s.str();
    }

private:
    std::string label;
    hpx::util::high_resolution_timer timer;
};

template <typename...>
struct WhichType;
		
namespace octotiger {

template <int NDIM>
static constexpr auto q_lowest_dimension_length = NDIM == 1 ? 3 : (NDIM == 2 ? 9 : 27);

//TODO avoid these duplicates
template <int NDIM>
static KOKKOS_INLINE_FUNCTION int flip_dim_device(const int d, int flip_dim) {
    int dims[NDIM];
    int k = d;
    for (int dim = 0; dim < NDIM; dim++) {
        dims[dim] = k % 3;
        k /= 3;
    }
    k = 0;
    dims[flip_dim] = 2 - dims[flip_dim];
    for (int dim = 0; dim < NDIM; dim++) {
        k *= 3;
        k += dims[NDIM - 1 - dim];
    }
    return k;
}

// max function
template<class T> 
KOKKOS_INLINE_FUNCTION const T& max(const T& a, const T& b)
{
    return (a < b) ? b : a;
}

template <int NDIM, int INX, typename twoDimensionalView, typename threeDimensionalView>
// F should be 0-initialized, things will be added to it and returned
// all other arguments are read-only
safe_real flux_kokkos(const int angmom_count, const int angmom_index,
    const twoDimensionalView U, const threeDimensionalView Q,
    threeDimensionalView F, const twoDimensionalView X,
    safe_real omega) {

    auto sTimer = scoped_timer("flux_kokkos");

    using geo = const cell_geometry<NDIM, INX>;

    static const auto nf = physics<NDIM>::field_count();

    // this is assuming that the indeces have the same length for each dimension
    // which we can do according to Dominic
    auto indices_size = static_cast<long int>(geo().get_indexes(3, geo::face_pts()[0][0]).size());

    // -------------------------------------------------------------------------------------------
    // this was asserts and preparations for reducing the size of the third dimension for fluxes and
    // F currently unused - third dimension is still == geo.H_N3
    //
    static const auto numGhostCellsUnitiger = geo::H_BW;    // in each direction
    static_assert(geo::H_BW == 3);

	// one ghost cell less for the flux calculation in one direction, 
	// assert all assumptions
    assert(NDIM == 3 && indices_size == (INX + 1) * INX * INX ||
        NDIM == 2 && indices_size == (INX + 1) * INX);
    // if (indices_size > geo::H_N3) throw std::runtime_error("indices size mismatch");
    static_assert((NDIM == 3 && geo::H_N3 == (geo::H_NX * geo::H_NX * geo::H_NX)) ||
        (NDIM == 2 && geo::H_N3 == (geo::H_NX * geo::H_NX)));
    static_assert((INX + 1) == (geo::H_NX - numGhostCellsUnitiger - numGhostCellsUnitiger + 1));
	
	// make it one larger, to allow for easier indexing
	auto flux_size_dim_3 = static_cast<int>(std::pow(INX + 1, NDIM));
    //--------------------------------------------------------------------------------------------

    Kokkos::View<safe_real[NDIM][geo::H_N3][nf][geo::NFACEDIR]> fluxes(
        Kokkos::ViewAllocateWithoutInitializing("fluxes"));

    Kokkos::View<int**> kokkosIndices(
        Kokkos::ViewAllocateWithoutInitializing("indices"), NDIM, indices_size);
	auto kokkosIhost = Kokkos::create_mirror_view(typename Kokkos::DefaultHostExecutionSpace::memory_space(), kokkosIndices);

    Kokkos::parallel_for("init_I", Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>({0, 0}, {indices_size, NDIM}),
        KOKKOS_LAMBDA(const int i, const int dim) { kokkosIhost(dim, i) = geo().get_indexes(3, geo::face_pts()[dim][0])[i]; });
    Kokkos::fence();
	Kokkos::deep_copy(kokkosIndices, kokkosIhost);

	safe_real amax = 0.0;
	
    // auto policy = Kokkos::MDRangePolicy<Kokkos::Experimental::HPX, Kokkos::Rank<2>>
    auto policy =
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {static_cast<int>(indices_size), NDIM});
        // Kokkos::MDRangePolicy<Kokkos::Serial, Kokkos::Rank<2>>({0, 0}, {static_cast<int>(indices_size), NDIM});

    Kokkos::fence();

    Kokkos::parallel_reduce("compute fluxes", policy,
        KOKKOS_LAMBDA(const int indexIteration, const int dim, safe_real& maxAmax) {
            // ----------------------------------------------------------------------------------------//
            // these arrays/references to static constexpr host variables are not allowed in device code
            // static constexpr auto faces = geo::face_pts(); 
            // static constexpr auto weights = geo::face_weight();
            // static constexpr auto xloc = geo::xloc();
            // static constexpr auto kdelta = geo::kronecker_delta();
            // hacky workaround here: copy them into the kernel, index the relevant dimensionality
            // TODO make this nicer
            static constexpr int lower_face_members_device[3][3][9] = { { { 0 } }, 
                        { { 3, 0, 6 }, { 1, 0, 2 } }, 
                        { { 12, 0, 3, 6, 9, 15, 18, 21, 24 }, { 10, 0, 1, 2, 9, 11,
			                        18, 19, 20 }, { 4, 0, 1, 2, 3, 5, 6, 7, 8 } } };
			static constexpr auto faces = lower_face_members_device[NDIM-1];

            static constexpr safe_real quad_weights_device[3][9] = { { 1.0 }, { 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0 }, { 16. / 36., 1. / 36., 4. / 36., 1. / 36., 4. / 36., 4.
                    / 36., 1. / 36., 4. / 36., 1. / 36. } };
            static constexpr auto weights = quad_weights_device[NDIM-1];

            static constexpr int face_locs_device[3][27][3] = {
            /**/{ { -1 }, { 0 }, { 1 } },

            /**/{
            /**/{ -1, -1 }, { +0, -1 }, { +1, -1 },
            /**/{ -1, +0 }, { +0, +0 }, { +1, +0 },
            /**/{ -1, +1 }, { +0, +1 }, { +1, +1 } },

            /**/{
            /**/{ -1, -1, -1 }, { +0, -1, -1 }, { +1, -1, -1 },
            /**/{ -1, +0, -1 }, { +0, +0, -1 }, { 1, +0, -1 },
            /**/{ -1, +1, -1 }, { +0, +1, -1 }, { +1, +1, -1 },
            /**/{ -1, -1, +0 }, { +0, -1, +0 }, { +1, -1, +0 },
            /**/{ -1, +0, +0 }, { +0, +0, +0 }, { +1, +0, +0 },
            /**/{ -1, +1, +0 }, { +0, +1, +0 }, { +1, +1, +0 },
            /**/{ -1, -1, +1 }, { +0, -1, +1 }, { +1, -1, +1 },
            /**/{ -1, +0, +1 }, { +0, +0, +1 }, { +1, +0, +1 },
            /**/{ -1, +1, +1 }, { +0, +1, +1 }, { +1, +1, +1 } } };
            static constexpr auto xloc = face_locs_device[NDIM-1];

            static constexpr int kdeltas_device[3][3][3][3] = { { { { } } }, { { { 0, 1 }, { -1, 0 } } }, { { { 0, 0, 0 }, { 0, 0, 1 }, { 0, -1, 0 } }, { { 0, 0, -1 }, { 0, 0,
                    0 }, { 1, 0, 0 } }, { { 0, 1, 0 }, { -1, 0, 0 }, { 0, 0, 0 } } } };
            static constexpr auto kdelta = kdeltas_device[NDIM - 1];

            // a copy of geo::H_DN, for device code
            static constexpr int H_DNX_device = NDIM == 3 ? H_NX * H_NX : (NDIM == 2 ? H_NX : 1);
            static constexpr int H_DNY_device = NDIM == 3 ? H_NX : 1;
            static constexpr int H_DNZ_device = 1;
           	static constexpr int H_DN_device[3] = { H_DNX_device, H_DNY_device, H_DNZ_device };

            static const auto nf_device = physics<NDIM>::field_count();  //--expt-relaxed-constexpr
            printf("nf %d ", nf_device);
            //--------------------------------------------------------------//

            const auto dx = X(0, H_DNX_device) - X(0, 0);

            safe_real this_flux[nf_device];

            // const auto& indices = geo::get_indexes(3, geo::face_pts()[dim][0]);
            const auto& i = kokkosIndices(dim,indexIteration);
            safe_real ap = 0.0, am = 0.0;
            safe_real this_ap, this_am;

            for (int fi = 0; fi < geo::NFACEDIR; fi++) {
				const auto d = faces[dim][fi];
                // UR0, UL0 are not needed for now
                // auto UR0 = Kokkos::subview(U, Kokkos::ALL, i);
                // auto UL0 = Kokkos::subview(U, Kokkos::ALL, i - geo::H_DN[dim]);
                auto UR = Kokkos::subview(Q, Kokkos::ALL, i, d);
                auto UL = Kokkos::subview(Q, Kokkos::ALL, i - H_DN_device[dim], flip_dim_device<NDIM>(d, dim));

                Kokkos::Array<safe_real, NDIM> vg;
                if CONSTEXPR(NDIM > 1) {
                    vg[0] = -omega * (X(1, i) + 0.5 * xloc[d][1] * dx);
                    vg[1] = +omega * (X(0, i) + 0.5 * xloc[d][0] * dx);
                    if CONSTEXPR(NDIM == 3) {
                        vg[2] = 0.0;
                    }
                }
                else {
                    vg[0] = 0.0;
                }
                // printf("%d ", fi);

                // physics<NDIM>::flux(UL, UR, UL0, UR0, this_flux, dim, this_am, this_ap, vg, dx);
                physics<NDIM>::flux(UL, UR, this_flux, dim, this_am, this_ap, vg, dx);  //--expt-relaxed-constexpr

                am = std::min(am, this_am); //--expt-relaxed-constexpr
                ap = std::max(ap, this_ap);  //--expt-relaxed-constexpr
                for (int f = 0; f < nf_device; f++) {
                    fluxes(dim, i, f, fi) = this_flux[f];
                }
            }
            
            maxAmax = std::max(ap, safe_real(-am));  //--expt-relaxed-constexpr

                // field update from fluxes
			for (int f = 0; f < nf_device; f++) {
                    F(dim, f, i) = 0.0;
                for (int fi = 0; fi < geo::NFACEDIR; fi++) {
                        const auto& w = weights[fi];
                    F(dim, f, i) += w * fluxes(dim, i, f, fi);
					}
				}
                // angular momentum update from linear momentum
			for (int angmom_pair = 0; angmom_pair < angmom_count; angmom_pair++) {
                const int sx_i = angmom_index + angmom_pair * (NDIM + geo::NANGMOM);
				const int zx_i = sx_i + NDIM;
                for (int n = 0; n < geo::NANGMOM; n++) {
                    F(dim, zx_i + n, i) = fluxes(dim, i, zx_i + n, 0);

					for (int m = 0; m < NDIM; m++) {
						if (dim != m) {
							for (int l = 0; l < NDIM; l++) {
                                for (int fi = 0; fi < geo::NFACEDIR; fi++) {
									const auto d = faces[dim][fi];
                                        F(dim, zx_i + n, i) += weights[fi] * kdelta[n][m][l] *
                                        xloc[d][m] * 0.5 * dx * fluxes(dim, i, sx_i + l, fi);
								}
							}
						}
					}
				}
			}
        },
        Kokkos::Max<safe_real>(amax));

    Kokkos::fence();
		return amax;
}    // flux_kokkos

// template <int NDIM, int INX>
// safe_real flux_kokkos_hpx(const int angmom_count, const int angmom_index,
//     const Kokkos::View<safe_real**> U, const Kokkos::View<safe_real***> Q,
//     Kokkos::View<safe_real***> F, const Kokkos::View<safe_real**> X, safe_real omega) {
//     auto a =
//         hpx::async([=]() -> safe_real { return flux_kokkos( angmom_count, angmom_index, U, Q, F, X, omega); });
//     return a.get();
// }

template <int NDIM, int INX>
safe_real compute_flux_kokkos(hydro_computer<NDIM, INX>& computer,
    Kokkos::View<safe_real**, Kokkos::HostSpace>& kokkosUhost,
    Kokkos::View<safe_real***, Kokkos::HostSpace> kokkosQhost,
    Kokkos::View<safe_real***, Kokkos::HostSpace> kokkosFhost,
    Kokkos::View<safe_real**, Kokkos::HostSpace> kokkosXhost,
    safe_real omega
    ) {

    auto kokkosF = Kokkos::create_mirror_view(typename Kokkos::DefaultExecutionSpace::memory_space(), kokkosFhost);
    auto kokkosU = Kokkos::create_mirror_view(typename Kokkos::DefaultExecutionSpace::memory_space(), kokkosUhost);
    auto kokkosX = Kokkos::create_mirror_view(typename Kokkos::DefaultExecutionSpace::memory_space(), kokkosXhost);
    auto kokkosQ = Kokkos::create_mirror_view(typename Kokkos::DefaultExecutionSpace::memory_space(), kokkosQhost);

    Kokkos::fence();

    Kokkos::deep_copy(kokkosU, kokkosUhost);
    Kokkos::deep_copy(kokkosX, kokkosXhost);
    Kokkos::deep_copy(kokkosQ, kokkosQhost);

    Kokkos::fence();
    
    auto amax = octotiger::flux_kokkos<NDIM, INX>(computer.getAngMomCount(),
        computer.getAngMomIndex(), kokkosU, kokkosQ, kokkosF, kokkosX, omega);

    Kokkos::fence();

    Kokkos::deep_copy(kokkosFhost, kokkosF);

    Kokkos::fence();

    return amax;
}

template <int NDIM, int INX>
safe_real compute_flux_kokkos(hydro_computer<NDIM, INX>& computer,
    const std::vector<std::vector<safe_real>>& U0, const std::vector<std::vector<safe_real>>& U,
    const std::vector<std::vector<std::array<safe_real, q_lowest_dimension_length<NDIM>>>>& q,
    std::vector<std::vector<std::vector<safe_real>>>& F,
    const std::vector<std::vector<safe_real>>& X, const safe_real omega, const int nf,
    const int H_N3) {
    // do this iteration in the kokkosified version
    // F is an output of the flux kernel, zero-initialize
    Kokkos::View<safe_real***, Kokkos::HostSpace> kokkosFhost("flux", NDIM, nf, H_N3);    

    Kokkos::View<safe_real**, Kokkos::HostSpace> kokkosUhost(Kokkos::ViewAllocateWithoutInitializing("state"), nf, H_N3);
    Kokkos::View<safe_real**, Kokkos::HostSpace> kokkosU0host(
        Kokkos::ViewAllocateWithoutInitializing("initial state"), nf, H_N3);
    Kokkos::parallel_for("init_U",
        Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>(
            {0, 0}, {nf, H_N3}),
        KOKKOS_LAMBDA(int j, int k) {
            kokkosUhost(j, k) = U[j][k];
            kokkosU0host(j, k) = U0[j][k];
        });

    Kokkos::View<safe_real**, Kokkos::HostSpace> kokkosXhost(
        Kokkos::ViewAllocateWithoutInitializing("cell center coordinates"), NDIM, H_N3);
    Kokkos::parallel_for("init_X",
        Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>(
            {0, 0}, {NDIM, H_N3}),
        KOKKOS_LAMBDA(int i, int k) { kokkosXhost(i, k) = X[i][k]; });

    Kokkos::View<safe_real* * [q_lowest_dimension_length<NDIM>], Kokkos::HostSpace> kokkosQhost(
        Kokkos::ViewAllocateWithoutInitializing("reconstruction"), nf, H_N3);
    // have this run in serial, takes forever otherwise
    Kokkos::parallel_for("init_Q",
        Kokkos::MDRangePolicy<Kokkos::Serial, Kokkos::Rank<3>>(
            {0, 0, 0}, {nf, H_N3, q_lowest_dimension_length<NDIM>}),
        KOKKOS_LAMBDA(int i, int j, int k) { kokkosQhost(i, j, k) = q[i][j][k]; }); //--expt-relaxed-constexpr

    // std::cout << kokkosQ.extent(0) << " " << kokkosQ.extent(1) << " " << kokkosQ.extent(2) << "Q"
    // << std::endl;
    Kokkos::fence();

    auto amax = compute_flux_kokkos(computer, kokkosUhost, kokkosQhost, kokkosFhost, kokkosXhost, omega);

    Kokkos::fence();

    // copy back the values obtained for F
    for (int i = 0; i < NDIM; ++i) {
        for (int j = 0; j < nf; ++j) {
            for (int k = 0; k < H_N3; ++k) {
                F[i][j][k] = kokkosFhost(i, j, k);
            }
        }
    }

    return amax;
}

}    // namespace octotiger
