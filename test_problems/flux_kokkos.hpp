//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <Kokkos_Core.hpp>
#include "octotiger/geometry.hpp"          // todo reasonably include this from physics
#include "octotiger/unitiger/hydro.hpp"    // todo reasonably include this from physics
#include "octotiger/unitiger/physics.hpp"

#include <hpx/util/high_resolution_timer.hpp>

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
template <int NDIM, int INX>
// safe_real flux_kokkos(hydro_computer<NDIM, INX>& hydroComputer, const hydro::state_type &U, const
// hydro::recon_type<NDIM> &Q, hydro::flux_type &F, hydro::x_type &X, 		safe_real omega) {

// F should be 0-initialized, things will be added to it and returned
// all other arguments are read-only
safe_real flux_kokkos(const hydro_computer<NDIM, INX>& hydroComputer,
    const Kokkos::View<safe_real**> U, const Kokkos::View<safe_real***> Q,
    Kokkos::View<safe_real***> F, const Kokkos::View<safe_real**> X, safe_real omega) {

    auto sTimer = scoped_timer("flux_kokkos");

		static const cell_geometry<NDIM, INX> geo;

		static constexpr auto nf = physics<NDIM>::field_count();
		const auto& angmom_count = hydroComputer.getAngMomCount();
		const auto& angmom_index = hydroComputer.getAngMomIndex();

    // this is assuming that the indeces have the same length for each dimension
    // which we can do according to Dominic
    auto indices_size = geo.get_indexes(3, geo.face_pts()[0][0]).size();

    // -------------------------------------------------------------------------------------------
    // this was asserts and preparations for reducing the size of the third dimension for fluxes and F
    // currently unused - third dimension is still == geo.H_N3
    //
	static constexpr auto numGhostCellsUnitiger = geo.H_BW; // in each direction
	static_assert(geo.H_BW == 3);

	// one ghost cell less for the flux calculation in one direction, 
	// assert all assumptions
	assert( NDIM == 3 && indices_size == (INX+1)*INX*INX || NDIM == 2 && indices_size == (INX+1)*INX );
	// if (indices_size > geo.H_N3) throw std::runtime_error("indices size mismatch");
	static_assert(( NDIM==3 && geo.H_N3 == (geo.H_NX * geo.H_NX * geo.H_NX)) || ( NDIM==2 && geo.H_N3 == (geo.H_NX * geo.H_NX)));
	static_assert((INX+1) == (geo.H_NX - numGhostCellsUnitiger - numGhostCellsUnitiger + 1));
	
	// make it one larger, to allow for easier indexing
	auto flux_size_dim_3 = static_cast<int>(std::pow(INX + 1, NDIM));
    //--------------------------------------------------------------------------------------------

    Kokkos::View<safe_real[NDIM][nf][geo.H_N3][geo.NFACEDIR]> fluxes(
        Kokkos::ViewAllocateWithoutInitializing("fluxes"));

		static constexpr auto faces = geo.face_pts();
		static constexpr auto weights = geo.face_weight();
		static constexpr auto xloc = geo.xloc();
		static constexpr auto kdelta = geo.kronecker_delta();

    const auto dx = X(0, geo.H_DNX) - X(0, 0);

		safe_real amax = 0.0;

	
    // auto policy = Kokkos::MDRangePolicy<Kokkos::Experimental::HPX, Kokkos::Rank<2>>
    auto policy =
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {static_cast<int>(indices_size), NDIM});
        // Kokkos::MDRangePolicy<Kokkos::Serial, Kokkos::Rank<2>>({0, 0}, {static_cast<int>(indices_size), NDIM});

    Kokkos::fence();

    Kokkos::parallel_reduce("process dimensions", policy,
        // "process dimensions", NDIM,
        KOKKOS_LAMBDA(const int indexIteration, const int dim, safe_real& maxAmax) {
			// printf("%d", hpx::get_worker_thread_num());
            std::array<safe_real, nf> this_flux;

            const auto& indices = geo.get_indexes(3, geo.face_pts()[dim][0]);
            const auto& i = indices[indexIteration];
				safe_real ap = 0.0, am = 0.0;
				safe_real this_ap, this_am;
				for (int fi = 0; fi < geo.NFACEDIR; fi++) {
					const auto d = faces[dim][fi];
                // UR0, UL0 are not needed for now
                // auto UR0 = Kokkos::subview(U, Kokkos::ALL, i);
                // auto UL0 = Kokkos::subview(U, Kokkos::ALL, i - geo.H_DN[dim]);
                    auto UR = Kokkos::subview(Q, Kokkos::ALL, i, d);
                auto UL = Kokkos::subview(Q, Kokkos::ALL, i - geo.H_DN[dim], geo.flip_dim(d, dim));
                    std::array<safe_real, NDIM> vg;
                    if
                        CONSTEXPR(NDIM > 1) {
                            vg[0] = -omega * (X(1, i) + 0.5 * xloc[d][1] * dx);
                            vg[1] = +omega * (X(0, i) + 0.5 * xloc[d][0] * dx);
                            if
                                CONSTEXPR(NDIM == 3) {
							vg[2] = 0.0;
						}
                        }
                    else {
						vg[0] = 0.0;
					}
                // physics<NDIM>::flux(UL, UR, UL0, UR0, this_flux, dim, this_am, this_ap, vg, dx);
                physics<NDIM>::flux(UL, UR, this_flux, dim, this_am, this_ap, vg, dx);
					am = std::min(am, this_am);
					ap = std::max(ap, this_ap);
					for (int f = 0; f < nf; f++) {
                        fluxes(dim, f, i, fi) = this_flux[f];
					}
				}
            

            maxAmax = std::max(ap, safe_real(-am));

                // field update from fluxes
			for (int f = 0; f < nf; f++) {
                    F(dim, f, i) = 0.0;
					for (int fi = 0; fi < geo.NFACEDIR; fi++) {
                        const auto& w = weights[fi];
                        F(dim, f, i) += w * fluxes(dim, f, i, fi);
					}
				}
                // angular momentum update from linear momentum
			for (int angmom_pair = 0; angmom_pair < angmom_count; angmom_pair++) {
				const int sx_i = angmom_index + angmom_pair * (NDIM + geo.NANGMOM);
				const int zx_i = sx_i + NDIM;
				for (int n = 0; n < geo.NANGMOM; n++) {
                        F(dim, zx_i + n, i) = fluxes(dim, zx_i + n, i, 0);

					for (int m = 0; m < NDIM; m++) {
						if (dim != m) {
							for (int l = 0; l < NDIM; l++) {
								for (int fi = 0; fi < geo.NFACEDIR; fi++) {
									const auto d = faces[dim][fi];
                                        F(dim, zx_i + n, i) += weights[fi] * kdelta[n][m][l] *
                                            xloc[d][m] * 0.5 * dx * fluxes(dim, sx_i + l, i, fi);
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
}    // namespace octotiger