//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef OCTOTIGER____FLUX____HPP123
#define OCTOTIGER____FLUX____HPP123

#include "octotiger/unitiger/physics.hpp"
#include "octotiger/unitiger/physics_impl.hpp"

template<int NDIM, int INX, class PHYS>
safe_real hydro_computer<NDIM, INX, PHYS>::flux(const hydro::state_type &U, const hydro::recon_type<NDIM> &Q, hydro::flux_type &F, hydro::x_type &X,
		safe_real omega) {

	PROFILE();

	static const cell_geometry<NDIM, INX> geo;

	static thread_local std::vector<std::vector<safe_real>> UR(geo.NFACEDIR, std::vector < safe_real > (nf_));
	static thread_local std::vector<std::vector<safe_real>> UL(geo.NFACEDIR, std::vector < safe_real > (nf_));
	static thread_local std::vector<std::vector<safe_real>> FR(geo.NFACEDIR, std::vector < safe_real > (nf_));
	static thread_local std::vector<std::vector<safe_real>> FL(geo.NFACEDIR, std::vector < safe_real > (nf_));
	static thread_local std::vector<safe_real> this_flux(nf_);

	static constexpr auto faces = geo.face_pts();
	static constexpr auto weights = geo.face_weight();
	static constexpr auto xloc = geo.xloc();
	static constexpr auto levi_civita = geo.levi_civita();

	const auto dx = X[0][geo.H_DNX] - X[0][0];

	safe_real amax = 0.0;
	for (int dim = 0; dim < NDIM; dim++) {

		const auto indices = geo.get_indexes(3, geo.face_pts()[dim][0]);

		// zero-initialize F
		for (int f = 0; f < nf_; f++) {
#pragma ivdep
			for (const auto &i : indices) {
				F[dim][f][i] = 0.0;
			}
		}

		for (const auto &i : indices) {
			safe_real ap = 0.0, am = 0.0;
			safe_real this_ap[geo.NFACEDIR], this_am[geo.NFACEDIR];
			for (int fi = 0; fi < geo.NFACEDIR; fi++) {
				const auto d = faces[dim][fi];
				for (int f = 0; f < nf_; f++) {
					UR[fi][f] = Q[f][d][i];
					UL[fi][f] = Q[f][geo::flip_dim(d, dim)][i - geo.H_DN[dim]];
				}
				std::array < safe_real, NDIM > x;
				std::array < safe_real, NDIM > vg;
				for (int dim = 0; dim < NDIM; dim++) {
					x[dim] = X[dim][i] + 0.5 * xloc[d][dim] * dx;
				}
				if constexpr (NDIM > 1) {
					vg[0] = -omega * (X[1][i] + 0.5 * xloc[d][1] * dx);
					vg[1] = +omega * (X[0][i] + 0.5 * xloc[d][0] * dx);
					if constexpr (NDIM == 3) {
						vg[2] = 0.0;
					}
				} else {
					vg[0] = 0.0;
				}

				safe_real amr, apr, aml, apl;

				PHYS::template physical_flux<INX>(UR[fi], FR[fi], dim, amr, apr, x, vg);
				PHYS::template physical_flux<INX>(UL[fi], FL[fi], dim, aml, apl, x, vg);
				this_ap[fi] = std::max(std::max(apr, apl), safe_real(0.0));
				this_am[fi] = std::min(std::min(amr, aml), safe_real(0.0));
			}
			double ap_max = 0.0, am_min = 0.0;
			if (experiment == 1) {
				for (int fi = 0; fi < geo.NFACEDIR; fi++) {
					ap_max = std::max(this_ap[fi], ap_max);
					am_min = std::min(this_am[fi], am_min);
				}
			}
			for (int fi = 0; fi < geo.NFACEDIR; fi++) {
				double ap, am;
				if (experiment == 1) {
					ap = ap_max;
					am = am_min;
				} else if (experiment == 2) {
					ap = this_ap[0];
					am = this_am[0];
				} else {
					ap = this_ap[fi];
					am = this_am[fi];
				}
#pragma ivdep
				for (int f = 0; f < nf_; f++) {
					if (this_ap - this_am != 0.0) {
						this_flux[f] = (ap * FL[fi][f] - am * FR[fi][f] + ap * am * (UR[fi][f] - UL[fi][f])) / (ap - am);
					} else {
						this_flux[f] = (FL[fi][f] + FR[fi][f]) / 2.0;
					}
				}
				am = std::min(am, this_am[fi]);
				ap = std::max(ap, this_ap[fi]);
#pragma ivdep
				for (int f = 0; f < nf_; f++) {
					// field update from flux
					F[dim][f][i] += weights[fi] * this_flux[f];
				}
			}
			const auto this_amax = std::max(ap, safe_real(-am));
			if (this_amax > amax) {
				amax = this_amax;
			}
		}
	}
	return amax;
}

#endif
