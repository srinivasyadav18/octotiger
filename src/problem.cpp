/*
 * problem.cpp
 *
 *  Created on: May 29, 2015
 *      Author: dmarce1
 */

#include "problem.hpp"
#include "grid.hpp"
#include "lane_emden.hpp"
#include <cmath>
#include "exact_sod.hpp"
#include "defs.hpp"

init_func_type problem = nullptr;
refine_test_type refine_test_function = refine_test;

bool refine_sod(integer level, integer max_level, real x, real y, real z, std::vector<real> U, std::array<std::vector<real>, NDIM> dudx) {
	for (integer i = 0; i != NDIM; ++i) {
		if (std::abs(dudx[i][rho_i] / U[rho_i]) > 0.1) {
			return level < max_level;
		}
	}
	return false;
}

bool refine_blast(integer level, integer max_level, real x, real y, real z, std::vector<real> U, std::array<std::vector<real>, NDIM> dudx) {
	for (integer i = 0; i != NDIM; ++i) {
		if (std::abs(dudx[i][rho_i] / U[rho_i]) > 0.1) {
			return level < max_level;
		}
		if (std::abs(dudx[i][tau_i]) > 0.01) {
			return level < max_level;
		}
	}
	return false;
}

bool refine_test(integer level, integer max_level, real x, real y, real z, std::vector<real> U, std::array<std::vector<real>, NDIM> dudx) {
	bool rc = false;
	real den_floor = 1.0e-4;
	integer test_level = max_level;
	for (integer this_test_level = test_level; this_test_level >= 1; --this_test_level) {
		if (U[rho_i] > den_floor) {
			rc = rc || (level < this_test_level);
		}
		if (rc) {
			break;
		}
		den_floor /= 8.0;
	}
	return rc;

}

bool refine_test_bibi(integer level, integer max_level, real x, real y, real z, std::vector<real> U, std::array<std::vector<real>, NDIM> dudx) {
	bool rc = false;
	real den_floor = 1.0e-4;
//	integer test_level = ((U[spc_de_i]+U[spc_dc_i]) < 0.5*U[rho_i] ? max_level  - 1 : max_level);
//	integer test_level = ((U[spc_ae_i]+U[spc_de_i]) > 0.5*U[rho_i] ? max_level  - 1 : max_level);
//	integer test_level = ((U[spc_dc_i]+U[spc_de_i]) > 0.5*U[rho_i] ? max_level  - 1 : max_level);
	integer test_level = max_level;
	for (integer this_test_level = test_level; this_test_level >= 1; --this_test_level) {
		if (U[rho_i] > den_floor) {
			rc = rc || (level < this_test_level);
		}
		if (rc) {
			break;
		}
		den_floor /= 8.0;
	}
	return rc;

}

void set_refine_test(const refine_test_type& rt) {
	refine_test_function = rt;
}

refine_test_type get_refine_test() {
	return refine_test_function;
}

void set_problem(const init_func_type& p) {
	problem = p;
}

init_func_type get_problem() {
	return problem;
}

std::vector<real> null_problem(real x, real y, real z, real dx) {
	std::vector < real > u(NF, real(0));
	return u;
}

std::vector<real> blast_wave(real x, real y, real z, real dx) {
	const real fgamma = grid::get_fgamma();
	x -= 0.453;
	y -= 0.043;
	std::vector < real > u(NF, real(0));
	u[spc_dc_i] = u[rho_i] = 1.0;
	const real a = std::sqrt(2.0) * dx;
	real r = std::sqrt(x * x + y * y + z * z);
	u[egas_i] = std::max(1.0e-10, exp(-r * r / a / a));
	u[tau_i] = std::pow(u[egas_i], ONE / fgamma);
	return u;
}

std::vector<real> sod_shock_tube(real x, real y, real z, real) {

	const real fgamma = grid::get_fgamma();
	std::vector < real > u(NF, real(0));
	if (x < real(0)) {
		u[spc_dc_i] = u[rho_i] = sod_init.rhol;
		u[egas_i] = sod_init.pl / (sod_init.gamma - 1.0);
	} else {
		u[spc_ac_i] = u[rho_i] = sod_init.rhor;
		u[egas_i] = sod_init.pr / (sod_init.gamma - 1.0);
	}

	u[tau_i] = std::pow(u[egas_i], ONE / fgamma);
	return u;
}

const real dxs = 0.0;
const real dys = -0.0;

std::vector<real> double_solid_sphere_analytic_phi(real x0, real y0, real z0) {
	std::vector < real > u(4, real(0));
	auto u1 = solid_sphere_analytic_phi(x0, y0, z0, dxs);
	auto u2 = solid_sphere_analytic_phi(x0, y0, z0, dys);
	for (integer f = 0; f != 4; ++f) {
		u[f] = u1[f] + u2[f];
	}
	return u;
}

const real ssr0 = 1.0 / 3.0;
std::vector<real> solid_sphere_analytic_phi(real x, real y, real z, real xshift) {
	const real r0 = ssr0;
	const real M = 1.0;
	std::vector < real > g(4);
	x -= xshift;
//	x0 -= -0.0444;
//	y0 -= +0.345;
//	z0 -= -.2565;
	const real r = std::sqrt(x * x + y * y + z * z);
	const real r3 = r * r * r;
	const real Menc = M * std::pow(std::min(r / r0, 1.0), 3);
	if (r < r0) {
		g[phi_i] = -M * (3.0 * r0 * r0 - r * r) / (2.0 * r0 * r0 * r0);
	} else {
		g[phi_i] = -M / r;
	}
	g[gx_i] = -Menc * x / r3;
	g[gy_i] = -Menc * y / r3;
	g[gz_i] = -Menc * z / r3;
	return g;
}

std::vector<real> double_solid_sphere(real x0, real y0, real z0, real dx) {
	std::vector < real > u(NF, real(0));
	auto u1 = solid_sphere(x0, y0, z0, dx, dxs);
	auto u2 = solid_sphere(x0, y0, z0, dx, dys);
	for (integer f = 0; f != NF; ++f) {
		u[f] = u1[f] + u2[f];
	}
	return u;
}

std::vector<real> solid_sphere(real x0, real y0, real z0, real dx, real xshift) {
	const integer N = 25;
	const real r0 = ssr0;
	const real rho_floor = 1.0e-50;
	const real V = 4.0 / 3.0 * M_PI * r0 * r0 * r0;
	const real drho = 1.0 / real(N * N * N) / V;
	std::vector < real > u(NF, real(0));
	x0 -= xshift;
//	x0 -= -0.0444;
//	y0 -= +0.345;
//	z0 -= -.2565;
	const auto mm = [](real a, real b) {
		if( a * b < ZERO ) {
			return ZERO;
		} else if( a > ZERO ) {
			return std::min(a,b);
		} else {
			return std::max(a,b);
		}
	};
	const real xmax = std::max(std::abs(x0 + dx / 2.0), std::abs(x0 - dx / 2.0));
	const real ymax = std::max(std::abs(y0 + dx / 2.0), std::abs(y0 - dx / 2.0));
	const real zmax = std::max(std::abs(z0 + dx / 2.0), std::abs(z0 - dx / 2.0));
	const real xmin = mm(x0 + dx / 2.0, x0 - dx / 2.0);
	const real ymin = mm(y0 + dx / 2.0, y0 - dx / 2.0);
	const real zmin = mm(z0 + dx / 2.0, z0 - dx / 2.0);
	if (xmax * xmax + ymax * ymax + zmax * zmax <= r0 * r0) {
		u[rho_i] += drho * N * N * N;
	} else if (xmin * xmin + ymin * ymin + zmin * zmin <= r0 * r0) {
		x0 -= dx / 2.0;
		y0 -= dx / 2.0;
		z0 -= dx / 2.0;
		const real d = dx / N;
		x0 += d / 2.0;
		y0 += d / 2.0;
		z0 += d / 2.0;
		const real r2 = r0 * r0;
		for (integer i = 0; i < N; ++i) {
			const real x2 = std::pow(x0 + real(i) * d, 2);
			for (integer j = 0; j < N; ++j) {
				const real y2 = std::pow(y0 + real(j) * d, 2);
#pragma GCC ivdep
				for (integer k = 0; k < N; ++k) {
					const real z2 = std::pow(z0 + real(k) * d, 2);
					if (x2 + y2 + z2 < r2) {
						u[rho_i] += drho;
					}
				}
			}
		}
	}
	u[rho_i] = std::max(u[rho_i], rho_floor);
	return u;
}

const real x0 = 0.0;
const real y0_ = 0.0;
const real z0 = 0.0;
const real rmax = 3.7;
const real dr = rmax / 128.0;
const real alpha = real(1) / real(5);

std::vector<real> star(real x, real y, real z, real) {
	const real fgamma = grid::get_fgamma();

	x -= 0.125;
	y -= 0.0;
	z -= 0.0;
	real menc;
	const real r = std::sqrt(x * x + y * y + z * z);
	std::vector < real > u(NF, real(0));
	real theta;
	const real n = real(1) / (fgamma - real(1));
	const real rho_min = 1.0e-3;
	const real theta_min = std::pow(rho_min, real(1) / n);
	const auto c0 = real(4) * real(M_PI) * alpha * alpha / (n + real(1));
	if (r <= rmax) {
		theta = lane_emden(r, dr);
		theta = std::max(theta, theta_min);
	} else {
		theta = theta_min;
	}
	u[rho_i] = std::pow(theta, n);
	u[egas_i] = std::pow(theta, fgamma * n) * c0 / (fgamma - real(1));
	if (theta <= theta_min) {
		u[egas_i] *= real(100);
	}
	u[tau_i] = std::pow(u[egas_i], (real(1) / real(fgamma)));
	u[sx_i] = -DEFAULT_OMEGA * y * u[rho_i];
	u[sy_i] = +DEFAULT_OMEGA * x * u[rho_i];
	return u;
}

std::vector<real> equal_mass_binary(real x, real y, real z, real) {
	const integer don_i = spc_ac_i;
	const integer acc_i = spc_dc_i;
	const real fgamma = grid::get_fgamma();

	real theta;
	real alpha = 1.0 / 15.0;
	const real n = real(1) / (fgamma - real(1));
	const real rho_min = 1.0e-12;
	std::vector < real > u(NF, real(0));
	const real d = 1.0 / 2.0;
	real x1 = x - d;
	real x2 = x + d;
	real y1 = y;
	real y2 = y;
	real z1 = z;
	real z2 = z;

	const real r1 = std::sqrt(x1 * x1 + y1 * y1 + z1 * z1) / alpha;
	const real r2 = std::sqrt(x2 * x2 + y2 * y2 + z2 * z2) / alpha;

	const real theta_min = std::pow(rho_min, real(1) / n);
	const auto c0 = real(4) * real(M_PI) * alpha * alpha / (n + real(1));

	if (r1 <= rmax || r2 <= rmax) {
		real r = std::min(r1, r2);
		theta = lane_emden(r, dr);
		theta = std::max(theta, theta_min);
	} else {
		theta = theta_min;
	}
	u[rho_i] = std::pow(theta, n);
	u[egas_i] = std::pow(theta, fgamma * n) * c0 / (fgamma - real(1));
	u[egas_i] = std::max(u[egas_i], ei_floor);
	u[tau_i] = std::pow(u[egas_i], (real(1) / real(fgamma)));
	u[sx_i] = -DEFAULT_OMEGA * y * u[rho_i];
	u[sy_i] = +DEFAULT_OMEGA * x * u[rho_i];
	u[egas_i] += HALF * DEFAULT_OMEGA * DEFAULT_OMEGA * (x * x + y * y) * u[rho_i];
	if (x < ZERO) {
		u[acc_i] = u[rho_i];
		u[don_i] = ZERO;
	} else {
		u[don_i] = u[rho_i];
		u[acc_i] = ZERO;
	}
	return u;
}
