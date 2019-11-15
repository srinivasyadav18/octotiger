//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef OCTOTIGER_UNITIGER_PHYSICS_HPP_
#define OCTOTIGER_UNITIGER_PHYSICS_HPP_

#include "octotiger/unitiger/safe_real.hpp"
#include "octotiger/test_problems/blast.hpp"
#include "octotiger/test_problems/exact_sod.hpp"

template<int NDIM>
struct physics {

	static constexpr int rho_i = 0;
	static constexpr int egas_i = 1;
	static constexpr int tau_i = 2;
	static constexpr int pot_i = 3;
	static constexpr int sx_i = 4;
	static constexpr int sy_i = 5;
	static constexpr int sz_i = 6;
	static constexpr int lx_i = 4 + NDIM;
	static constexpr int ly_i = 5 + NDIM;
	static constexpr int lz_i = 6 + NDIM;
	static constexpr int spc_i = 4 + NDIM + (NDIM == 1 ? 0 : std::pow(3, NDIM - 2));
	static safe_real de_switch_1;
	static safe_real de_switch_2;

	enum test_type {
		SOD, BLAST, KH, CONTACT
	};

	static std::string get_test_type_string(test_type t) {
		switch (t) {
		case SOD:
			return "SOD";
		case BLAST:
			return "BLAST";
		case KH:
			return "KH";
		case CONTACT:
			return "CONTACT";
		default:
			return "OCTOTIGER";
		}
	}

	static int field_count();

	static void set_fgamma(safe_real fg);

	static void to_prim(std::vector<safe_real> u, safe_real &p, safe_real &v, int dim);

	static void physical_flux(const std::vector<safe_real> &U, std::vector<safe_real> &F, int dim, safe_real &am, safe_real &ap, std::array<safe_real, NDIM> &x,
			std::array<safe_real, NDIM> &vg);

	template<int INX>
	static void post_process(hydro::state_type &U, safe_real dx);

	template<int INX>
	static void source(hydro::state_type &dudt, const hydro::state_type &U, const hydro::flux_type &F, const hydro::x_type X, safe_real omega, safe_real dx);

	/*** Reconstruct uses this - GPUize****/
	template<int INX>
	static const hydro::state_type& pre_recon(const hydro::state_type &U, const hydro::x_type X, safe_real omega, bool angmom);
	/*** Reconstruct uses this - GPUize****/
	template<int INX>
	static void post_recon(std::vector<std::vector<std::vector<safe_real>>> &Q, const hydro::x_type X, safe_real omega, bool angmom);
	template<int INX>
	using comp_type = hydro_computer<NDIM, INX, physics<NDIM>>;

	template<int INX>
	std::vector<typename comp_type<INX>::bc_type> initialize(test_type t, hydro::state_type &U, hydro::x_type &X);

	template<int INX>
	static void analytic_solution(test_type test, hydro::state_type &U, const hydro::x_type &X, safe_real time);

	template<int INX>
	static const std::vector<std::vector<double>>& find_contact_discs(const hydro::state_type &U);

	static void set_n_species(int n);

	static void set_dual_energy_switches(safe_real one, safe_real two) {
		de_switch_1 = one;
		de_switch_2 = two;
	}

	static int get_angmom_index() {
		return sx_i;
	}

	inline safe_real x_deg(safe_real rho) {
		return std::pow(rho / B_, 1.0 / 3.0);
	}

	inline safe_real H_deg(safe_real x) {
		return 8.0 * A_ / B_ * std::sqrt(x * x + 1.0);
	}

	inline safe_real P_deg(safe_real h, safe_real x) {
		if (x < 0.01) {
			return 1.6 * A_ * std::pow(x, 5);
		} else {
			return B_ * h / 8.0 * x * (2 * x * x - 3) + 3 * A_ * asinh(x);
		}
	}

	inline safe_real E_deg(safe_real h, safe_real p, safe_real x) {
		if (x < 0.01) {
			return 2.4 * A_ * std::pow(x, 5);
		} else {
			return B_ * std::pow(x, 3) * h - p;
		}

	}

	inline safe_real dP_deg_drho(safe_real h, safe_real x) {
		return (64.0 / 3.0) * std::pow(A_ / B_, 2) * x * x / h;
	}

	inline safe_real thermal_energy(safe_real rho, safe_real egas, safe_real tau, safe_real ek) {
		const auto x = x_deg(rho);
		const auto h_deg = H_deg(x);
		const auto p_deg = P_deg(h_deg, x);
		const auto e_deg = E_deg(h_deg, p_deg, x);
		auto etherm = egas - e_deg - ek;
		if (etherm < de_switch_1 * egas) {
			return std::pow(tau, fgamma_);
		} else {
			return etherm;
		}
	}

	inline safe_real pressure(safe_real rho, safe_real egas, safe_real tau, safe_real ek) {
		const auto x = x_deg(rho);
		const auto h_deg = H_deg(x);
		const auto p_deg = P_deg(h_deg, x);
		const auto e_deg = E_deg(h_deg, p_deg, x);
		auto etherm = egas - e_deg - ek;
		if (etherm < de_switch_1 * egas) {
			etherm = std::pow(tau, fgamma_);
		}
		return p_deg + (fgamma_ - 1.0) * etherm;
	}

	inline std::pair<safe_real, safe_real> pressure_and_sound_speed(safe_real rho, safe_real egas, safe_real tau, safe_real ek) {
		const auto x = x_deg(rho);
		const auto h_deg = H_deg(x);
		const auto p_deg = P_deg(h_deg, x);
		const auto e_deg = E_deg(h_deg, p_deg, x);
		auto etherm = egas - e_deg - ek;
		if (etherm < de_switch_1 * egas) {
			etherm = std::pow(tau, fgamma_);
		}
		const auto p = p_deg + (fgamma_ - 1.0) * etherm;
		const auto dp_deps = (fgamma_ - 1.0) * rho;
		const auto dp_drho = dPdeg_drho(h_deg, x) + (fgamma_ - 1.0) * etherm / rho;
		auto c = dp_deps * p / (rho * rho) + dp_drho;
		return std::make_pair(p, std::sqrt(std::max(0.0, c)));
	}

private:
	static safe_real A_;
	static safe_real B_;
	static int nf_;
	static int n_species_;
	static safe_real fgamma_;

};

template<int NDIM>
safe_real physics<NDIM>::A_ = 0.0;

template<int NDIM>
safe_real physics<NDIM>::B_ = 1.0;

template<int NDIM>
safe_real physics<NDIM>::de_switch_1 = 1e-3;

template<int NDIM>
safe_real physics<NDIM>::de_switch_2 = 1e-1;

template<int NDIM>
int physics<NDIM>::nf_ = (4 + NDIM + (NDIM == 1 ? 0 : std::pow(3, NDIM - 2))) + physics<NDIM>::n_species_;

template<int NDIM>
int physics<NDIM>::n_species_ = 5;

template<int NDIM>
safe_real physics<NDIM>::fgamma_ = 7. / 5.;

#endif /* OCTOTIGER_UNITIGER_PHYSICS_HPP_ */
