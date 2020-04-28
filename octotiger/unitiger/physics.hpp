//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef OCTOTIGER_UNITIGER_PHYSICS_HPP_
#define OCTOTIGER_UNITIGER_PHYSICS_HPP_

#include "octotiger/unitiger/safe_real.hpp"
#include "octotiger/test_problems/blast.hpp"
#include "octotiger/test_problems/exact_sod.hpp"
#include "octotiger/helmholtz.hpp"
#include <functional>

template<int NDIM>
struct physics {
	static constexpr char const *field_names3[] =
			{ "rho", "egas", "ein", "pot", "sx", "sy", "sz", "zx", "zy", "zz", "spc_1", "spc_2", "spc_3", "spc_4", "spc_5" };
	static constexpr char const *field_names2[] = { "rho", "egas", "ein", "pot", "sx", "sy", "zz", "spc_1", "spc_2", "spc_3", "spc_4", "spc_5" };
	static constexpr char const *field_names1[] = { "rho", "egas", "ein", "pot", "sx", "spc_1", "spc_2", "spc_3", "spc_4", "spc_5" };
	static constexpr int rho_i = 0;
	static constexpr int egas_i = 1;
	static constexpr int ein_i = 2;
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
		SOD, BLAST, KH, CONTACT, KEPLER
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

	static bool contact_field(int f) {
		return (f == rho_i || (f >= spc_i && f < spc_i + n_species_));
	}

	static void set_fgamma(safe_real fg);

	static void to_prim(std::vector<safe_real> u, safe_real &p, safe_real &v, safe_real &c, int dim);

	static void enforce_outflows(hydro::state_type &U, const hydro::x_type &X, int face) {

	}

	template<int INX>
	static void physical_flux(const std::vector<safe_real> &U, std::vector<safe_real> &F, int dim, safe_real &am, safe_real &ap, std::array<safe_real, NDIM> &x,
			std::array<safe_real, NDIM> &vg);

	template<int INX>
	static void post_process(hydro::state_type &U, const hydro::x_type &X, safe_real dx);

	template<int INX>
	static void source(hydro::state_type &dudt, const hydro::state_type &U, const hydro::flux_type &F, const hydro::x_type X, safe_real omega, safe_real dx);

	template<int INX>
	static void derivative_source(hydro::state_type &dudt, const hydro::state_type &U, const std::vector<std::vector<std::vector<safe_real>>> &Q,
			const hydro::x_type X, safe_real omega, safe_real dx);

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

	static void set_dual_energy_switches(safe_real one, safe_real two);

	static void set_central_force(safe_real GM) {
		GM_ = GM;
	}
	static int get_angmom_index() {
		return sx_i;
	}

	template<int INX>
	static void enforce_outflow(hydro::state_type &U, int dim, int dir);

	using eos_func = std::function<double(double,double,double, double)>;
	using eos_func2 = std::function<std::pair<double,double>(double,double,double, double)>;

	static void set_code_units(double g, double cm, double s) {
		code_to_g = g;
		code_to_cm = cm;
		code_to_s = s;
		mh = 1.6605402e-24 / g;
		//		mh = 1.6733e-24 / g;
		kb = 1.380658e-16 * s * s / (g * cm * cm);
		const auto c = 2.99792458e10 * s / cm;
		const auto me = 9.1093897e-28 / g;
		const auto h = 6.6260755e-27 * s / (g * cm * cm);
		B_ = 8.0 * M_PI * mh / 3.0 * std::pow(me * c / h, 3);
		A_ = M_PI * std::pow(me * c, 4) * c / 3 / std::pow(h, 3);
	}

	static void set_helmholtz_eos() {
		pressure_from_energy = helmholtz_pressure_from_energy;
		pressure_and_soundspeed = helmholtz_pressure_and_soundspeed;
		T_from_energy = helmholtz_T_from_energy;
		helmholtz_set_cgs_units(code_to_cm, code_to_g, code_to_s, 1);

	}

	static void set_segretain_eos() {
		pressure_from_energy = segretain_pressure_from_energy;
		pressure_and_soundspeed = segretain_pressure_and_soundspeed;
		T_from_energy = segretain_T_from_energy;
	}

	static eos_func pressure_from_energy;
	static eos_func T_from_energy;
	static eos_func2 pressure_and_soundspeed;

	static void set_atomic_data(const std::vector<double> a, const std::vector<double> z) {
		for (int s = 0; s < n_species_; s++) {
			A[s] = a[s];
			Z[s] = z[s];
		}
	}
	static std::pair<double,double> ztwd_pressure_and_energy(double rho,  double abar, double zbar) {
		std::pair<double,double> rc;
		const auto mu = abar / zbar;
		const auto x = std::pow(rho / B_ / mu, 1.0 / 3.0);
		double Pdeg;
		double Edeg;
		if (x < 0.01) {
			Pdeg = 1.6 * A_ * x * x * x * x * x;
			Edeg = 2.4 * A_ * x * x * x * x * x;
		} else {
			Pdeg = A_ * (x * (2 * x * x - 3) * std::sqrt(x * x + 1) + 3 * asinh(x));
			const auto hdeg = 8 * A_ / (mu*B_) * (std::sqrt(1 + x * x) - 1);
			Edeg = rho * hdeg - Pdeg;
		}
		rc.first = Pdeg;
		rc.second = Edeg;
		return rc;
	}

private:
	static safe_real rho_sink_radius_;
	static safe_real rho_sink_floor_;
	static int nf_;
	constexpr static int n_species_ = 5;
	static safe_real fgamma_;
	static safe_real GM_;
	static std::array<safe_real, n_species_> A;
	static std::array<safe_real, n_species_> Z;

	static double code_to_g, code_to_cm, code_to_s, mh, kb, A_, B_;

	static double ideal_pressure_from_energy(double rho, double e, double, double) {
		return (fgamma_ - 1.0) * e;
	}

	static std::pair<double, double> ideal_pressure_and_soundspeed(double rho, double e, double abar, double zbar) {
		std::pair<double, double> pcs;
		pcs.first = ideal_pressure_from_energy(rho, e, abar, zbar);
		pcs.second = std::sqrt(std::max(fgamma_ * pcs.first / rho, 0.0));
		return pcs;
	}

	static double ideal_T_from_energy(double rho, double e, double abar, double zbar) {
		return (fgamma_ - 1.0) * e * abar * mh / (rho * kb * (zbar + 1));
	}


	static double segretain_pressure_from_energy(double rho, double e, double abar, double zbar) {
		const auto tmp = ztwd_pressure_and_energy(rho, abar, zbar);
		const double Pdeg = tmp.first;
		const double Edeg = tmp.second;
		return Pdeg + (fgamma_ - 1.0) * (e - Edeg);
	}

	static std::pair<double, double> segretain_pressure_and_soundspeed(double rho, double e, double abar, double zbar) {
		std::pair<double, double> rc;
		const auto mu = abar / zbar;
		const auto tmp = ztwd_pressure_and_energy(rho, abar, zbar);
		const double Pdeg = tmp.first;
		const double Edeg = tmp.second;
		const auto x = std::pow(rho / B_ / mu, 1.0 / 3.0);
		const auto dPdeg_drho = 8 * A_ * x * x / B_ / mu / std::sqrt(x * x + 1);
		const auto cs2 = dPdeg_drho + fgamma_ * (fgamma_ - 1) * (e - Edeg) / rho;
		rc.first = Pdeg + (fgamma_-1)*(e-Edeg);
		rc.second = std::sqrt(std::max(cs2, 0.0));
		return rc;
	}

	static double segretain_T_from_energy(double rho, double e, double abar, double zbar) {
		const auto tmp = ztwd_pressure_and_energy(rho, abar, zbar);
		const double Edeg = tmp.second;
		return (e - Edeg) * abar * mh / (rho * kb * (zbar + 1));
	}

};

template<int NDIM>
std::array<safe_real, physics<NDIM>::n_species_> physics<NDIM>::A = { 1, 1, 1, 1, 1 };

template<int NDIM>
std::array<safe_real, physics<NDIM>::n_species_> physics<NDIM>::Z = { 1, 1, 1, 1, 1 };

template<int NDIM>
typename physics<NDIM>::eos_func physics<NDIM>::pressure_from_energy = ideal_pressure_from_energy;

template<int NDIM>
typename physics<NDIM>::eos_func physics<NDIM>::T_from_energy = ideal_T_from_energy;

template<int NDIM>
typename physics<NDIM>::eos_func2 physics<NDIM>::pressure_and_soundspeed = ideal_pressure_and_soundspeed;

template<int NDIM>
double physics<NDIM>::mh = 1.0;

template<int NDIM>
double physics<NDIM>::A_ = 1.0;

template<int NDIM>
double physics<NDIM>:: B_ = 1.0;

template<int NDIM>
double physics<NDIM>::kb = 1.0;



template<int NDIM>
double physics<NDIM>::code_to_g;

template<int NDIM>
double physics<NDIM>::code_to_cm;

template<int NDIM>
double physics<NDIM>::code_to_s;

//definitions of the declarations (and initializations) of the static constexpr variables
template<int NDIM>
constexpr char const *physics<NDIM>::field_names1[];
template<int NDIM>
constexpr char const *physics<NDIM>::field_names2[];
template<int NDIM>
constexpr char const *physics<NDIM>::field_names3[];

template<int NDIM>
safe_real physics<NDIM>::rho_sink_radius_ = 0.0;

template<int NDIM>
safe_real physics<NDIM>::rho_sink_floor_ = 0.0;

template<int NDIM>
safe_real physics<NDIM>::GM_ = 0.0;

template<int NDIM>
safe_real physics<NDIM>::de_switch_1 = 1e-3;

template<int NDIM>
safe_real physics<NDIM>::de_switch_2 = 1e-1;

template<int NDIM>
int physics<NDIM>::nf_ = (4 + NDIM + (NDIM == 1 ? 0 : std::pow(3, NDIM - 2))) + physics<NDIM>::n_species_;

template<int NDIM>
safe_real physics<NDIM>::fgamma_ = 7. / 5.;

#endif /* OCTOTIGER_UNITIGER_PHYSICS_HPP_ */
