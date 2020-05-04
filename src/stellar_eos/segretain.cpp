#include "octotiger/stellar_eos/segretain_eos.hpp"
#include <algorithm>
#include <cmath>
;

void segretain_eos::set_units(double g, double cm, double s) {
	stellar_eos::set_units(g, cm, s);
	const auto c = 2.99792458e10 * s / cm;
	const auto me = 9.1093897e-28 / g;
	const auto h = 6.6260755e-27 * s / (g * cm * cm);
	B = 8.0 * M_PI * amu / 3.0 * std::pow(me * c / h, 3);
	A = M_PI * std::pow(me * c, 4) * c / 3 / std::pow(h, 3);
//	abort();
}

std::pair<double, double> segretain_eos::ztwd_pressure_and_energy(double rho, double abar, double zbar) {
	std::pair<double, double> rc;
	const auto mu = abar / zbar;
	const auto x = std::pow(rho / B / mu, 1.0 / 3.0);
	double Pdeg;
	double Edeg;
	if (x < 0.01) {
		const auto x2 = x * x;
		const auto x3 = x2 * x;
		const auto x5 = x3 * x2;
		Pdeg = 1.6 * A * x5;
		Edeg = 2.4 * A * x5;
	} else {
		Pdeg = A * (x * (2 * x * x - 3) * std::sqrt(x * x + 1) + 3 * asinh(x));
		const auto hdeg = 8 * A / (mu * B) * (std::sqrt(1 + x * x) - 1);
		Edeg = rho * hdeg - Pdeg;
	}
	rc.first = Pdeg;
	rc.second = Edeg;
	return rc;
}

double segretain_eos::pressure_from_energy(double rho, double e, double abar, double zbar) {
	const auto tmp = ztwd_pressure_and_energy(rho, abar, zbar);
	const double Pd = tmp.first;
	const double Ed = tmp.second;
//	e = std::max(e, Ed);
//	printf( "%e %e %e\n", Pdeg, Edeg, e);
	return sqrt(Pd * Pd + pow((fgamma - 1.0) * (e - Ed), 2));
}

std::pair<double, double> segretain_eos::pressure_and_soundspeed(double rho, double e, double abar, double zbar) {
	std::pair<double, double> rc;
	const auto mu = abar / zbar;
	const auto tmp = ztwd_pressure_and_energy(rho, abar, zbar);
	const double Pd = tmp.first;
	const double Ed = tmp.second;
//	e = std::max(e, Ed);
	const auto x = std::pow(rho / B / mu, 1.0 / 3.0);
	const auto dPddx = 8.0 * A * std::pow(x, 4) / std::sqrt(x * x + 1);
	const auto dEddx = x < 0.01 ? 12.0 * A * std::pow(x, 4) : 24.0 * A * x * x * (std::sqrt(x * x + 1) - 1);
	const auto dPddrho = dPddx / (3 * B * mu * x * x);
	const auto dEddrho = dEddx / (3 * B * mu * x * x);
	const auto P = sqrt(Pd * Pd + pow((fgamma - 1.0) * (e - Ed), 2));
	auto dPdrho = (dPddrho * Pd + pow(fgamma - 1, 2) * (e / rho - dEddrho) * (e - Ed)) / P;
	const auto dPdeps = pow(fgamma - 1, 2) * rho * (e - Ed) / P;
	const auto cs2 = dPdrho + P / (rho * rho) * dPdeps;
//	if (cs2 < 0.0) {
//		printf("%e %e %e %e %e %e %e\n", cs2, x, e, Ed, dPdrho, abar, zbar);
//	}
	rc.first = P;
	rc.second = std::sqrt(std::max(cs2, 0.0));
	return rc;
}

double segretain_eos::T_from_energy(double rho, double e, double abar, double zbar) {
	const auto tmp = ztwd_pressure_and_energy(rho, abar, zbar);
	const double Edeg = tmp.second;
	return (e - Edeg) * abar * amu / (rho * kb * (zbar + 1));
}

double segretain_eos::energy_from_T(double rho, double T, double abar, double zbar) {
	const auto tmp = ztwd_pressure_and_energy(rho, abar, zbar);
	const double Edeg = tmp.second;
	return Edeg + 1.0 / (fgamma - 1) * rho * T * kb * (zbar + 1) / (amu * abar);
}
