#include "octotiger/grid.hpp"

#include <hpx/runtime/threads/run_as_os_thread.hpp>

#include "octotiger/test_problems/blast.hpp"

#include <algorithm>
#include <functional>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>

#if !defined(OCTOTIGER_HAVE_BOOST_MULTIPRECISION)
#include <quadmath.h>
using sed_real = __float128;
#else
#include <boost/multiprecision/cpp_bin_float.hpp>
using sed_real = boost::multiprecision::cpp_bin_float_quad;
#endif




/*extern "C" {*/
/* Subroutine */int sed_1d__(sed_real *time, int *nstep,
		sed_real * xpos, sed_real *eblast, sed_real *omega_in__,
		sed_real * xgeom_in__, sed_real *rho0, sed_real *vel0,
		sed_real *ener0, sed_real *pres0, sed_real *cs0, sed_real *gam0,
		sed_real *den, sed_real *ener, sed_real *pres, sed_real *vel,
		sed_real *cs);
//}

constexpr real blast_wave_t0 = 7e-4;


std::vector<real> blast_wave_analytic(real x, real y, real z, real t) {
	real r = std::sqrt(x * x + y * y + z * z);
	t += blast_wave_t0;
	real rmax = 3.0 * opts().xscale;
	real d, v, p;
	sedov::solution(t, r, rmax, d, v, p);
	std::vector<real> u(opts().n_fields, 0.0);
	u[rho_i] = u[spc_i] = std::max(d,1.0e-20);
	real s = d * v;
	u[sx_i] = s * x / r;
	u[sy_i] = s * y / r;
	u[sz_i] = s * z / r;
	real e = std::max(p / (grid::get_fgamma() - 1),1.0e-20);
	u[egas_i] = e + s * v * 0.5;
	u[tau_i] = std::pow(e, 1 / grid::get_fgamma());
	return u;
}

std::vector<real> blast_wave(real x, real y, real z, real dx) {
	return blast_wave_analytic(x,y,z,0.0);


}
