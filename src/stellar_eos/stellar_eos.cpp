/*
 * stellar_eos.cpp
 *
 *  Created on: May 3, 2020
 *      Author: dmarce1
 */

#include "octotiger/stellar_eos/stellar_eos.hpp"

stellar_eos::stellar_eos() {
	fgamma = 5.0 / 3.0;
	g_to_code = s_to_code = cm_to_code = 1.0;
	amu = 1.66053878283e-24;
	kb = 1.380650424e-16;
}

void stellar_eos::set_units(double g, double cm, double s) {
	g_to_code = g;
	cm_to_code = cm;
	s_to_code = s;
	amu = 1.66053878283e-24 / g_to_code;
	kb = 1.380650424e-16 / g_to_code / cm_to_code / cm_to_code * s_to_code * s_to_code;
}

void stellar_eos::set_fgamma(double fg) {
	fgamma = fg;
}
