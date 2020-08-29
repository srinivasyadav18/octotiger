//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

#include "octotiger/config/export_definitions.hpp"
#include "octotiger/defs.hpp"
#include "octotiger/interaction_types.hpp"
#include "octotiger/options_enum.hpp"
#include "octotiger/real.hpp"

#include <hpx/include/naming.hpp>

#include <cstddef>
#include <string>
#include <vector>

/* Must look like this - no spaces
 COMMAND_LINE_ENUM(problem_type,DWD,SOD,BLAST,NONE,SOLID_SPHERE,STAR,MOVING_STAR,RADIATION_TEST,ROTATING_STAR,MARSHAK,AMR_TEST);

 COMMAND_LINE_ENUM(eos_type,IDEAL,WD);
 */

COMMAND_LINE_ENUM(problem_type, DWD, SOD, BLAST, NONE, SOLID_SPHERE, STAR, MOVING_STAR, RADIATION_TEST, ROTATING_STAR, MARSHAK, AMR_TEST, ADVECTION);

COMMAND_LINE_ENUM(eos_type, IDEAL, WD);

class options {
public:

	//! Turns inflow boundary condition on. Default = off
	//! This is a direct copy from the outermost interior cell to all the exterior cells, allowing inflow.
	bool inflow_bc;

	//! Turns on relfecting boundary conditions. Default = off. Does not work with gravity
	bool reflect_bc;

	//! Do not use.
	int experiment;

	//! Turns on contact discontinuity detection. Default = on
	bool cdisc_detect;

	//! Forces refinement of all nodes to max_level, effecting a unigrid. default = off
	bool unigrid;

	//! Turns off the binary diagnostics. default = off
	bool disable_diagnostics;

	// ! Enables benchmark test. Deprecated. Default = off
	bool bench;

	//! Disabled SILO output. Default = off
	bool disable_output;

	//! Adds an extra level of refinement to a donor or accretor core
	bool core_refine;

	//! Enables gravity. Default = on
	bool gravity;

	//! Enables hydrodynamics. Default = on
	bool hydro;

	//! Enables radiation hydrodynamics. Default = off.
	bool radiation;

	//! Sets the reduced speed of light factor. Default = 1.0. Valid ranges > 0 && <= 1.0
	real clight_reduce;

	//! Enables v1309 SCF. Default = off
	bool v1309;

	//! Enables the implicit part of the radiation step. Default = on
	bool rad_implicit;

	//! Enables the immediate re-writing of a checkpoint file on startup. Default = off
	bool rewrite_silo;

	//! Deprecated. DO NOT USE
	bool correct_am_grav;

	//! Enables the gravitational angular momentum correction. Default = on
	bool correct_am_hydro;

	//! For the rotating star problem this creates an AMR boundary through the middle of the star. Default = off
	bool rotating_star_amr;

	//! Enables the output of idle rates in the SILO files Default = on
	bool idle_rates;

	//! Frequency to output SILO when running SCF. Default = 25
	integer scf_output_frequency;

	//! Number of parallel files for SILO output. Must not be greater than the number of processors. Default = number of HPX localities
	integer silo_num_groups;

	//! Deprecated
	integer amrbnd_order;

	//! Number of additional regrid iterations to perform on startup. Can be used to increase the resolution of a checkpoint by
	//! refinement. Default = 0.
	integer extra_regrid;

	//! Number of extra levels of refinement for the accretor. Default = 0
	integer accretor_refine;

	//! Number of extra levels of refinement for the donor. Default = 0.
	integer donor_refine;

	//! Minimum number of refinement levels. Default = 1
	integer min_level;

	//! Maximum number of refinement levels. Default = 1
	integer max_level;

	//! For the binary (dwd) problem, target this number of total subgrids by adjusting density floor. Set to -1 to dissable. Default = -1
	integer ngrids;

	//! Stop the evolution after stop_step number of steps. Set to -1 to disable. Default = -1.
	integer stop_step;

	//! Number of x subgrid cells to offset SILO output by. Can be used to expose AMR cells. Default = 0;
	integer silo_offset_x;

	//! Number of y subgrid cells to offset SILO output by. Can be used to expose AMR cells. Default = 0;

	integer silo_offset_y;

	//! Number of z subgrid cells to offset SILO output by. Can be used to expose AMR cells. Default = 0;
	integer silo_offset_z;

	//! Deprecated
	integer future_wait_time;

	//! Maximum allowed percentage change for rho and tau per timestep. Default = 0.33333333333
	real dt_max;

	//! Energy for blast wave problem. Default = 1.0
	real eblast0;

	//! x location to place center of star for rotating star problem. Default = 0.0
	real rotating_star_x;

	//! Dual energy update switch. Default = 0.1. Set to 1.0 for polytropic EOS.
	real dual_energy_sw2;

	//! Dual energy presure switch. Default = 0.001. Set to 1.0 for polytropic EOS.
	real dual_energy_sw1;

	//! Deprecated
	real hard_dt;

	//! Driving rate in relative amount of angular momentum to remove every orbit for the binary (dwd) problem. Default = 0.0
	real driving_rate;

	//! Time to drive angular momentum in binary (dwd) problem. Default = 0.0;
	real driving_time;

	//! Deprecated
	real entropy_driving_rate;

	//! Deprecated
	real entropy_driving_time;

	//! Fixed rotoational rate of the grid. Default = 0.0;
	real omega;

	//! Output frequency. For the binary (dwd) problem this is in units of initial orbital periods. For all others it is in absolute time. Default = 0.01
	real output_dt;

	//! refinement density floor in code units. Used tor density based refinment criteria. Default = 1.0e-3
	real refinement_floor;

	//! Stop the simulation when it reaches this time in code units. Set to -1 to diable. Default =
	real stop_time;

	//! Opening criteria for FMM solver. For 8x8x8 subgrids valid values are between 0.34 and 0.5. Default=0.5
	real theta;

	//! The computational domain runs from [-xscale:xscale] in each dimension, in code units. Default = 1.0
	real xscale;

	//! TODO: Dominic document these after cleaning them up
	real code_to_g;
	real code_to_s;
	real code_to_cm;

	//! CFL factor. Default = 0.4
	real cfl;

	//! density floor in code units. Set to < 0.0 to disable. Default = -1
	real rho_floor;

	//! tau floor in code units. Set to < 0.0 to disable. Default = -1.
	real tau_floor;

	// right density for Sod problem Default = 1.0
	real sod_rhol;
	// left density for Sod problem Default = 0.125
	real sod_rhor;
	// left pressure for Sod problem Default = 1.0
	real sod_pl;
	// right pressure for Sod = Default = 0.1
	real sod_pr;

	// TODO: Sagiv please document these:
	real sod_theta;
	real sod_phi;
	real sod_gamma;

	real solid_sphere_xcenter;
	real solid_sphere_ycenter;
	real solid_sphere_zcenter;
	real solid_sphere_radius;
	real solid_sphere_mass;
	real solid_sphere_rho_min;

	real star_xcenter;
	real star_ycenter;
	real star_zcenter;
	real star_rmax;
	real star_alpha;
	real star_rho_out;
	real star_dr;
	real star_n;
	real star_rho_center;

	real moving_star_xvelocity;
	real moving_star_yvelocity;
	real moving_star_zvelocity;

	// TODO: Gregor please document these
	size_t cuda_streams_per_locality;
	size_t cuda_streams_per_gpu;
	size_t cuda_scheduling_threads;

	//! Input IC file for rotating star problem. Default = empty
	std::string input_file;

	//! Config file with commandline parameters. Default = empty
	std::string config_file;

	//! Directory for SILO and text output
	std::string data_dir;

	//! When this string is nonempty, Octo-Tiger will start-up but then immediately write to output_filename before the first time-step and exit. Default is empty
	std::string output_filename;

	//! Restart checkpoint SILO filename. Default is empty. When empty Octo-Tiger will initialized based on the problem selected
	std::string restart_filename;

	//! Number of species. Must be at least 1. Must be at least 5 for the binary (dwd) problem. Default = 5
	integer n_species;

	//! Number of fiels. Do not change this, set by Octo-Tiger.
	integer n_fields;

	//! Equation of state type. "ideal" for ideal gas, "wd" for cold wdeos + ideal gas. Default=ideal
	eos_type eos;

	//! Problem selection.

	//! Supported choices for problem are
	//! DWD - construct a binary using the SCF and run it
	//! SOD - Sod shock tube
	//! BLAST -Sedov-Taylor blast wave
	//! SOLID_SPHERE - Gravity of a solid sphere
	//! MOVING_STAR - A moving polytrope
	//! ROTATING_STAR - A rotating polytrope
	//! MARSHAK - Radiation Marshak test
	problem_type problem;

	// TODO: Gregor
	interaction_kernel_type m2m_kernel_type;
	interaction_kernel_type p2m_kernel_type;
	interaction_kernel_type p2p_kernel_type;

	// Specify atomic mass of each species. Used for temperature and opacity calculations. Default for each = solar
	std::vector<real> atomic_mass;

	// Specify atomic number of each species. Used for temperature and opacity calculations. Default for each = solar
	std::vector<real> atomic_number;

	// Specify hydrogen fraction for each species. Used for opacity calculations. Default for each = solar
	std::vector<real> X;

	// Specify metallicity fraction for each species. Used for opacity calculations. Default for each = solar
	std::vector<real> Z;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & eblast0;
		arc & rho_floor;
		arc & tau_floor;
		arc & sod_rhol;
		arc & sod_rhor;
		arc & sod_pl;
		arc & sod_pr;
		arc & sod_theta;
		arc & sod_phi;
		arc & sod_gamma;
		arc & solid_sphere_xcenter;
		arc & solid_sphere_ycenter;
		arc & solid_sphere_zcenter;
		arc & solid_sphere_radius;
		arc & solid_sphere_mass;
		arc & solid_sphere_rho_min;
		arc & star_xcenter;
		arc & star_ycenter;
		arc & star_zcenter;
		arc & star_rmax;
		arc & star_alpha;
		arc & star_dr;
		arc & star_n;
		arc & star_rho_center;
		arc & star_rho_out;
		arc & moving_star_xvelocity;
		arc & moving_star_yvelocity;
		arc & moving_star_zvelocity;
		arc & inflow_bc;
		arc & reflect_bc;
		arc & cdisc_detect;
		arc & experiment;
		arc & unigrid;
		arc & rotating_star_amr;
		arc & rotating_star_x;
		arc & future_wait_time;
		arc & silo_offset_x;
		arc & silo_offset_y;
		arc & silo_offset_z;
		arc & scf_output_frequency;
		arc & silo_num_groups;
		arc & amrbnd_order;
		arc & dual_energy_sw1;
		arc & dual_energy_sw2;
		arc & hard_dt;
		arc & correct_am_grav;
		arc & correct_am_hydro;
		arc & rewrite_silo;
		arc & rad_implicit;
		arc & n_fields;
		arc & n_species;
		arc & input_file;
		arc & config_file;
		arc & hydro;
		arc & gravity;
		arc & bench;
		arc & radiation;
		arc & m2m_kernel_type;
		arc & p2m_kernel_type;
		arc & p2p_kernel_type;
		arc & entropy_driving_rate;
		arc & entropy_driving_time;
		arc & driving_rate;
		arc & driving_time;
		arc & refinement_floor;
		arc & ngrids;
		arc & v1309;
		arc & clight_reduce;
		arc & stop_time;
		arc & min_level;
		arc & max_level;
		arc & xscale;
		arc & dt_max;
		arc & cfl;
		arc & omega;
		arc & restart_filename;
		arc & output_filename;
		arc & output_dt;
		arc & stop_step;
		arc & disable_diagnostics;
		arc & disable_output;
		arc & theta;
		arc & core_refine;
		arc & donor_refine;
		arc & extra_regrid;
		arc & accretor_refine;
		arc & idle_rates;
		int tmp = problem;
		arc & tmp;
		problem = static_cast<problem_type>(tmp);
		tmp = eos;
		arc & tmp;
		eos = static_cast<eos_type>(tmp);
		arc & data_dir;
		arc & m2m_kernel_type;
		arc & p2p_kernel_type;
		arc & p2m_kernel_type;
		arc & cuda_streams_per_locality;
		arc & cuda_streams_per_gpu;
		arc & cuda_scheduling_threads;
		arc & atomic_mass;
		arc & atomic_number;
		arc & X;
		arc & Z;
		arc & code_to_g;
		arc & code_to_s;
		arc & code_to_cm;
	}

	OCTOTIGER_EXPORT bool process_options(int argc, char *argv[]);

	static OCTOTIGER_EXPORT std::vector<hpx::id_type> all_localities;
};

OCTOTIGER_EXPORT options& opts();

template<class T = real>
struct hydro_state_t: public std::vector<T> {
	hydro_state_t() :
			std::vector<T>(opts().n_fields) {
	}
};

#endif /* OPTIONS_HPP_ */
