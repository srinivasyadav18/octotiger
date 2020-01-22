//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <fenv.h>
#include <time.h>
#include <type_traits>
#include <typeinfo>

#define OCTOTIGER_GRIDDIM 8    // usually set in build scripts

#include <Kokkos_Core.hpp>

#include <hpx/hpx_main.hpp>

#include "octotiger/geometry.hpp"
#include "octotiger/unitiger/hydro.hpp"
#include "octotiger/unitiger/physics.hpp"

#include "flux_kokkos.hpp"
#include "octotiger/unitiger/safe_real.hpp"

static constexpr double tmax = 2.49e-5;
static constexpr safe_real dt_out = tmax / 249;

// #define H_BW 3
// #define H_NX (INX + 2 * H_BW)
// #define H_N3 std::pow(INX+2*H_BW,NDIM)

// why is there two different INX?
// one is constexpr in hydro_defs.hpp - used with unitiger
// one the template parameter used here 
// ask Dominic again - how are they different?


hydro::recon_type<2> first_q;

template <typename Viewtype>
void printContents(Viewtype& printView){
  Kokkos::parallel_for("print some contents",
      Kokkos::MDRangePolicy<typename Viewtype::execution_space, Kokkos::Rank<2>>(
          {0, 0}, {printView.extent(0), printView.extent(1)}),
      KOKKOS_LAMBDA(int j, int k) {
        // if((j*k)%1000 == 1){
          printf("%d,%d, %f; ", j, k, printView(j, k));
        // }
      });
}

template<int NDIM>
void printState(
	std::vector<std::vector<safe_real>> U,
	std::vector<std::vector<safe_real>> U0,
	std::vector<std::vector<std::vector<safe_real>>> F,
	hydro::x_type X,
	const safe_real omega,
	hydro::recon_type<NDIM> q,
	const safe_real a){
		for (auto & u0 : U){
			for (auto & u : u0){
				printf("%e; ", u);
			}
		}
     	printf("\n");
     	printf("\n");
		for (auto & u0 : U0){
			for (auto & u : u0){
				printf("%e; ", u);
			}
		}
     	printf("\n");
     	printf("\n");
		for (auto & f0 : F){
			for (auto & f : f0){
				for (auto & u : f){
					printf("%e; ", u);
				}
		     	printf("\n");
			}
	     	printf("\n");
		}
     	printf("\n");
     	printf("\n");
		for (auto & x : X){
			for (auto & u : x){
				printf("%e; ", u);
			}
		}
     	printf("\n");
     	printf("\n");
     	printf("%f; \n", omega);
		for (auto & q0 : q){
			for (auto & x : q0){
				for (auto & u : x){
					printf("%e; ", u);
				}
			}
		}
     	printf("\n");
     	printf("\n");
     	printf("\n");
}


template<int NDIM, int INX, bool with_kokkos>
void run_test(typename physics<NDIM>::test_type problem, bool with_correction) {
	static constexpr int H_BW = 3;
	// static constexpr auto H_NX = INX + 2 * H_BW;
	// static constexpr int H_N3 = std::pow(INX+2*H_BW,NDIM);
    static constexpr auto H_N3 = static_cast<int>(PowerNDIM(INX+2*H_BW));
	printf("test H_N3 %d, H_BW %d, NDIM %d, INX %d \n", H_N3, H_BW, NDIM, INX);
	// static_assert(std::pow(206, NDIM) == PowerNDIM(206));
	// static_assert(H_N3 == PowerNDIM(INX+2*H_BW));
	static constexpr safe_real CFL = (0.4 / NDIM);
	hydro_computer<NDIM, INX> computer;
	if (with_correction) {
		computer.use_angmom_correction(physics<NDIM>::variables_.sx_i, 1);
	}
	const auto nf = physics<NDIM>::field_count();
	std::vector<std::vector<std::vector<safe_real>>> F(NDIM, std::vector<std::vector<safe_real>>(nf, std::vector<safe_real>(H_N3)));
	std::vector<std::vector<safe_real>> U(nf, std::vector<safe_real>(H_N3));
	std::vector<std::vector<safe_real>> U0(nf, std::vector<safe_real>(H_N3));
	hydro::x_type X(NDIM);
	for (int dim = 0; dim < NDIM; dim++) {
		X[dim].resize(H_N3);
	}

	safe_real t = 0.0;
	int iter = 0;
	int oter = 0;
	physics<NDIM> phys;
	computer.set_bc(phys.template initialize<INX>(problem, U, X));
	const safe_real dx = X[0][cell_geometry<NDIM, INX>::H_DNX] - X[0][0];
	computer.output(U, X, oter++, 0);
//	const safe_real omega = 2.0 * M_PI / tmax / 10.0;
	const safe_real omega = 0.0;
	printf("omega = %e\n", (double) omega);

	const auto tstart = time(NULL);
	while (t < tmax) {
		U0 = U;
		auto q = computer.reconstruct(U, X, omega);
		first_q = q;
		safe_real a;
		if (with_kokkos){
	        a = octotiger::compute_flux_kokkos<NDIM, INX>(computer, U0, U, q, F, X, omega, nf, H_N3);
		}
		else{
			a = computer.flux(U, q, F, X, omega);
		}
		printState<NDIM>(U, U0, F, X, omega, q, a);
		safe_real dt = CFL * dx / a;
		dt = std::min(double(dt), tmax - t + 1.0e-20);
		computer.advance(U0, U, F, X, dx, dt, 1.0, omega);
		computer.boundaries(U);
		q = computer.reconstruct(U, X, omega);
		if (with_kokkos){
	        octotiger::compute_flux_kokkos<NDIM, INX>(computer, U0, U, q, F, X, omega, nf, H_N3);
		}
		else{
			computer.flux(U, q, F, X, omega);
		}
		computer.advance(U0, U, F, X, dx, dt, 0.25, omega);
		computer.boundaries(U);
		// printState<NDIM>(U, U0, F, X, omega, q, a);
		q = computer.reconstruct(U, X, omega);
		if (with_kokkos){
	        octotiger::compute_flux_kokkos<NDIM, INX>(computer, U0, U, q, F, X, omega, nf, H_N3);
		}
		else{
			computer.flux(U, q, F, X, omega);
		}
		computer.advance(U0, U, F, X, dx, dt, 2.0 / 3.0, omega);
		computer.boundaries(U);
		computer.post_process(U, dx);
		t += dt;
		computer.boundaries(U);
		if (int(t / dt_out) != int((t - dt) / dt_out))
			computer.output(U, X, oter++, t);
		iter++;
		printf("%i %e %e %e\n", iter, double(a), double(t), double(dt));
	}
	const auto tstop = time(NULL);
	U0 = U;
	physics<NDIM>::template analytic_solution<INX>(problem, U, X, t);
	computer.output(U, X, iter++, t);

	phys.template pre_recon<INX>(U0, X, omega, with_correction);
	phys.template pre_recon<INX>(U, X, omega, with_correction);
	std::vector<safe_real> L1(nf);
	std::vector<safe_real> L2(nf);
	std::vector<safe_real> Linf(nf);
	for (int f = 0; f < nf; f++) {
		L1[f] = L2[f] = Linf[f];
		for (int i = 0; i < H_N3; i++) {
			L1[f] += std::abs(U0[f][i] - U[f][i]);
			L2[f] += std::pow(U0[f][i] - U[f][i], 2);
			Linf[f] = max_device((double) Linf[f], std::abs(U0[f][i] - U[f][i]));
		}
		L2[f] = sqrt(L2[f]);
		L1[f] /= INX * INX;
		L2[f] /= INX * INX;
	}

	FILE *fp1 = fopen("L1.dat", "at");
	FILE *fp2 = fopen("L2.dat", "at");
	FILE *fpinf = fopen("Linf.dat", "at");
	fprintf(fp1, "%i ", INX);
	fprintf(fp2, "%i ", INX);
	fprintf(fpinf, "%i ", INX);
	for (int f = 0; f < nf; f++) {
		fprintf(fp1, "%e ", (double) L1[f]);
		fprintf(fp2, "%e ", (double) L2[f]);
		fprintf(fpinf, "%e ", (double) Linf[f]);
	}
	fprintf(fp1, "\n");
	fprintf(fp2, "\n");
	fprintf(fpinf, "\n");
	fclose(fp1);
	fclose(fp2);
	fclose(fpinf);
	FILE* fp = fopen( "time.dat", "wt");
	fprintf( fp, "%i %li\n", INX, tstop -tstart);
	fclose(fp);
}


template<int NDIM, int INX>
void compare(typename physics<NDIM>::test_type problem, bool with_correction) {
	static constexpr int H_BW = 3;
	// static constexpr auto H_NX = INX + 2 * H_BW;
	// static constexpr int H_N3 = std::pow(INX+2*H_BW,NDIM);
    static constexpr auto H_N3 = static_cast<int>(PowerNDIM(INX+2*H_BW));
	printf("test H_N3 %d, H_BW %d, NDIM %d, INX %d \n", H_N3, H_BW, NDIM, INX);
	// static_assert(std::pow(206, NDIM) == PowerNDIM(206));
	// static_assert(H_N3 == PowerNDIM(INX+2*H_BW));
	// assert(std::pow(206, NDIM) == PowerNDIM(206));
	// assert(H_N3 == PowerNDIM(INX+2*H_BW));

	hydro_computer<NDIM, INX> computer;
	if (with_correction) {
		computer.use_angmom_correction(physics<NDIM>::variables_.sx_i, 1);
	}
	const auto nf = physics<NDIM>::field_count();
	std::vector<std::vector<std::vector<safe_real>>> F_ref(NDIM, std::vector<std::vector<safe_real>>(nf, std::vector<safe_real>(H_N3)));
	auto F = F_ref;
	std::vector<std::vector<safe_real>> U(nf, std::vector<safe_real>(H_N3));
	std::vector<std::vector<safe_real>> U0(nf, std::vector<safe_real>(H_N3));
	hydro::x_type X(NDIM);
	for (int dim = 0; dim < NDIM; dim++) {
		X[dim].resize(H_N3);
	}

	int oter = 0;
	physics<NDIM> phys;
	computer.set_bc(phys.template initialize<INX>(problem, U, X));
	const safe_real dx = X[0][cell_geometry<NDIM, INX>::H_DNX] - X[0][0];
	computer.output(U, X, oter++, 0);
//	const safe_real omega = 2.0 * M_PI / tmax / 10.0;
	const safe_real omega = 0.0;
	printf("omega = %e\n", (double) omega);

	U0 = U;
	auto q = computer.reconstruct(U, X, omega);
	first_q = q;
	safe_real a_ref, a;
	a = octotiger::compute_flux_kokkos<NDIM, INX>(computer, U0, U, q, F, X, omega, nf, H_N3);
	a_ref = computer.flux(U, q, F_ref, X, omega);
	printf("a: %e %e \n", a, a_ref);

	// printState<NDIM>(U, U0, F, X, omega, q, a);
	for (int i=0; i < F.size(); ++i){
		for (int j=0; j < F[i].size(); ++j){
			for (int k=0; k < F[i][j].size(); ++k){
				auto diff = abs(F[i][j][k] - F_ref[i][j][k]);
				bool is_same = ((diff/F[i][j][k] < 1e-5) && (diff/F_ref[i][j][k] < 1e-5)) || diff < 1e-15;
				if (!is_same){
					printf("%d %d %d: %e %e; ", i, j, k, F[i][j][k], diff);
				}
			}
		}
    	printf("\n");
	}

	// static constexpr safe_real CFL = (0.4 / NDIM);
	// safe_real t = 0.;
	// safe_real dt = CFL * dx / a;
	// dt = std::min(double(dt), tmax - t + 1.0e-20);
	// computer.advance(U0, U, F, X, dx, dt, 1.0, omega);
	// computer.boundaries(U);
	// q = computer.reconstruct(U, X, omega);
}

// for timing purposes, test with random numbers
template <int NDIM, int INX>
void test_random_numbers() {
	static constexpr int H_BW = 3;
    static constexpr auto H_N3 = static_cast<int>(PowerNDIM(INX+2*H_BW));
	printf("test H_N3 %d, H_BW %d, NDIM %d, INX %d \n", H_N3, H_BW, NDIM, INX);
    const auto nf = physics<NDIM>::field_count();
    hydro_computer<NDIM, INX> computer;
    computer.use_angmom_correction(physics<NDIM>::variables_.sx_i, 1);
    const safe_real omega = 0.0;

    // use random data for the view fields
    hydro::x_type X(NDIM);
    for (int dim = 0; dim < NDIM; dim++) {
        X[dim].resize(H_N3);
    }

    std::vector<std::vector<std::vector<safe_real>>> F(
        NDIM, std::vector<std::vector<safe_real>>(nf, std::vector<safe_real>(H_N3)));
	auto F_comp = F;
    std::vector<std::vector<safe_real>> U(nf, std::vector<safe_real>(H_N3));
    std::vector<std::vector<safe_real>> U0(nf, std::vector<safe_real>(H_N3));
    std::vector<std::vector<std::array<safe_real, octotiger::q_lowest_dimension_length<NDIM>>>> Q(nf,
        std::vector<std::array<safe_real, octotiger::q_lowest_dimension_length<NDIM>>>(
            H_N3, std::array<safe_real, octotiger::q_lowest_dimension_length<NDIM>>()));

    // fill with non-zero entries to avoid divide by 0 error
    for (auto& vec : Q) {
        for (auto& arr : vec) {
            arr.fill(1.0);
        }
	}
	printf("compute flux with random numbers for %d, %d \n", NDIM, INX);
    for (int i = 0; i < 10; ++i) {
		for (auto & u : U){
			std::generate_n(u.begin(), H_N3, [](){
				return (std::rand() % 10) / 10. ; 
			});
		}
		for (auto & u : U0){
			std::generate_n(u.begin(), H_N3, [](){
				return (std::rand() % 10) / 10. ; 
			});
		}
		physics<NDIM> phys;
		computer.set_bc(phys.template initialize<INX>(physics<NDIM>::BLAST, U, X));
		Q = computer.reconstruct(U, X, omega);
		auto X_comp = X;

        auto a = octotiger::compute_flux_kokkos<NDIM, INX>(computer, U0, U, Q, F, X, omega, nf, H_N3);
		auto a_comp = computer.flux(U, Q, F_comp, X_comp, omega);

#if !defined(__CUDA_ARCH__)
		if (a_comp != a){
			printf("ERROR a %f %f \n", a_comp, a);
		}
		printf("%.5e; ", a);
		for (int i=0; i < NDIM; ++i){
			for (int j=0; j < nf; ++j){
				for (int k=0; k < H_N3; ++k){
					if( std::abs( F[i][j][k] - F_comp[i][j][k] ) > 1e-14 ){
						printf("F %.5e %.5e; ", F[i][j][k], F[i][j][k] - F_comp[i][j][k]);
					}
				}
			}
		}
		printf("\n");
		for (int i=0; i < NDIM; ++i){
			for (int k=0; k < H_N3; ++k){
				if( std::abs( X[i][k] - X_comp[i][k] ) > 1e-14 ){
					printf("X %f %f; ", X[i][k], X_comp[i][k]);
				}
			}
		}	
#endif 
	}
}

int main(int argc, char** argv) {

	Kokkos::initialize(argc, argv);
    // Kokkos::print_configuration(std::cout);

	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);

	// test_random_numbers<2, 10>();

#if !defined(__CUDA_ARCH__)
	// run_test<2, 200>(physics<2>::BLAST, true);
	// // 1 1.384974e+03 7.220354e-07 7.220354e-07
	// // 2 7.828761e+02 1.999377e-06 1.277341e-06
	// // 3 5.861212e+02 3.705509e-06 1.706132e-06
	// // 4 4.868369e+02 5.759584e-06 2.054076e-06
	// // 5 4.255020e+02 8.109750e-06 2.350165e-06
	// // 6 3.830523e+02 1.072036e-05 2.610609e-06
	// // 7 3.514781e+02 1.356549e-05 2.845127e-06
	// // 8 3.268058e+02 1.662541e-05 3.059921e-06
	// // 9 3.068242e+02 1.988460e-05 3.259195e-06
	// // 10 2.901991e+02 2.333051e-05 3.445910e-06
	// // 11 2.760700e+02 2.490000e-05 1.569487e-06

//	// run_test<3, 50>(physics<3>::BLAST, true);
//	// 1 1.309447e+03 2.036484e-06 2.036484e-06
//	// 2 5.083718e+02 7.281989e-06 5.245505e-06
//	// 3 3.312983e+02 1.533113e-05 8.049141e-06
//	// 4 2.520857e+02 2.490000e-05 9.568871e-06
#endif

	// <3, 50> apparently takes a lot of memory (or it is unevenly distributed)
 	// `mmap() failed to allocate thread stack due to insufficient resources, increase /proc/sys/vm/max_map_count 
	//  or add -Ihpx.stacks.use_guard_pages=0 to the command line: HPX(unhandled_exception)`
	// test_random_numbers<3, 50>();
	// test_random_numbers<2, 200>();
	// run_tests<3, 50, true>(physics<3>::BLAST, true);
	// run_test<2, 8, false>(physics<2>::BLAST, true);
	// run_test<2, 8, true>(physics<2>::BLAST, true); // identical
	// compare<2, 8>(physics<2>::BLAST, true);
	compare<2, 10>(physics<2>::BLAST, true);
	// compare<2, 20>(physics<2>::BLAST, true);
	// compare<2, 100>(physics<2>::BLAST, true);
	// compare<2, 200>(physics<2>::BLAST, true);
	// compare<2, 200>(physics<2>::BLAST, true);
	// run_test<2, 10, false>(physics<2>::BLAST, true);
	// run_test<2, 10, true>(physics<2>::BLAST, true);
	// flux_kokkos: 0.00442272 seconds
	
    Kokkos::finalize();
	return 0;
}

