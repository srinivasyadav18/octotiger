//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <fenv.h>
#include <time.h>
#include <type_traits>
#include <typeinfo>

#define OCTOTIGER_GRIDDIM 8    // usually set in build scripts

#include "flux_kokkos.hpp"

#include <hpx/hpx_main.hpp>

#include "octotiger/unitiger/safe_real.hpp"
#include "octotiger/unitiger/unitiger.hpp"

static constexpr double tmax = 2.49e-5;
static constexpr safe_real dt_out = tmax / 249;

// #define H_BW 3
// #define H_NX (INX + 2 * H_BW)
// #define H_N3 std::pow(INX+2*H_BW,NDIM)

// why is there two different INX?
// one is constexpr in hydro_defs.hpp - used with unitiger
// one the template parameter used here 
// ask Dominic again - how are they different?

template<int NDIM, int INX>
void run_test(typename physics<NDIM>::test_type problem, bool with_correction) {
	static constexpr auto H_BW = 3;
	static constexpr auto H_NX = INX + 2 * H_BW;
	static constexpr auto H_N3 = std::pow(INX+2*H_BW,NDIM);

	static constexpr safe_real CFL = (0.4 / NDIM);
	hydro_computer<NDIM, INX> computer;
	if (with_correction) {
		computer.use_angmom_correction(physics<NDIM>::sx_i, 1);
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
		auto a = computer.flux(U, q, F, X, omega);
		safe_real dt = CFL * dx / a;
		dt = std::min(double(dt), tmax - t + 1.0e-20);
		computer.advance(U0, U, F, X, dx, dt, 1.0, omega);
		computer.boundaries(U);
		q = computer.reconstruct(U, X, omega);
		computer.flux(U, q, F, X, omega);
		computer.advance(U0, U, F, X, dx, dt, 0.25, omega);
		computer.boundaries(U);
		q = computer.reconstruct(U, X, omega);
		computer.flux(U, q, F, X, omega);
		computer.advance(U0, U, F, X, dx, dt, 2.0 / 3.0, omega);
		computer.boundaries(U);
		computer.post_process(U, dx);
		t += dt;
		computer.boundaries(U);
		if (int(t / dt_out) != int((t - dt) / dt_out))
			computer.output(U, X, oter++, t);
		iter++;
		printf("%i %e %e\n", iter, double(t), double(dt));
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
			Linf[f] = std::max((double) Linf[f], std::abs(U0[f][i] - U[f][i]));
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

template <int NDIM>
static constexpr auto q_lowest_dimension_length = NDIM == 1 ? 3 : (NDIM == 2 ? 9 : 27);

template <int NDIM, int INX>
safe_real compute_flux_kokkos(const hydro_computer<NDIM, INX>& computer,
    const std::vector<std::vector<safe_real>>& U0, const std::vector<std::vector<safe_real>>& U,
    const std::vector<std::vector<std::array<safe_real, q_lowest_dimension_length<NDIM>>>>& q,
    std::vector<std::vector<std::vector<safe_real>>>& F,
    const std::vector<std::vector<safe_real>>& X, const safe_real omega, const int nf,
    const int H_N3) {
    // do this iteration in the kokkosified version
    // F is an output of the flux kernel, zero-initialize
    Kokkos::View<safe_real***> kokkosF("flux", NDIM, nf, H_N3);    // TODO ask dominic if these are the right dimensions
	auto kokkosFhost = Kokkos::create_mirror_view(kokkosF);

    Kokkos::View<safe_real**> kokkosU(Kokkos::ViewAllocateWithoutInitializing("state"), nf, H_N3);
    Kokkos::View<safe_real**> kokkosU0(
        Kokkos::ViewAllocateWithoutInitializing("initial state"), nf, H_N3);
	auto kokkosUhost = Kokkos::create_mirror_view(kokkosU);
	auto kokkosU0host = Kokkos::create_mirror_view(kokkosU0);
    Kokkos::parallel_for("init_U", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nf, H_N3}),
        KOKKOS_LAMBDA(int j, int k) {
            kokkosUhost(j, k) = U[j][k];
            kokkosU0host(j, k) = U0[j][k];
        });

    Kokkos::View<safe_real**> kokkosX(
        Kokkos::ViewAllocateWithoutInitializing("cell center coordinates"), NDIM, H_N3);
	auto kokkosXhost = Kokkos::create_mirror_view(kokkosX);
    Kokkos::parallel_for("init_X", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {NDIM, H_N3}),
        KOKKOS_LAMBDA(int i, int k) { kokkosXhost(i, k) = X[i][k]; });

    Kokkos::View<safe_real* * [q_lowest_dimension_length<NDIM>]> kokkosQ(
        Kokkos::ViewAllocateWithoutInitializing("reconstruction"), nf, H_N3);
	auto kokkosQhost = Kokkos::create_mirror_view(kokkosQ);
	// have this run in serial, takes forever otherwise
    Kokkos::parallel_for("init_Q",
        Kokkos::MDRangePolicy<Kokkos::Serial, Kokkos::Rank<3>>({0, 0, 0}, {nf, H_N3, q_lowest_dimension_length<NDIM>}),
        KOKKOS_LAMBDA(int i, int j, int k) { kokkosQhost(i, j, k) = q[i][j][k]; });

    // std::cout << kokkosQ.extent(0) << " " << kokkosQ.extent(1) << " " << kokkosQ.extent(2) << "Q" << std::endl;
	Kokkos::fence();

	Kokkos::deep_copy(kokkosUhost, kokkosU);
	Kokkos::deep_copy(kokkosU0host, kokkosU0);
	Kokkos::deep_copy(kokkosXhost, kokkosX);
	Kokkos::deep_copy(kokkosQhost, kokkosQ);

    Kokkos::fence();

    // computer.flux(U, q, F, X, omega);
    auto amax = octotiger::flux_kokkos_hpx(computer, kokkosU, kokkosQ, kokkosF, kokkosX, omega);

	Kokkos::fence();
	
	Kokkos::deep_copy (kokkosF, kokkosFhost);

    Kokkos::fence();
    // copy back the values obtained for F
    for (int i = 0; i < NDIM; ++i) {
        for (int j = 0; j < nf; ++j) {
            for (int k = 0; k < H_N3; ++k) {
                F[i][j][k] = kokkosFhost(i, j, k);
            }
        }
    }
    return amax;
}

// for timing purposes, test with random numbers
template <int NDIM, int INX>
void test_random_numbers() {
    static constexpr auto H_N3 = static_cast<int>(std::pow(INX + 2 * H_BW, NDIM));
    const auto nf = physics<NDIM>::field_count();
    hydro_computer<NDIM, INX> computer;
    computer.use_angmom_correction(physics<NDIM>::sx_i, 1);
    const safe_real omega = 0.0;

    // use random data for the view fields
    hydro::x_type X(NDIM);
    for (int dim = 0; dim < NDIM; dim++) {
        X[dim].resize(H_N3);
    }
    std::vector<std::vector<std::vector<safe_real>>> F(
        NDIM, std::vector<std::vector<safe_real>>(nf, std::vector<safe_real>(H_N3)));
    std::vector<std::vector<safe_real>> U(nf, std::vector<safe_real>(H_N3));
    std::vector<std::vector<safe_real>> U0(nf, std::vector<safe_real>(H_N3));
    std::vector<std::vector<std::array<safe_real, q_lowest_dimension_length<NDIM>>>> Q(nf,
        std::vector<std::array<safe_real, q_lowest_dimension_length<NDIM>>>(
            H_N3, std::array<safe_real, q_lowest_dimension_length<NDIM>>()));

    // fill with non-zero entries to avoid divide by 0 error
    for (auto& vec : Q) {
        for (auto& arr : vec) {
            arr.fill(1.0);
        }
    }

    for (int i = 0; i < 10; ++i) {
        compute_flux_kokkos(computer, U0, U, Q, F, X, omega, nf, H_N3);
    }
}

template <int NDIM, int INX>
void run_test_kokkos(typename physics<NDIM>::test_type problem, bool with_correction) {
	static constexpr auto H_BW = 3;
	static constexpr auto H_NX = INX + 2 * H_BW;
    static constexpr auto H_N3 = static_cast<int>(std::pow(INX + 2 * H_BW, NDIM));

	static constexpr safe_real CFL = (0.4 / NDIM);
	hydro_computer<NDIM, INX> computer;
	if (with_correction) {
		computer.use_angmom_correction(physics<NDIM>::sx_i, 1);
	}
	const auto nf = physics<NDIM>::field_count();
	// assert all our indices are integral
	static_assert(std::is_integral<decltype(NDIM)>::value);
	static_assert(std::is_integral<decltype(nf)>::value);
	static_assert(std::is_integral<decltype(H_N3)>::value);
	
    std::vector<std::vector<std::vector<safe_real>>> F(
        NDIM, std::vector<std::vector<safe_real>>(nf, std::vector<safe_real>(H_N3)));
	std::vector<std::vector<safe_real>> U(nf, std::vector<safe_real>(H_N3));
	std::vector<std::vector<safe_real>> U0(nf, std::vector<safe_real>(H_N3));

	// x_type = std::vector<std::vector<safe_real>>;
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
		// this is the return type of q
		// recon_type =std::vector<std::vector<std::array<safe_real, NDIM == 1 ? 3 : (NDIM == 2 ? 9 : 27)>>>;
		auto q = computer.reconstruct(U, X, omega);
		auto a = computer.flux(U, q, F, X, omega);
		safe_real dt = CFL * dx / a;
		dt = std::min(double(dt), tmax - t + 1.0e-20);
		computer.advance(U0, U, F, X, dx, dt, 1.0, omega);
		computer.boundaries(U);
		q = computer.reconstruct(U, X, omega);
        compute_flux_kokkos(computer, U0, U, q, F, X, omega, nf, H_N3);
		computer.advance(U0, U, F, X, dx, dt, 0.25, omega);
		computer.boundaries(U);
		q = computer.reconstruct(U, X, omega);
		computer.flux(U, q, F, X, omega);
		computer.advance(U0, U, F, X, dx, dt, 2.0 / 3.0, omega);
		computer.boundaries(U);
		computer.post_process(U, dx);
		t += dt;
		computer.boundaries(U);
		if (int(t / dt_out) != int((t - dt) / dt_out))
			computer.output(U, X, oter++, t);
		iter++;
		printf("%i %e %e\n", iter, double(t), double(dt));
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
			Linf[f] = std::max((double) Linf[f], std::abs(U0[f][i] - U[f][i]));
		}
		L2[f] = sqrt(L2[f]);
		L1[f] /= INX * INX;
		L2[f] /= INX * INX;
	}

	FILE *fp1 = fopen("L1_k.dat", "at");
	FILE *fp2 = fopen("L2_k.dat", "at");
	FILE *fpinf = fopen("Linf_k.dat", "at");
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

int main(int argc, char** argv) {
    Kokkos::initialize(argc, argv);
    Kokkos::print_configuration(std::cout);

	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);

	test_random_numbers<3, 50>();
	run_test_kokkos<3, 50>(physics<3>::BLAST, true);
	// run_test<2, 200>(physics<2>::BLAST, true);

    Kokkos::finalize();
	return 0;
}