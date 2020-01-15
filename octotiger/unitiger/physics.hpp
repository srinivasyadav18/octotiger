//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../cuda_util/cuda_global_def.hpp"
#include "../test_problems/blast.hpp"
#include "../test_problems/exact_sod.hpp"
#include "./safe_real.hpp"
#include "./util.hpp"

#ifndef OCTOTIGER_UNITIGER_PHYSICS_HPP_
#define OCTOTIGER_UNITIGER_PHYSICS_HPP_

template <int NDIM>
struct physics_variables
{
    const int rho_i = 0;
    const int egas_i = 1;
    const int tau_i = 2;
    const int pot_i = 3;
    const int sx_i = 4;
    const int sy_i = 5;
    const int sz_i = 6;
    const int zx_i = 4 + NDIM;
    const int zy_i = 5 + NDIM;
    const int zz_i = 6 + NDIM;
    // static constexpr int spc_i = 4 + NDIM + (NDIM == 1 ? 0 : std::pow(3, NDIM - 2));
    const int spc_i = 4 + NDIM + (NDIM == 1 ? 0 : (NDIM == 2 ? 1 : (NDIM == 3 ? 3 : -1)));
    bool angmom_ = false;
    safe_real fgamma_ = 7. / 5.;
    static const int n_species_ = 5;
    // static const int nf_ = (4 + NDIM + (NDIM == 1 ? 0 : std::pow(3, NDIM - 2))) + n_species_;
    static const int nf_ =
        (4 + NDIM + (NDIM == 1 ? 0 : (NDIM == 2 ? 1 : (NDIM == 3 ? 3 : -1)))) + n_species_;
};

template <int NDIM>
struct physics
{
    static physics_variables<NDIM> variables_;

    enum test_type
    {
        SOD,
        BLAST,
        KH,
        CONTACT
    };

    static constexpr int field_count() {
        return variables_.nf_;
    }

    static void set_fgamma(safe_real fg) {
        variables_.fgamma_ = fg;
    }

    // static void set_n_species(int n) {
    // 	n_species_ = n;
    // }

    static void set_angmom() {
        variables_.angmom_ = true;
    }

    static void to_prim(
        const std::vector<safe_real>& u, safe_real& p, safe_real& v, int dim, safe_real dx) {
        const auto rho = u[variables_.rho_i];
        const auto rhoinv = safe_real(1.) / rho;
        safe_real ek = 0.0;
        for (int dim = 0; dim < NDIM; dim++) {
            ek += pow(u[variables_.sx_i + dim], 2) * rhoinv * safe_real(0.5);
        }
        auto ein = u[variables_.egas_i] - ek;
        if (ein < safe_real(0.001) * u[variables_.egas_i]) {
            ein = pow(u[variables_.tau_i], variables_.fgamma_);
        }

        v = u[variables_.sx_i + dim] * rhoinv;
        p = (variables_.fgamma_ - 1.0) * ein;
    }

    static void physical_flux(const std::vector<safe_real>& U, std::vector<safe_real>& F, int dim,
        safe_real& am, safe_real& ap, const std::array<safe_real, NDIM>& vg, safe_real dx) {
        safe_real p, v, v0, c;
        to_prim(U, p, v0, dim, dx);
        v = v0 - vg[dim];
        c = std::sqrt(variables_.fgamma_ * p / U[variables_.rho_i]);
        am = v - c;
        ap = v + c;
        for (int f = 0; f < variables_.nf_; f++) {
            F[f] = v * U[f];
        }
        F[variables_.sx_i + dim] += p;
        F[variables_.egas_i] += v0 * p;
    }

    static void flux(const std::vector<safe_real>& UL, const std::vector<safe_real>& UR,
        const std::vector<safe_real>& UL0, const std::vector<safe_real>& UR0,
        std::vector<safe_real>& F, int dim, safe_real& am, safe_real& ap,
        std::array<safe_real, NDIM>& vg, safe_real dx) {
        // safe_real pr, vr, pl, vl, vr0, vl0, amr, apr, aml, apl;
        safe_real amr, apr, aml, apl;

        static thread_local std::vector<safe_real> FR(variables_.nf_), FL(variables_.nf_);

        physical_flux(UR, FR, dim, amr, apr, vg, dx);
        physical_flux(UL, FL, dim, aml, apl, vg, dx);
        ap = std::max(std::max(apr, apl), safe_real(0.0));
        am = std::min(std::min(amr, aml), safe_real(0.0));
        for (int f = 0; f < variables_.nf_; f++) {
            F[f] = (ap * FL[f] - am * FR[f] + ap * am * (UR[f] - UL[f])) / (ap - am);
        }
    }

    template <int INX>
    static void post_process(hydro::state_type& U, safe_real dx) {
        static const cell_geometry<NDIM, INX> geo;
        constexpr auto dir = geo.direction();
        const static auto is = geo.find_indices(geo.H_BW, geo.H_NX - geo.H_BW);
        for (auto i : is) {
            safe_real ek = 0.0;
            for (int dim = 0; dim < NDIM; dim++) {
                ek += U[variables_.sx_i + dim][i] * U[variables_.sx_i + dim][i];
            }
            ek *= 0.5 * INVERSE(U[variables_.rho_i][i]);
            auto egas_max = U[variables_.egas_i][i];
            for (int d = 0; d < geo.NDIR; d++) {
                egas_max = std::max(egas_max, U[variables_.egas_i][i + dir[d]]);
            }
            safe_real ein = U[variables_.egas_i][i] - ek;
            if (ein > 0.1 * egas_max) {
                U[variables_.tau_i][i] = POWER(ein, 1.0 / variables_.fgamma_);
            }
        }
    }

    template <int INX>
    static void source(hydro::state_type& dudt, const hydro::state_type& U,
        const hydro::flux_type& F, const hydro::x_type X, safe_real omega, safe_real dx) {
        static const cell_geometry<NDIM, INX> geo;
        for (int dim = 0; dim < NDIM; dim++) {
            static constexpr auto kdelta = geo.kronecker_delta();
            for (int n = 0; n < geo.NANGMOM; n++) {
                const auto m = dim;
                for (int l = 0; l < NDIM; l++) {
                    for (const auto& i : geo.find_indices(geo.H_BW, geo.H_NX - geo.H_BW)) {
                        const auto fr = F[dim][variables_.sx_i + l][i + geo.H_DN[dim]];
                        const auto fl = F[dim][variables_.sx_i + l][i];
                        dudt[variables_.zx_i + n][i] -= kdelta[n][m][l] * 0.5 * (fr + fl);
                    }
                }
            }
        }
        for (const auto& i : geo.find_indices(geo.H_BW, geo.H_NX - geo.H_BW)) {
            if
                CONSTEXPR(NDIM == 3) {
                    dudt[variables_.zx_i][i] -= omega * X[2][i] * U[variables_.sx_i][i];
                    dudt[variables_.zy_i][i] -= omega * X[2][i] * U[variables_.sy_i][i];
                }
            if
                CONSTEXPR(NDIM >= 2) {
                    dudt[variables_.zx_i][i] +=
                        omega * (X[0][i] * U[variables_.sx_i][i] + X[1][i] * U[variables_.sy_i][i]);
                }
        }
        for (const auto& i : geo.find_indices(geo.H_BW, geo.H_NX - geo.H_BW)) {
            dudt[variables_.sx_i][i] += U[variables_.sy_i][i] * omega;
            dudt[variables_.sy_i][i] -= U[variables_.sx_i][i] * omega;
        }
    }

    template <int INX>
    static void pre_angmom(const hydro::state_type& U, const hydro::recon_type<NDIM>& Q,
        std::array<safe_real, cell_geometry<NDIM, INX>::NANGMOM>& Z,
        std::array<std::array<safe_real, cell_geometry<NDIM, INX>::NDIR>, NDIM>& S, int i,
        double dx) {
        static const cell_geometry<NDIM, INX> geo;
        static const auto indices = geo.find_indices(1, geo.H_NX - 1);
        for (int d = 0; d < geo.NDIR; d++) {
            if (d != geo.NDIR / 2) {
                const auto rho = Q[variables_.rho_i][i][d];
                for (int f = 0; f < NDIM; f++) {
                    S[f][d] *= rho;
                }
            }
        }
        for (int f = 0; f < geo.NANGMOM; f++) {
            const auto rho = U[variables_.rho_i][i];
            Z[f] *= rho;
        }
    }

    template <int INX>
    static void post_angmom(const hydro::state_type& U, const hydro::recon_type<NDIM>& Q,
        std::array<safe_real, cell_geometry<NDIM, INX>::NANGMOM>& Z,
        std::array<std::array<safe_real, cell_geometry<NDIM, INX>::NDIR>, NDIM>& S, int i,
        double dx) {
        static const cell_geometry<NDIM, INX> geo;
        static const auto indices = geo.find_indices(1, geo.H_NX - 1);
        for (int d = 0; d < geo.NDIR; d++) {
            if (d != geo.NDIR / 2) {
                const auto rho = Q[variables_.rho_i][i][d];
                for (int f = 0; f < NDIM; f++) {
                    S[f][d] /= rho;
                }
            }
        }
        for (int f = 0; f < geo.NANGMOM; f++) {
            const auto rho = U[variables_.rho_i][i];
            Z[f] /= rho;
        }
    }

    template <int INX>
    static const hydro::state_type pre_recon(
        const hydro::state_type& U, const hydro::x_type X, safe_real omega, bool angmom) {
        static const cell_geometry<NDIM, INX> geo;
        static const auto indices = geo.find_indices(0, geo.H_NX);
        auto V = U;
        const auto dx = X[0][geo.H_DNX] - X[0][0];
        for (const auto& i : indices) {
            const auto rho = V[variables_.rho_i][i];
            const auto rhoinv = 1.0 / rho;
            for (int dim = 0; dim < NDIM; dim++) {
                auto& s = V[variables_.sx_i + dim][i];
                V[variables_.egas_i][i] -= 0.5 * s * s * rhoinv;
                s *= rhoinv;
            }
            for (int si = 0; si < variables_.n_species_; si++) {
                V[variables_.spc_i + si][i] *= rhoinv;
            }
            V[variables_.pot_i][i] *= rhoinv;
        }
        return V;
    }

    template <int INX>
    static hydro::recon_type<NDIM> post_recon(
        const hydro::recon_type<NDIM>& P, const hydro::x_type X, safe_real omega, bool angmom) {
        static const cell_geometry<NDIM, INX> geo;
        static const auto indices = geo.find_indices(2, geo.H_NX - 2);
        auto Q = P;
        const auto dx = X[0][geo.H_DNX] - X[0][0];
        for (const auto& i : indices) {
            for (int d = 0; d < geo.NDIR; d++) {
                if (d != geo.NDIR / 2) {
                    const auto rho = Q[variables_.rho_i][i][d];
                    for (int dim = 0; dim < NDIM; dim++) {
                        auto& v = Q[variables_.sx_i + dim][i][d];
                        Q[variables_.egas_i][i][d] += 0.5 * v * v * rho;
                        v *= rho;
                    }
                    Q[variables_.pot_i][i][d] *= rho;
                    safe_real w = 0.0;
                    for (int si = 0; si < variables_.n_species_; si++) {
                        w += Q[variables_.spc_i + si][i][d];
                        Q[variables_.spc_i + si][i][d] *= rho;
                    }
                    if (w == 0.0) {
                        printf("NO SPECIES %i\n", i);
                        //		sleep(10);
                        abort();
                    }
                    w = 1.0 / w;
                    for (int si = 0; si < variables_.n_species_; si++) {
                        Q[variables_.spc_i + si][i][d] *= w;
                    }
                }
            }
        }
        return Q;
    }

    template <int INX>
    static void analytic_solution(
        test_type test, hydro::state_type& U, const hydro::x_type& X, safe_real time) {
        static const cell_geometry<NDIM, INX> geo;
        static safe_real rmax = 0.0;
        static std::once_flag one;

        std::call_once(one, [&X]() {
            for (int i = 0; i < geo.H_N3; i++) {
                safe_real r = 0.0;
                for (int dim = 0; dim < NDIM; dim++) {
                    r += X[dim][i] * X[dim][i];
                }
                r = sqrt(r);
                rmax = std::max(r, rmax);
            }
            rmax *= 2.0;
        });

        for (int f = 0; f < variables_.nf_; f++) {
            for (auto& u : U[f]) {
                u = 0.0;
            }
        }

        for (int i = 0; i < geo.H_N3; i++) {
            safe_real r = 0.0;
            safe_real rsum = 0.0;
            for (int dim = 0; dim < NDIM; dim++) {
                r += X[dim][i] * X[dim][i];
                rsum += X[dim][i];
            }
            r = sqrt(r);
            double den = 0, vel = 0, pre = 0;

            if (test == BLAST) {
                sedov::solution(time + 7e-4, r, rmax, den, vel, pre, NDIM);
                for (int dim = 0; dim < NDIM; dim++) {
                    U[variables_.sx_i + dim][i] = den * vel * X[dim][i] / r;
                }
            } else if (test == SOD) {
                sod_state_t sod_state;
                exact_sod(&sod_state, &sod_init, rsum / std::sqrt(NDIM), time);
                den = sod_state.rho;
                vel = sod_state.v;
                for (int dim = 0; dim < NDIM; dim++) {
                    U[variables_.sx_i + dim][i] = den * vel / std::sqrt(NDIM);
                }
            } else if (test == CONTACT) {
                pre = 1.0;
                vel = 10.0;
                den = 1.0 + 1.0e-6 * sin(2.0 * M_PI * (X[0][i] - vel * time));
            }

            U[variables_.rho_i][i] = den;
            U[variables_.tau_i][i] =
                pow(pre / (variables_.fgamma_ - 1.0), 1.0 / variables_.fgamma_);
            U[variables_.egas_i][i] = pre / (variables_.fgamma_ - 1.0) + 0.5 * den * vel * vel;
            U[variables_.spc_i][i] = den;
        }
    }

    template <int INX>
    std::vector<typename hydro_computer<NDIM, INX>::bc_type> initialize(
        test_type t, hydro::state_type& U, hydro::x_type& X) {
        static const cell_geometry<NDIM, INX> geo;

        std::vector<typename hydro_computer<NDIM, INX>::bc_type> bc(2 * NDIM);

        for (int i = 0; i < 2 * NDIM; i++) {
            bc[i] = hydro_computer<NDIM, INX>::OUTFLOW;
        }

        switch (t) {
        case SOD:
        case BLAST:
            break;
        case KH:
        case CONTACT:
            for (int i = 0; i < 2 * NDIM; i++) {
                bc[i] = hydro_computer<NDIM, INX>::PERIODIC;
            }
            break;
        }
        U.resize(variables_.nf_);
        for (int dim = 0; dim < NDIM; dim++) {
            X[dim].resize(geo.H_N3);
        }
        for (int f = 0; f < variables_.nf_; f++) {
            U[f].resize(geo.H_N3, 0.0);
        }

        const safe_real dx = 1.0 / INX;

        for (int i = 0; i < geo.H_N3; i++) {
            int k = i;
            int j = 0;
            for (int dim = 0; dim < NDIM; dim++) {
                X[NDIM - 1 - dim][i] = (((k % geo.H_NX) - geo.H_BW) + 0.5) * dx - 0.5;
                k /= geo.H_NX;
                j++;
            }
        }

        for (int i = 0; i < geo.H_N3; i++) {
            for (int f = 0; f < variables_.nf_; f++) {
                U[f][i] = 0.0;
            }
            // const auto xlocs = geo.xloc();
            // const auto weights = geo.volume_weight();
            std::array<safe_real, NDIM> x;
            double rho = 0, vx = 0, vy = 0, vz = 0, p = 0, r;
            safe_real x2, xsum;
            for (int dim = 0; dim < NDIM; dim++) {
                x[dim] = X[dim][i];
            }
            xsum = 0.0;
            for (int dim = 0; dim < NDIM; dim++) {
                xsum += x[dim];
            }
            switch (t) {
            case CONTACT:
                p = 1.0;
                vx = 10.0;
                rho = 1.0 + 1.0e-6 * sin(2.0 * M_PI * x[0]);
                break;
            case SOD:
                if (xsum < 0) {
                    rho = 1.0;
                    p = 1.0;
                } else {
                    rho = 0.125;
                    p = 0.1;
                }
                break;
            case BLAST:

                x2 = 0.0;
                for (int dim = 0; dim < NDIM; dim++) {
                    x2 += x[dim] * x[dim];
                }
                r = sqrt(x2);
                double v;
                sedov::solution(7e-4, r, std::sqrt(3) + 5.0 * dx, rho, v, p, NDIM);
                p = std::max((variables_.fgamma_ - 1.0) * 1.0e-20, p);
                //				/***************/
                //				if( r < 1.5*dx ) {
                //					p = 1.0e+3;
                //				} else {
                //					p = 1.0e-3;
                //				}
                //				rho = 1.0;
                //				/**************/
                vx = v * X[0][i] / r;
                if
                    CONSTEXPR(NDIM >= 2) {
                        vy = v * X[1][i] / r;
                    }
                if
                    CONSTEXPR(NDIM == 3) {
                        vz = v * X[2][i] / r;
                    }
                break;
            case KH:
                const auto eps = []() { return (rand() + 0.5) / RAND_MAX * 1.0e-3; };

                U[variables_.tau_i][i] = 1.0;
                p = 1.0;
                if (x[1] < 0.0) {
                    rho = 1.0 + eps();
                    vx = -0.5;
                } else {
                    rho = 2.0 + eps();
                    vx = +0.5;
                }
                break;
            }
            U[variables_.rho_i][i] += rho;
            U[variables_.spc_i][i] += rho;
            U[variables_.sx_i][i] += (rho * vx);
            U[variables_.egas_i][i] += (p / (variables_.fgamma_ - 1.0) + 0.5 * rho * vx * vx);
            U[variables_.tau_i][i] +=
                (std::pow(p / (variables_.fgamma_ - 1.0), 1.0 / variables_.fgamma_));
            if
                CONSTEXPR(NDIM >= 2) {
                    U[variables_.sy_i][i] += rho * vy;
                    U[variables_.egas_i][i] += 0.5 * rho * vy * vy;
                }
            if
                CONSTEXPR(NDIM >= 3) {
                    U[variables_.sz_i][i] += rho * vz;
                    U[variables_.egas_i][i] += 0.5 * rho * vz * vz;
                }
        }

        return bc;
    }
};

// copy of the three member functions for kokkos view data type and a data struct to be handed over
template <int NDIM, typename ViewSlice>
KOKKOS_INLINE_FUNCTION void to_prim(const physics_variables<NDIM> var, const ViewSlice u,
    safe_real& p, safe_real& v, int dim, safe_real dx) {
    const auto rho = u[var.rho_i];
    const auto rhoinv = safe_real(1.) / rho;
    safe_real ek = 0.0;
    for (int dim = 0; dim < NDIM; dim++) {
        ek += pow(u[var.sx_i + dim], 2) * rhoinv * safe_real(0.5);
    }
    auto ein = u[var.egas_i] - ek;
    if (ein < safe_real(0.001) * u[var.egas_i]) {
        ein = pow(u[var.tau_i], var.fgamma_);
    }

    v = u[var.sx_i + dim] * rhoinv;
    p = (var.fgamma_ - 1.0) * ein;
}

template <int NDIM, typename ViewSlice>
KOKKOS_INLINE_FUNCTION void physical_flux(const physics_variables<NDIM> var,
    const ViewSlice U, safe_real* /*safe_real[field_count()]*/ F, int dim, safe_real& am,
    safe_real& ap, const safe_real* /*safe_real[NDIM]*/ vg, safe_real dx) {
    safe_real p, v, v0, c;
    to_prim(var, U, p, v0, dim, dx);
    v = v0 - vg[dim];
    c = std::sqrt(var.fgamma_ * p / U[var.rho_i]);
    am = v - c;
    ap = v + c;
    for (int f = 0; f < var.nf_; f++) {
        F[f] = v * U[f];
    }
    F[var.sx_i + dim] += p;
    F[var.egas_i] += v0 * p;
}

template <int NDIM, typename ViewSlice>
KOKKOS_INLINE_FUNCTION void flux(const physics_variables<NDIM> var, const ViewSlice UL,
    const ViewSlice UR,    // const ViewSlice UL0, const ViewSlice UR0,
    safe_real* /*safe_real[field_count()]*/ F, int dim, safe_real& am, safe_real& ap,
    const safe_real* /*safe_real[NDIM]*/ vg, safe_real dx) {
    // safe_real pr, vr, pl, vl, vr0, vl0, amr, apr, aml, apl;
    safe_real amr, apr, aml, apl;

    // fill fr and fl
	//TODO test whether this does the same thing
    safe_real FR[var.nf_];
    safe_real FL[var.nf_];

    physical_flux(var, UR, FR, dim, amr, apr, vg, dx);
    physical_flux(var, UL, FL, dim, aml, apl, vg, dx);
    ap = max_device(max_device(apr, apl), safe_real(0.0));
    am = min_device(min_device(amr, aml), safe_real(0.0));
    for (int f = 0; f < var.nf_; f++) {
        F[f] = (ap * FL[f] - am * FR[f] + ap * am * (UR[f] - UL[f])) / (ap - am);
    }
}

#endif /* OCTOTIGER_UNITIGER_PHYSICS_HPP_ */
