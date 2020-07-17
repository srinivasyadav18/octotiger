//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#ifdef OCTOTIGER_HAVE_CUDA
#include "octotiger/common_kernel/interaction_constants.hpp"
#include "octotiger/common_kernel/multiindex.hpp"

#include <buffer_manager.hpp>

namespace octotiger {
namespace fmm {
    namespace monopole_interactions {

        // extern __constant__ double device_four_constants[4 * FULL_STENCIL_SIZE];
        // void kokkos_p2p_interactions_kernel(
        //     const double (&local_monopoles)[NUMBER_LOCAL_MONOPOLE_VALUES],
        //     double (&potential_expansions)[NUMBER_POT_EXPANSIONS_SMALL],
        //     const double theta, const double dx);
        void kokkos_p2p_interactions(std::vector<real, recycler::recycle_std<real>>& buffer,
            std::vector<real, recycler::recycle_std<real>>& output, const double dx,
            const double theta, bool *masks);
    }    // namespace monopole_interactions
}    // namespace fmm
}    // namespace octotiger
#endif