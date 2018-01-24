#pragma once

#include "../common_kernel/multiindex.hpp"

#include <vector>

namespace octotiger {
namespace fmm {
    namespace multipole_interactions {

        two_phase_stencil calculate_stencil(void);

    }    // namespace multipole_interactions
}    // namespace fmm
}    // namespace octotiger
