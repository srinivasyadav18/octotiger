#pragma once

#include "octotiger/cuda_util/cuda_global_def.hpp"

#define OCTOTIGER_WITH_KOKKOS //TODO remove

#if defined(OCTOTIGER_WITH_KOKKOS)
  #include <Kokkos_Macros.hpp>
  #define OCTOTIGER_FUNCTION KOKKOS_FUNCTION
#else /* defined(OCTOTIGER_WITH_KOKKOS) */
  #define OCTOTIGER_FUNCTION CUDA_CALLABLE_METHOD
#endif /* defined(OCTOTIGER_WITH_KOKKOS) */