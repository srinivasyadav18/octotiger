//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef OCTOTIGER_UNITIGER_UTI1L_HPP_
#define OCTOTIGER_UNITIGER_UTI1L_HPP_

#include "octotiger/unitiger/safe_real.hpp"

// to compile with nvcc, cannot use std::pow in constexpr; this is a workaround
#define PowerNDIM(num) (NDIM == 1 ? num : (NDIM == 2 ? num*num : (NDIM == 3 ? num*num*num : -1)))

template<int NDIM, int NX>
std::array<int, NDIM> index_to_dims(int i) {
	std::array<int, NDIM> dims;
	for (int j = 0; j < NDIM; j++) {
		dims[NDIM - 1 - j] = i % NX;
		i /= NX;
	}
	return dims;
}

template<int a, int b>
static constexpr int int_pow() {
	int c = 1;
	for( int i = 0; i < b; i++) {
		c *= a;
	}
	return a;
}

// #ifdef OCTOTIGER_WITH_KOKKOS

#include <Kokkos_Core.hpp>

// max function to use from CUDA and host
template<class T> 
KOKKOS_INLINE_FUNCTION const T& max_device(const T& a, const T& b)
{
    return (a < b) ? b : a; //TODO use CUDA intrinsic max if __CUDA_ARCH__
}

template<class T> 
KOKKOS_INLINE_FUNCTION const T& min_device(const T& a, const T& b)
{
    return (a < b) ? a : b;
}

// #endif

template<class T>
static inline void limit_slope(T &ql, T q0, T &qr) {
	const T tmp1 = qr - ql;
	const T tmp2 = qr + ql;

	if (bool(qr < q0) != bool(q0 < ql)) {
		qr = ql = q0;
		return;
	}
	const T tmp3 = tmp1 * tmp1 / 6.0;
	const T tmp4 = tmp1 * (q0 - 0.5 * tmp2);
	if (tmp4 > tmp3) {
		ql = 3.0 * q0 - 2.0 * qr;
	} else if (-tmp3 > tmp4) {
		qr = 3.0 * q0 - 2.0 * ql;
	}
}

static inline safe_real minmod(safe_real a, safe_real b) {
	return (std::copysign(0.5, a) + std::copysign(0.5, b)) * std::min(std::abs(a), std::abs(b));
}

static inline safe_real minmod_theta(safe_real a, safe_real b, safe_real c) {
	return minmod(c * minmod(a, b), 0.5 * (a + b));
}

#endif /* OCTOTIGER_UNITIGER_UTIL_HPP_ */
