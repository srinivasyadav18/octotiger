#ifdef OCTOTIGER_HAVE_CUDA

#include <buffer_manager.hpp>
#include <cuda_buffer_util.hpp>
#include "octotiger/options.hpp"
#include "octotiger/cuda_util/cuda_helper.hpp"
#include <cuda_runtime.h>
#include <stream_manager.hpp>

#include "octotiger/unitiger/hydro_impl/flux_kernel_interface.hpp"


// TODO Duplicated from cell_geometry class - come up with a better way to get this to the device
// Includes flip_dim, faces, xloc and quad_weights
// Move cell geometry class to device?

__device__ inline int flip_dim(const int d, const int flip_dim) {
		int dims[3];
		int k = d;
		for (int dim = 0; dim < 3; dim++) {
			dims[dim] = k % 3;
			k /= 3;
		}
		k = 0;
		dims[flip_dim] = 2 - dims[flip_dim];
		for (int dim = 0; dim < 3; dim++) {
			k *= 3;
			k += dims[2 - dim];
		}
		return k;
}

__device__ const int faces[3][9] = { { 12, 0, 3, 6, 9, 15, 18, 21, 24 }, { 10, 0, 1, 2, 9, 11,
			18, 19, 20 }, { 4, 0, 1, 2, 3, 5, 6, 7, 8 } };

__device__ const int xloc[27][3] = {
	/**/{ -1, -1, -1 }, { +0, -1, -1 }, { +1, -1, -1 },
	/**/{ -1, +0, -1 }, { +0, +0, -1 }, { 1, +0, -1 },
	/**/{ -1, +1, -1 }, { +0, +1, -1 }, { +1, +1, -1 },
	/**/{ -1, -1, +0 }, { +0, -1, +0 }, { +1, -1, +0 },
	/**/{ -1, +0, +0 }, { +0, +0, +0 }, { +1, +0, +0 },
	/**/{ -1, +1, +0 }, { +0, +1, +0 }, { +1, +1, +0 },
	/**/{ -1, -1, +1 }, { +0, -1, +1 }, { +1, -1, +1 },
	/**/{ -1, +0, +1 }, { +0, +0, +1 }, { +1, +0, +1 },
	/**/{ -1, +1, +1 }, { +0, +1, +1 }, { +1, +1, +1 } };

__device__ const double quad_weights[9] = { 16. / 36., 1. / 36., 4. / 36., 1. / 36., 4. / 36., 4.
			/ 36., 1. / 36., 4. / 36., 1. / 36. };

__device__ const int offset = 0;
__device__ const int compressedH_DN[3] = {100, 10, 1};
__device__ const int face_offset = 27 * 1000;
__device__ const int dim_offset = 1000;

__global__ void
__launch_bounds__(128, 2)
 flux_cuda_kernel(const double * __restrict__ q_combined, const double * __restrict__ x_combined, double * __restrict__ f_combined,
    double * amax, int * amax_indices, int * amax_d, const bool * __restrict__ masks, const double omega, const double dx, const double A_, const double B_, const int nf, const double fgamma, const double de_switch_1) {
  __shared__ double sm_amax[128];
  __shared__ int sm_d[128];
  __shared__ int sm_i[128];

//  const int nf = 15;

  // Set during cmake step with -DOCTOTIGER_WITH_MAX_NUMBER_FIELDS
  double local_f[OCTOTIGER_MAX_NUMBER_FIELDS];
	// assumes maximal number (given by cmake) of species in a simulation.
  // Not the most elegant solution and rather old-fashion but one that works.
  // May be changed to a more flexible sophisticated object.
  for (int f = 0; f < nf; f++) {   
      local_f[f] = 0.0;
  }
  double local_x[3] = {0.0, 0.0, 0.0};
  double local_vg[3] = {0.0, 0.0, 0.0};

  double current_amax = 0.0;
  int current_d = 0;

  // 3 dim 1000 i workitems
  const int dim = blockIdx.z;
  const int tid = threadIdx.x * 64 + threadIdx.y * 8 + threadIdx.z;
  const int index = blockIdx.y * 128 + tid + 104;
  for (int f = 0; f < nf; f++) {
      f_combined[dim * nf * 1000 + f * 1000 + index] = 0.0;
  }
  if (index < 1000) {
    double mask = masks[index + dim * dim_offset];
    if(mask != 0.0) {
      for (int fi = 0; fi < 9; fi++) {    // 9
        double this_ap = 0.0, this_am = 0.0;    // tmps
        const int d = faces[dim][fi];
        const int flipped_dim = flip_dim(d, dim);
        for (int dim = 0; dim < 3; dim++) {
            local_x[dim] = x_combined[dim * 1000 + index] + (0.5 * xloc[d][dim] * dx);
        }
        local_vg[0] = -omega * (x_combined[1000 + index] + 0.5 * xloc[d][1] * dx);
        local_vg[1] = +omega * (x_combined[index] + 0.5 * xloc[d][0] * dx);
        local_vg[2] = 0.0;
        inner_flux_loop2<double>(omega, nf, A_, B_, q_combined, local_f, local_x, local_vg,
          this_ap, this_am, dim, d, dx, fgamma, de_switch_1,
          dim_offset * d + index, dim_offset * flipped_dim - compressedH_DN[dim] + index, face_offset);
        this_ap *= mask;
        this_am *= mask;
        const double amax_tmp = max_wrapper(this_ap, (-this_am));
        if (amax_tmp > current_amax) {
          current_amax = amax_tmp;
          current_d = d;
        }
        for (int f = 1; f < nf; f++) {
          f_combined[dim * nf * 1000 + f * 1000 + index] += quad_weights[fi] * local_f[f];
        }
      }
    }
    for (int f = 10; f < nf; f++) {
      f_combined[dim * nf * 1000 + index] += f_combined[dim * nf * 1000 + f * 1000 + index];
    }
  }
  // Find maximum:
  sm_amax[tid] = current_amax;
  sm_d[tid] = current_d;
  sm_i[tid] = index;
  __syncthreads();
  // Max reduction with multiple warps
  for (int tid_border = 64; tid_border >= 32; tid_border /= 2) {
    if(tid < tid_border) {
      if (sm_amax[tid + tid_border] > sm_amax[tid]) {
        sm_amax[tid] = sm_amax[tid + tid_border];
        sm_d[tid] = sm_d[tid + tid_border];
        sm_i[tid] = sm_i[tid + tid_border];
      }
    }
    __syncthreads();
  }
  // Max reduction within one warps
  for (int tid_border = 16; tid_border >= 1; tid_border /= 2) {
  if(tid < tid_border) {
      if (sm_amax[tid + tid_border] > sm_amax[tid]) {
        sm_amax[tid] = sm_amax[tid + tid_border];
        sm_d[tid] = sm_d[tid + tid_border];
        sm_i[tid] = sm_i[tid + tid_border];
      }
    }
  }

  if (tid == 0) {
    //printf("Block %i %i TID %i %i \n", blockIdx.y, blockIdx.z, tid, index);
    const int block_id = blockIdx.y + dim * 7;
    amax[block_id] = sm_amax[0];
    amax_indices[block_id] = sm_i[0];
    amax_d[block_id] = sm_d[0];

    // Save face to the end of the amax buffer
    const int flipped_dim = flip_dim(sm_d[0], dim);
    for (int f = 0; f < nf; f++) {
      amax[21 + block_id * 2 * nf + f] = q_combined[sm_i[0] + f * face_offset + dim_offset * sm_d[0]];
      amax[21 + block_id * 2 * nf + nf + f] = q_combined[sm_i[0] - compressedH_DN[dim] + f * face_offset +
          dim_offset * flipped_dim];
    }
  }
  return;
}
void launch_flux_cuda_kernel_post(stream_interface<hpx::cuda::experimental::cuda_executor, pool_strategy>& executor,
    dim3 const grid_spec, dim3 const threads_per_block, void *args[]) {
    executor.post(
    cudaLaunchKernel<decltype(flux_cuda_kernel)>,
    flux_cuda_kernel, grid_spec, threads_per_block, args, 0);
}

#endif