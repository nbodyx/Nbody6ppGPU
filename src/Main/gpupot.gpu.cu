//#include <iostream>
#include <cstdio>
// #include <cutil.h>
#ifdef CUDA_5
#  include <helper_cuda.h>
#  define CUDA_SAFE_CALL checkCudaErrors
#else
#  include <cutil.h>
#endif
#include "cuda_pointer.h"
#include "cuda_share.h"
#define NTHREAD 128


struct Particle{
	float2 pos[3];
	float mass;
	float pad;

	Particle(double x[3], double m){
		pos[0] = float2_split(x[0]);
		pos[1] = float2_split(x[1]);
		pos[2] = float2_split(x[2]);
		mass = (float)m;
	}
	Particle(int){
		pos[0].x = pos[0].y = pos[1].x = pos[1].y = pos[2].x = pos[2].y = mass = pad = 0.f;
	}
	__device__ Particle() {}
};

__global__ void pot_kernel(int n, int istart, Particle *ptcl, float2 *phi){
	__shared__ Particle jpbuf[NTHREAD];
	int i = NTHREAD * blockIdx.x + threadIdx.x;
	Particle ip = ptcl[i+istart-1];
	float2 phii = make_float2(0.f, 0.f);
	for(int j=0; j<n; j+= NTHREAD){
		__syncthreads();
		jpbuf[threadIdx.x] = ptcl[j + threadIdx.x];
		__syncthreads();
#pragma unroll 4
		for(int jj=0; jj<NTHREAD; jj++){
			// if(j+jj == i) continue;
			Particle &jp = jpbuf[jj];
			float dx = (jp.pos[0].x - ip.pos[0].x) + (jp.pos[0].y - ip.pos[0].y);
			float dy = (jp.pos[1].x - ip.pos[1].x) + (jp.pos[1].y - ip.pos[1].y);
			float dz = (jp.pos[2].x - ip.pos[2].x) + (jp.pos[2].y - ip.pos[2].y);
			float r2 = dx*dx + dy*dy + dz*dz;
			// if(r2==0.f) continue;
			float pij = jp.mass * rsqrtf(r2);
			// phii = float2_accum(phii, pij);
			if(r2 > 0.f) phii = float2_accum(phii, pij);
		}
		phii = float2_regularize(phii);
	}
	phi[i] = phii;
}

extern "C"  void gpunb_devinit_(int *irank);

void gpupot(
        int *rank,
        int istart,
        int ni,
		int n,
		double m[],
		double x[][3],
		double pot[]){
	gpunb_devinit_(rank);

	double t0 = get_wtime();
	cudaPointer <float2> phi;
	cudaPointer <Particle> ptcl;
	int ng = NTHREAD * (ni/NTHREAD + (ni%NTHREAD ? 1 : 0));
	int ntg = NTHREAD * (n/NTHREAD + (n%NTHREAD ? 1 : 0));
	if (ntg<ng+istart) ntg += NTHREAD;

	phi.allocate(ng);
	ptcl.allocate(ntg);

    //    std::cout << n << " " << ng << " "<< ntg << std::endl;
	for(int i=0; i<n; i++){
		// ptcl_h[i] = Particle(x[i], m[i]);
		ptcl[i] = Particle(x[i], m[i]);
	}
	for(int i=n; i<ntg; i++){
		// ptcl_h[i] = Particle(0);
		ptcl[i] = Particle(0);
	}

	// cudaMemcpy(ptcl_d, ptcl_h, ng * sizeof(Particle), cudaMemcpyHostToDevice);
	ptcl.htod(ntg);
	
	dim3 grid(ng/NTHREAD, 1, 1);
	dim3 threads(NTHREAD, 1, 1);
	int sharedMemSize = NTHREAD * sizeof(Particle);
	// pot_kernel <<<grid, threads, sharedMemSize >>> (n, ptcl_d, phi_d);
	pot_kernel <<<grid, threads, sharedMemSize >>> (n, istart, ptcl, phi);

	// cudaMemcpy(phi_h, phi_d, n * sizeof(float2), cudaMemcpyDeviceToHost);
	phi.dtoh(ni);
	for(int i=0; i<ni; i++){
		// pot[i] = (double)phi_h[i].x + (double)phi_h[i].y;
		pot[i] = (double)phi[i].x + (double)phi[i].y;
	}

	phi.free();
	ptcl.free();
	double t1 = get_wtime();
#ifdef PROFILE
//  R: rank; Ni: input i particle; NTOT: total j particle; pot: calculation time
	fprintf(stderr, "[R.%d GPU Pot.A] Ni %d  NTOT %d  pot(s) %f\n",*rank,ni,n,t1 - t0);
#endif
}

extern "C"{
	void gpupot_(
            int *irank,
            int *istart,
            int *ni,
			int *n,
			double m[],
			double x[][3],
			double pot[]){
      gpupot(irank, *istart, *ni, *n, m, x, pot);
	}
}
