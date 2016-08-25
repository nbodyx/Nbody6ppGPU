#include <cmath>
#include <stdio.h>
#include <unistd.h>
#include <cassert>
#include <cstdlib>
#include <omp.h>
#ifdef CUDA_5
#  include <helper_cuda.h>
#  define CUDA_SAFE_CALL checkCudaErrors
#else
#  include <cutil.h>
#endif
#include "cuda_pointer.h"
//#include "cuda_share.h"

#ifndef POTENTIAL
#define POTENTIAL
#endif

#define NTHREAD 64 // 64 or 128
// #define NJBLOCK 14 // for GTX 470
#define NJBLOCK 28 // for GTX660Ti 
#define NIBLOCK 32 // 16 or 32 
#define NIMAX (NTHREAD * NIBLOCK) // 2048

#define NXREDUCE 32 // must be >NJBLOCK
#define NYREDUCE  8

#define NNB_PER_BLOCK 256 // NNB per block, must be power of 2
#define NB_BUF_SIZE (1<<21)
// #define NNB_MAX       384 // total NNB at reduced

#define MAX_CPU 40
#define MAX_GPU 4

// for clearity, for myself
#define __out

#define NAN_CHECK(val) assert((val) == (val));

#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#else
static double get_wtime(){
	return 0.0;
}
#endif

typedef unsigned short uint16;

struct Jparticle{
	float3 pos;
	float  mass;
	float3 vel;
	float  pad;
	Jparticle() {}
	Jparticle(double mj, double xj[3], double vj[3]){
		pos.x = xj[0];
		pos.y = xj[1];
		pos.z = xj[2];
		mass  = mj;
		vel.x = vj[0];
		vel.y = vj[1];
		vel.z = vj[2];

		NAN_CHECK(xj[0]);
		NAN_CHECK(xj[1]);
		NAN_CHECK(xj[2]);
		NAN_CHECK(mj);
		NAN_CHECK(vj[0]);
		NAN_CHECK(vj[1]);
		NAN_CHECK(vj[2]);
	}
	__device__
	Jparticle(const float4 *buf){
		float4 tmp1 = buf[0];
		float4 tmp2 = buf[1];
		pos.x = tmp1.x;
		pos.y = tmp1.y;
		pos.z = tmp1.z;
		mass  = tmp1.w;
		vel.x = tmp2.x;
		vel.y = tmp2.y;
		vel.z = tmp2.z;
	}
};
struct Iparticle{
	float3 pos;
	float  h2;
	float3 vel;
	float  dtr;
	Iparticle() {}
	Iparticle(double h2i, double dtri, double xi[3], double vi[3]){
		pos.x = xi[0];
		pos.y = xi[1];
		pos.z = xi[2];
		h2    = h2i;
		vel.x = vi[0];
		vel.y = vi[1];
		vel.z = vi[2];
		dtr   = dtri;

		NAN_CHECK(xi[0]);
		NAN_CHECK(xi[1]);
		NAN_CHECK(xi[2]);
		NAN_CHECK(h2i);
		NAN_CHECK(vi[0]);
		NAN_CHECK(vi[1]);
		NAN_CHECK(vi[2]);
	}
};
struct Force{
	float3 acc;
	float  pot;
	float3 jrk;
	int    nnb;          //  8 words
	__device__  void clear(){
		acc.x = acc.y = acc.z = 0.f;
		jrk.x = jrk.y = jrk.z = 0.f;
		pot = 0.f;
		nnb = 0;
	}
	__device__ void operator+=(const Force &rhs){
		acc.x += rhs.acc.x;
		acc.y += rhs.acc.y;
		acc.z += rhs.acc.z;
		pot   += rhs.pot;
		jrk.x += rhs.jrk.x;
		jrk.y += rhs.jrk.y;
		jrk.z += rhs.jrk.z;
		if(nnb>=0 && rhs.nnb>=0){
			nnb += rhs.nnb;
		}else{
          nnb = -abs(nnb) - abs(rhs.nnb);
		}
	}
#if __CUDA_ARCH__ >= 300
	__device__ void reduce_with(const int mask){
		acc.x += __shfl_xor(acc.x, mask);
		acc.y += __shfl_xor(acc.y, mask);
		acc.z += __shfl_xor(acc.z, mask);
		pot   += __shfl_xor(pot  , mask);
		jrk.x += __shfl_xor(jrk.x, mask);
		jrk.y += __shfl_xor(jrk.y, mask);
		jrk.z += __shfl_xor(jrk.z, mask);
		int ntmp = __shfl_xor(nnb, mask);
		if(nnb>=0 && ntmp>=0){
			nnb += ntmp;
		}else{
          nnb = -abs(nnb) - abs(ntmp);
		}
	}
#endif
};

__device__ void dev_gravity(
		const int        jidx,
		const Iparticle &ip, 
		const Jparticle &jp, 
		__out Force     &fo,
		__out uint16     nblist[]){
	float dx = jp.pos.x - ip.pos.x;
	float dy = jp.pos.y - ip.pos.y;
	float dz = jp.pos.z - ip.pos.z;
	float dvx = jp.vel.x - ip.vel.x;
	float dvy = jp.vel.y - ip.vel.y;
	float dvz = jp.vel.z - ip.vel.z;

	float r2 = dx*dx + dy*dy + dz*dz;
#if 1
	float dxp = dx + ip.dtr * dvx;
	float dyp = dy + ip.dtr * dvy;
	float dzp = dz + ip.dtr * dvz;
	float r2p = dxp*dxp + dyp*dyp + dzp*dzp;
#else
	float r2p = r2;
#endif
	float rv = dx*dvx + dy*dvy + dz*dvz;
	float rinv1 = rsqrtf(r2);
    //	if(min(r2, r2p)  < jp.mass * ip.h2){
	if(min(r2, r2p)  < ip.h2){
		// fo.neib[fo.nnb++ % NBMAX] = j;
		nblist[fo.nnb & (NNB_PER_BLOCK-1)] = (uint16)jidx;
		fo.nnb++;
		rinv1 = 0.f;
	}
	float rinv2 = rinv1 * rinv1;
	float mrinv1 = jp.mass * rinv1;
	float mrinv3 = mrinv1 * rinv2;
	rv *= -3.f * rinv2;
	
#ifdef POTENTIAL
	fo.pot += mrinv1;
#endif
	fo.acc.x += mrinv3 * dx;
	fo.acc.y += mrinv3 * dy;
	fo.acc.z += mrinv3 * dz;
	// fo.acc.z += 1.0;
	fo.jrk.x += mrinv3 * (dvx + rv * dx);
	fo.jrk.y += mrinv3 * (dvy + rv * dy);
	fo.jrk.z += mrinv3 * (dvz + rv * dz);
}

//Use mass ratio as the factor to determine neighbor membership
__device__ void dev_gravity_m(
		const int        jidx,
		const Iparticle &ip, 
		const Jparticle &jp, 
		__out Force     &fo,
		__out uint16     nblist[]){
	float dx = jp.pos.x - ip.pos.x;
	float dy = jp.pos.y - ip.pos.y;
	float dz = jp.pos.z - ip.pos.z;
	float dvx = jp.vel.x - ip.vel.x;
	float dvy = jp.vel.y - ip.vel.y;
	float dvz = jp.vel.z - ip.vel.z;

	float r2 = dx*dx + dy*dy + dz*dz;
#if 1
	float dxp = dx + ip.dtr * dvx;
	float dyp = dy + ip.dtr * dvy;
	float dzp = dz + ip.dtr * dvz;
	float r2p = dxp*dxp + dyp*dyp + dzp*dzp;
#else
	float r2p = r2;
#endif
	float rv = dx*dvx + dy*dvy + dz*dvz;
	float rinv1 = rsqrtf(r2);
    if(min(r2, r2p)  < jp.mass * ip.h2){
    //	if(min(r2, r2p)  < ip.h2){
		// fo.neib[fo.nnb++ % NBMAX] = j;
		nblist[fo.nnb & (NNB_PER_BLOCK-1)] = (uint16)jidx;
		fo.nnb++;
		rinv1 = 0.f;
	}
	float rinv2 = rinv1 * rinv1;
	float mrinv1 = jp.mass * rinv1;
	float mrinv3 = mrinv1 * rinv2;
	rv *= -3.f * rinv2;
	
#ifdef POTENTIAL
	fo.pot += mrinv1;
#endif
	fo.acc.x += mrinv3 * dx;
	fo.acc.y += mrinv3 * dy;
	fo.acc.z += mrinv3 * dz;
	// fo.acc.z += 1.0;
	fo.jrk.x += mrinv3 * (dvx + rv * dx);
	fo.jrk.y += mrinv3 * (dvy + rv * dy);
	fo.jrk.z += mrinv3 * (dvz + rv * dz);
}

__global__ void gravity_kernel(
		const int       nbody,
		const Iparticle ipbuf[],
		const Jparticle jpbuf[],
		__out Force     fobuf[][NJBLOCK],
		__out uint16    nbbuf[][NJBLOCK][NNB_PER_BLOCK]){
	int ibid = blockIdx.x;
	int jbid = blockIdx.y;
	int tid = threadIdx.x;
	int iaddr = tid + blockDim.x * ibid;
	int jstart = (nbody * (jbid  )) / NJBLOCK;
	int jend   = (nbody * (jbid+1)) / NJBLOCK;

	Iparticle ip = ipbuf[iaddr];
	Force fo;
	fo.clear();
	uint16 *nblist = nbbuf[iaddr][jbid];
#if __CUDA_ARCH__ >= 300 // just some trial
	for(int j=jstart; j<jend; j+=32){
		__shared__ Jparticle jpshare[32];
		__syncthreads();
		float4 *src = (float4 *)&jpbuf[j];
		float4 *dst = (float4 *)jpshare;
		dst[tid] = src[tid];
		__syncthreads();
		if(jend-j < 32){
#pragma unroll 4
			for(int jj=0; jj<jend-j; jj++){
				const Jparticle jp = jpshare[jj];
				// const Jparticle jp( (float4 *)jpshare + 2*jj);
				dev_gravity(j-jstart+jj, ip, jp, fo, nblist);
			}
		}else{
#pragma unroll 8
			for(int jj=0; jj<32; jj++){
				const Jparticle jp = jpshare[jj];
				// const Jparticle jp( (float4 *)jpshare + 2*jj);
				dev_gravity(j-jstart+jj, ip, jp, fo, nblist);
			}
		}
	}
#else
	for(int j=jstart; j<jend; j+=NTHREAD){
		__shared__ Jparticle jpshare[NTHREAD];
		__syncthreads();
		float4 *src = (float4 *)&jpbuf[j];
		float4 *dst = (float4 *)jpshare;
		dst[        tid] = src[        tid];
		dst[NTHREAD+tid] = src[NTHREAD+tid];
		__syncthreads();

		if(jend-j < NTHREAD){
#pragma unroll 4
			for(int jj=0; jj<jend-j; jj++){
				Jparticle jp = jpshare[jj];
				dev_gravity(j-jstart+jj, ip, jp, fo, nblist);
			}
		}else{
#pragma unroll 8
			for(int jj=0; jj<NTHREAD; jj++){
				Jparticle jp = jpshare[jj];
				dev_gravity(j-jstart+jj, ip, jp, fo, nblist);
			}
		}
	}
#endif
	if(fo.nnb > NNB_PER_BLOCK) fo.nnb = -fo.nnb;
	fobuf[iaddr][jbid] = fo;
}

//Use mass ratio as neighbor radius criterion
__global__ void gravity_kernel_m(
		const int       nbody,
		const Iparticle ipbuf[],
		const Jparticle jpbuf[],
		__out Force     fobuf[][NJBLOCK],
		__out uint16    nbbuf[][NJBLOCK][NNB_PER_BLOCK]){
	int ibid = blockIdx.x;
	int jbid = blockIdx.y;
	int tid = threadIdx.x;
	int iaddr = tid + blockDim.x * ibid;
	int jstart = (nbody * (jbid  )) / NJBLOCK;
	int jend   = (nbody * (jbid+1)) / NJBLOCK;

	Iparticle ip = ipbuf[iaddr];
	Force fo;
	fo.clear();
	uint16 *nblist = nbbuf[iaddr][jbid];
#if __CUDA_ARCH__ >= 300 // just some trial
	for(int j=jstart; j<jend; j+=32){
		__shared__ Jparticle jpshare[32];
		__syncthreads();
		float4 *src = (float4 *)&jpbuf[j];
		float4 *dst = (float4 *)jpshare;
		dst[tid] = src[tid];
		__syncthreads();
		if(jend-j < 32){
#pragma unroll 4
			for(int jj=0; jj<jend-j; jj++){
				const Jparticle jp = jpshare[jj];
				// const Jparticle jp( (float4 *)jpshare + 2*jj);
				dev_gravity_m(j-jstart+jj, ip, jp, fo, nblist);
			}
		}else{
#pragma unroll 8
			for(int jj=0; jj<32; jj++){
				const Jparticle jp = jpshare[jj];
				// const Jparticle jp( (float4 *)jpshare + 2*jj);
				dev_gravity_m(j-jstart+jj, ip, jp, fo, nblist);
			}
		}
	}
#else
	for(int j=jstart; j<jend; j+=NTHREAD){
		__shared__ Jparticle jpshare[NTHREAD];
		__syncthreads();
		float4 *src = (float4 *)&jpbuf[j];
		float4 *dst = (float4 *)jpshare;
		dst[        tid] = src[        tid];
		dst[NTHREAD+tid] = src[NTHREAD+tid];
		__syncthreads();

		if(jend-j < NTHREAD){
#pragma unroll 4
			for(int jj=0; jj<jend-j; jj++){
				Jparticle jp = jpshare[jj];
				dev_gravity_m(j-jstart+jj, ip, jp, fo, nblist);
			}
		}else{
#pragma unroll 8
			for(int jj=0; jj<NTHREAD; jj++){
				Jparticle jp = jpshare[jj];
				dev_gravity_m(j-jstart+jj, ip, jp, fo, nblist);
			}
		}
	}
#endif
	if(fo.nnb > NNB_PER_BLOCK) fo.nnb = -fo.nnb;
	fobuf[iaddr][jbid] = fo;
}

#if __CUDA_ARCH__ >= 300
__device__ void warp_reduce_int(int inp, int *out){
	inp += __shfl_xor(inp, 1);
	inp += __shfl_xor(inp, 2);
	inp += __shfl_xor(inp, 4);
	inp += __shfl_xor(inp, 8);
# if NXREDUCE==32
	inp += __shfl_xor(inp, 16);
# endif
	*out = inp;
}
__device__ void warp_reduce_float8(float4 inp1, float4 inp2, float *out){
	const int tid = threadIdx.x;
	float4 tmp4L = (4&tid) ? inp2 : inp1;
	float4 tmp4R = (4&tid) ? inp1 : inp2;
	tmp4L.x += __shfl_xor(tmp4R.x, 4);
	tmp4L.y += __shfl_xor(tmp4R.y, 4);
	tmp4L.z += __shfl_xor(tmp4R.z, 4);
	tmp4L.w += __shfl_xor(tmp4R.w, 4);
	float4 tmp4;
	tmp4.x = (2&tid) ? tmp4L.z : tmp4L.x;
	tmp4.y = (2&tid) ? tmp4L.w : tmp4L.y;
	tmp4.z = (2&tid) ? tmp4L.x : tmp4L.z;
	tmp4.w = (2&tid) ? tmp4L.y : tmp4L.w;
	tmp4.x += __shfl_xor(tmp4.z, 2);
	tmp4.y += __shfl_xor(tmp4.w, 2);
	float2 tmp2;
	tmp2.x = (1&tid) ? tmp4.y : tmp4.x;
	tmp2.y = (1&tid) ? tmp4.x : tmp4.y;
	tmp2.x += __shfl_xor(tmp2.y, 1);

	tmp2.x += __shfl_xor(tmp2.x, 8);
# if NXREDUCE==32
	tmp2.x += __shfl_xor(tmp2.x, 16);
# endif
	if(tid < 8){
		out[tid] = tmp2.x;
	}
}
#endif

__global__ void force_reduce_kernel(
		const int ni,
		const Force fpart[][NJBLOCK],
		__out Force ftot []){
	const int xid = threadIdx.x;
	const int yid = threadIdx.y;
	const int bid = blockIdx.x;
	const int iaddr = yid + blockDim.y * bid;

#if __CUDA_ARCH__ >= 300
	Force f;
	if(xid < NJBLOCK){
		f = fpart[iaddr][xid];
	}else{
		f.clear();
	}
# if 0
# pragma unroll
	for(int mask=1; mask<NXREDUCE; mask*=2){
		f.reduce_with(mask);
	}
	if(iaddr < ni && xid == 0){
		ftot[iaddr] = f;
	}
# else
	if(iaddr < ni){
		const float4 tmp1 = make_float4(f.acc.x, f.acc.y, f.acc.z, f.pot);
		const float4 tmp2 = make_float4(f.jrk.x, f.jrk.y, f.jrk.z, 0.0f);
		const int    itmp = f.nnb;
		float *dst  = (float *)(ftot + iaddr);
		int   *idst = (int *)(dst + 7);
		warp_reduce_float8(tmp1, tmp2, dst);
		warp_reduce_int(itmp, idst);
	}
# endif
#else
	__shared__ Force fshare[NYREDUCE][NXREDUCE];
	if(xid < NJBLOCK){
		fshare[yid][xid] = fpart[iaddr][xid];
	}else{
		fshare[yid][xid].clear();
	}
	Force *fs = fshare[yid];
#if NXREDUCE==32
	if(xid < 16) fs[xid] += fs[xid + 16];
#endif
	if(xid < 8) fs[xid] += fs[xid + 8];
	if(xid < 4) fs[xid] += fs[xid + 4];
	if(xid < 2) fs[xid] += fs[xid + 2];
	if(xid < 1) fs[xid] += fs[xid + 1];
	
	if(iaddr < ni){
		ftot[iaddr] = fs[0];
	}
#endif
}

__global__ void gather_nb_kernel(
		const int    ni,
		const int    nj,
		const int    joff,
		const Force  fpart[][NJBLOCK],
		const Force  ftot [],
		const int    nboff[],
		const uint16 nbpart[][NJBLOCK][NNB_PER_BLOCK],
		__out   int  nblist[])
{
	const int xid = threadIdx.x;
	const int yid = threadIdx.y;
	const int bid = blockIdx.x;
	const int iaddr = yid + blockDim.y * bid;
	if(iaddr >= ni) return;
	if(ftot[iaddr].nnb < 0) return;

	const int mynnb = (xid < NJBLOCK) ? fpart[iaddr][xid].nnb
	                                  : 0;

	// now performe prefix sum
#if __CUDA_ARCH__ >= 300
	int ix = mynnb;
#pragma unroll
	for(int ioff=1; ioff<NXREDUCE; ioff*=2){
		int iy = __shfl_up(ix, ioff);
		if(xid>=ioff) ix += iy;
	}
	int iz = __shfl_up(ix, 1);
	const int off = (xid == 0) ? 0 : iz;
#else
	__shared__ int ishare[NYREDUCE][NXREDUCE];
	ishare[yid][xid] = mynnb;
	volatile int *ish = ishare[yid];
	if(xid>=1)  ish[xid] += ish[xid-1];
	if(xid>=2)  ish[xid] += ish[xid-2];
	if(xid>=4)  ish[xid] += ish[xid-4];
	if(xid>=8)  ish[xid] += ish[xid-8];
#if NXREDUCE==32
	if(xid>=16)  ish[xid] += ish[xid-16];
#endif
	const int off = (xid == 0) ? 0 
	                           : ish[xid-1];
#endif
	int *nbdst = nblist + nboff[iaddr] + off;

	const int jstart = (nj * xid) / NJBLOCK;
	if(xid < NJBLOCK){
		for(int k=0; k<mynnb; k++){
			const int nbid = (joff + jstart) + int(nbpart[iaddr][xid][k]);
			// const int nbid = iaddr * 1000 + k;
			nbdst[k] = nbid;
		}
	}
}


// Host Part

static double time_send, time_grav, time_reduce, time_out;
static long long numInter;
static int icall,ini,isend;
static cudaPointer <Jparticle> jpbuf[MAX_GPU];
static cudaPointer <Iparticle> ipbuf[MAX_GPU];
static cudaPointer <Force[NJBLOCK]> fpart[MAX_GPU];
static cudaPointer <Force>          ftot [MAX_GPU];
static cudaPointer <uint16[NJBLOCK][NNB_PER_BLOCK]> nbpart[MAX_GPU];
static cudaPointer <int> nblist [MAX_GPU];
static cudaPointer <int> nboff  [MAX_GPU];
static int numCPU, numGPU;
static int joff[MAX_GPU + 1];
static int nbody, nbodymax;
static int devid[MAX_GPU];
static bool is_open = false;
static bool devinit = false;

void GPUNB_devinit(int irank){
	if(devinit) return;

	assert(NXREDUCE >= NJBLOCK);
	assert(NXREDUCE <= 32);

	cudaGetDeviceCount(&numGPU);
	assert(numGPU <= MAX_GPU);
	char *gpu_list = getenv("GPU_LIST");
	if(gpu_list){
		// get GPU list from environment variable
		numGPU = 0;
		char *p = strtok(gpu_list, " ");
		while(p){
			devid[numGPU++] = atoi(p);
			p = strtok(NULL, " ");
			assert(numGPU <= MAX_GPU);
		}
	}else{
		// use all GPUs
		for(int i=0; i<numGPU; i++){
			devid[i] = i;
		}
	}
	
	// numGPU = 1;
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid == 0) numCPU = omp_get_num_threads();
	}
	assert(numCPU <= MAX_CPU);
	assert(numGPU <= numCPU);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid < numGPU){
			cudaSetDevice(devid[tid]);
		}
	}

#ifdef PROFILE
    char hostname[150];
    memset(hostname,0,150);
    gethostname(hostname,150);
	for(int i=0; i<numGPU; i++){
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, devid[i]);
      fprintf(stderr, "# GPU initialization - rank: %d; HOST %s; NGPU %d; device: %d %s\n", irank, hostname, numGPU, devid[i], prop.name);
	}
#endif
    devinit = true;
}

void GPUNB_open(int nbmax,int irank){
	time_send = time_grav = time_reduce = time_out = 0.0;
	numInter = 0;
    icall = ini = isend = 0;
	nbodymax = nbmax;

	GPUNB_devinit(irank);

	if(is_open){
		fprintf(stderr, "gpunb: it is already open\n");
		return;
	}
	is_open = true;


	for(int id=0; id<numGPU + 1; id++){
		joff[id] = (id * nbmax) / numGPU;
	}

	// omp_set_num_threads(numGPU);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid < numGPU){
			cudaSetDevice(devid[tid]);
			int nj = joff[tid+1] - joff[tid];
			jpbuf [tid].allocate(nj + NTHREAD);
			ipbuf [tid].allocate(NIMAX);
			fpart [tid].allocate(NIMAX);
			ftot  [tid].allocate(NIMAX);
			nbpart[tid].allocate(NIMAX);
			nblist[tid].allocate(NB_BUF_SIZE); // total ganged nblist
			nboff [tid].allocate(NIMAX+1);
		}
	}
#ifdef PROFILE
    fprintf(stderr, "# Open GPU regular force - rank: %d; nbmax: %d\n", irank, nbmax);
#endif
}

void GPUNB_close(){
	if(!is_open){
		fprintf(stderr, "gpunb: it is already close\n");
		return;
	}
	is_open = false;
	// omp_set_num_threads(numGPU);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid < numGPU){
			jpbuf [tid].free();
			ipbuf [tid].free();
			fpart [tid].free();
			ftot  [tid].free();
			nbpart[tid].free();
			nblist[tid].free();
			nboff [tid].free();
		}
	}
	// omp_set_num_threads(numCPU);
	nbodymax = 0;

// #ifdef PROFILE
// 	fprintf(stderr, "***********************\n");
// 	fprintf(stderr, "Closed NBODY6/GPU library\n");
// 	fprintf(stderr, "time send   : %f sec\n", time_send);
// 	fprintf(stderr, "time grav   : %f sec\n", time_grav);
// 	fprintf(stderr, "time reduce : %f sec\n", time_reduce);
// 	fprintf(stderr, "time regtot : %f sec\n", time_send + time_grav + time_reduce);
// 	fprintf(stderr, "%f Gflops (gravity part only)\n", 60.e-9 * numInter / time_grav);
// 	fprintf(stderr, "***********************\n");
// #endif
}

void GPUNB_send(
		int _nbody,
		double mj[],
		double xj[][3],
		double vj[][3]){
	assert(is_open);
	nbody = _nbody;
	assert(nbody <= nbodymax);
	time_send -= get_wtime();
    isend++;
	for(int id=0; id<numGPU + 1; id++){
		joff[id] = (id * nbody) / numGPU;
	}
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid < numGPU){
			int nj = joff[tid+1] - joff[tid];
			for(int j=0; j<nj; j++){
				int jj = j + joff[tid];
				jpbuf[tid][j] = Jparticle(mj[jj], xj[jj], vj[jj]);
			}
			jpbuf[tid].htod(nj);
		}
	}
	time_send += get_wtime();
}

void GPUNB_regf(
		int ni,
		double h2[],
		double dtr[],
		double xi[][3],
		double vi[][3],
		double acc[][3],
		double jrk[][3],
		double pot[],
		int lmax,
		int nnbmax,
		int *listbase,
        int m_flag){
	assert(is_open);

	time_grav -= get_wtime();
	numInter += ni * nbody;
    ini +=ni;
    icall++;
	assert(0 < ni && ni <= NIMAX);

	// omp_set_num_threads(numGPU);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid < numGPU){
			// cudaSetDevice(device_id[tid]);
			for(int i=0; i<ni; i++){
				ipbuf[tid][i] = Iparticle(h2[i], dtr[i], xi[i], vi[i]);
			}
			// set i-particles
			ipbuf[tid].htod(ni);

			// gravity kernel
			int niblock = 1 + (ni-1) / NTHREAD;
			dim3 grid(niblock, NJBLOCK, 1);
			dim3 threads(NTHREAD, 1, 1);
			int nj = joff[tid+1] - joff[tid];
            if(m_flag) {
              gravity_kernel_m <<< grid, threads >>> 
				(nj, ipbuf[tid], jpbuf[tid], fpart[tid], nbpart[tid]);
            }
            else {
              gravity_kernel <<< grid, threads >>> 
				(nj, ipbuf[tid], jpbuf[tid], fpart[tid], nbpart[tid]);
            }
			// CUDA_SAFE_THREAD_SYNC();

#if 0
			dim3 rgrid(niblock, 1, 1);
			reduce_kernel <<< rgrid, threads >>>
				(nj, joff[tid], fpart[tid], nbpart[tid], ftot[tid], nbtot[tid]);
#else
			const int ni8 = 1 + (ni-1) / NYREDUCE;
			dim3 rgrid   (ni8, 1, 1);
			dim3 rthreads(NXREDUCE, NYREDUCE, 1);
			force_reduce_kernel <<< rgrid, rthreads >>>
				(ni, fpart[tid], ftot[tid]);
#endif
			// CUDA_SAFE_THREAD_SYNC();
			ftot [tid].dtoh(ni);

			// now make prefix sum
			int nbsum = 0;
			for(int i=0; i<ni; i++){
				nboff[tid][i] = nbsum;
				const int nnb = ftot[tid][i].nnb;
				// assert(nnb >= 0);
				if(nnb >= 0) nbsum += nnb;
			}
			assert(nbsum <= NB_BUF_SIZE);
			nboff[tid].htod(ni);

			// debugging
			// for(int k=0; k<nbsum; k++) nblist[tid][k] = -1;
			// nblist[tid].htod(nbsum);

			gather_nb_kernel <<< rgrid, rthreads>>>
				(ni, nj, joff[tid], fpart[tid], ftot[tid], 
				 nboff[tid], nbpart[tid], nblist[tid]);
			// CUDA_SAFE_THREAD_SYNC();
			nblist[tid].dtoh(nbsum);
		}
	}

	const double wt = get_wtime();
	time_grav   += wt;
	time_reduce -= wt;

	// reduction phase
	// omp_set_num_threads(numCPU);
#pragma omp parallel for
	for(int i=0; i<ni; i++){
		double ax=0.0, ay=0.0, az=0.0;
		double jx=0.0, jy=0.0, jz=0.0;
		double po=0.0;

		for(int id=0; id<numGPU; id++){
			Force &fo = ftot[id][i];
			ax += fo.acc.x;
			ay += fo.acc.y;
			az += fo.acc.z;
			jx += fo.jrk.x;
			jy += fo.jrk.y;
			jz += fo.jrk.z;
			po += fo.pot;
		}
		acc[i][0] = ax;
		acc[i][1] = ay;
		acc[i][2] = az;
		jrk[i][0] = jx;
		jrk[i][1] = jy;
		jrk[i][2] = jz;
		pot[i]    = po;
	}
#pragma omp parallel for
	for(int i=0; i<ni; i++){
		bool overflow = false;
		int *nnbp = listbase + lmax * i;
		int *nblistp = nnbp + 1;
		int nnb = 0;
		for(int id=0; id<numGPU; id++){
			const int nnb_part = ftot[id][i].nnb;
			if(nnb_part < 0){
				overflow = true;
//				fprintf(stderr, "!!!overflow : i=%d, id=%d, nnb_part=%d\n", i, id, nnb_part);
			}
			// assert(!overflow);
			nnb += nnb_part;
			if(nnb > nnbmax){
				overflow = true;
//				fprintf(stderr, "!!!overflow : i=%d, id=%d, nnb_tot =%d, nnbmax=%d\n", i, id, nnb, nnbmax);
			}
			// assert(!overflow);
			if(!overflow){
				const int off = nboff[id][i]; 
				for(int k=0; k<nnb_part; k++){
					*nblistp++ = nblist[id][off + k];
				}
			}
		}
		if(overflow){
			// *nnbp = -1;
			*nnbp = nnb ? -abs(nnb) : -9999;
		}else{
			*nnbp = nnb;
		}
	}
	time_reduce += get_wtime();
}

void GPUNB_profile(int irank) {
#ifdef PROFILE
  if(icall) {
    // R: rank; D: GPU device number;
    // Nsend: number of call gpunb_send between two profile check points (adjust time interval);
    // Ngrav: number of call gpunb_ref between two profile check points;
    // <Ni>:  averaged i particle number per regular block;
    // send: j particle sending time;
    // grav: force calculation time;
    // nb:   gather neighbor and force time;
    // out:  output data time;
    // Perf: performance for gpu regular calculation
    fprintf(stderr,"[R.%d GPU Reg.F ] Nsend %d  Ngrav %d  <Ni> %d   send(s) %f grav(s) %f  nb(s) %f  out(s) %f  Perf.(Gflops) %f\n",irank,isend,icall,ini/isend,time_send,time_grav,time_reduce,time_out,60.e-9*numInter/time_grav);
  }
  time_send = time_grav = time_reduce = time_out = 0.0;
  numInter = 0;
  icall = ini = isend= 0;
#else
  return;
#endif
}

extern "C" {
	void gpunb_devinit_(int *irank){
		GPUNB_devinit(*irank);
	}
	void gpunb_open_(int *nbmax,int *irank){
		GPUNB_open(*nbmax,*irank);
	}
	void gpunb_close_(){
		GPUNB_close();
	}
	void gpunb_send_(
			int *nj,
			double mj[],
			double xj[][3],
			double vj[][3]){
		GPUNB_send(*nj, mj, xj, vj);
	}
	void gpunb_regf_(
			int *ni,
			double h2[],
			double dtr[],
			double xi[][3],
			double vi[][3],
			double acc[][3],
			double jrk[][3],
			double pot[],
			int *lmax,
			int *nbmax,
			int *list,
            int *m_flag ){ // list[][lmax]
         GPUNB_regf(*ni, h2, dtr, xi, vi, acc, jrk, pot, *lmax, *nbmax, list, *m_flag);
	}
    void gpunb_profile_(int *irank){
      GPUNB_profile(*irank);
    }
}

