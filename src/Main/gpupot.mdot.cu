//#include <iostream>
#include <cstdio>
// #include <cutil.h>
#include <cassert>
#ifdef CUDA_5
#  include <helper_cuda.h>
#  define CUDA_SAFE_CALL checkCudaErrors
#else
#  include <cutil.h>
#endif
#include "cuda_pointer.h"
#include "cuda_share.h"

#define NTHREAD 64 // 64, 96, 128 or 192; should be same as the one in gpunb.gpu.cu
#define NJBLOCK 28 // 8800GTS/512 has 16
#define NIBLOCK 32 // 16 or 32 
#define NIMAX (NTHREAD * NIBLOCK) // 2048

#define NXREDUCE 32 // must be >NJBLOCK
#define NYREDUCE 8

//#define NAN_CHECK(val) assert((val) == (val));

//from gpunb.gpu.cu=================================//

struct Particle{
	float2 pos[3];
	float mass;
	float pad;

	Particle(double x[3], double m){
		pos[0] = float2_split(x[0]);
		pos[1] = float2_split(x[1]);
		pos[2] = float2_split(x[2]);
		mass = (float)m;

        NAN_CHECK(x[0]);
        NAN_CHECK(x[1]);
        NAN_CHECK(x[2]);
        NAN_CHECK(m);
	}
	Particle(int){
		pos[0].x = pos[0].y = pos[1].x = pos[1].y = pos[2].x = pos[2].y = mass = pad = 0.f;
	}
	__device__ Particle() {}
};

__global__ void pot_reduce_kernel(
		const int ni,
		const float2 phipart[][NJBLOCK],
        float2 phi[]){
  //thread x * y + block x============================//
  //thread x for NJBLOCK==============================//
	const int xid = threadIdx.x;
	const int yid = threadIdx.y;
	const int bid = blockIdx.x;
    //thread y & block x for active particle============//
	const int iaddr = yid + blockDim.y * bid;

	__shared__ float2 phishare[NYREDUCE][NXREDUCE];

    __syncthreads();
	if(xid < NJBLOCK){
      phishare[yid][xid] = phipart[iaddr][xid];
	}else{
      phishare[yid][xid] = make_float2(0.f,0.f);
	}
    __syncthreads();
	float2 *phis = phishare[yid];
    
#if NXREDUCE==32
	if(xid < 16) phis[xid] = float2_add(phis[xid],phis[xid + 16]);
#endif
	if(xid < 8) phis[xid] = float2_add(phis[xid],phis[xid + 8]);
	if(xid < 4) phis[xid] = float2_add(phis[xid],phis[xid + 4]);
	if(xid < 2) phis[xid] = float2_add(phis[xid],phis[xid + 2]);
	if(xid < 1) phis[xid] = float2_add(phis[xid],phis[xid + 1]);
	
	if(iaddr < ni){
      phi[iaddr] = float2_regularize(phis[0]);
	}
}

__global__ void pot_kernel_float(
          int ni,                           
          int n,
          Jparticle *ipbuf,
          Jparticle *ptcl,
          float2 phipart[][NJBLOCK]){
	int i = NTHREAD * blockIdx.x + threadIdx.x;
    int jbid = blockIdx.y;
    int jstart = (n * (jbid  )) / NJBLOCK;
    int jend   = (n * (jbid+1)) / NJBLOCK;

    Jparticle ip=Jparticle();
    if(i<ni) ip = ipbuf[i];
	float2 phii = make_float2(0.f, 0.f);
	for(int j=jstart; j<jend; j+= NTHREAD){
      __shared__ Jparticle jpbuf_f[NTHREAD];

      __syncthreads();
      float4 *src = (float4 *)&ptcl[j];
      float4 *dst = (float4 *)jpbuf_f;
      dst[threadIdx.x]         = src[threadIdx.x];
      dst[NTHREAD+threadIdx.x] = src[NTHREAD+threadIdx.x];
      __syncthreads();

      if(jend-j < NTHREAD) {
#pragma unroll 4        
        for(int jj=0; jj<jend-j; jj++){
			// if(j+jj == i) continue;
            Jparticle &jp = jpbuf_f[jj];
			float dx = jp.pos.x - ip.pos.x;
			float dy = jp.pos.y - ip.pos.y;
			float dz = jp.pos.z - ip.pos.z;
			float r2 = dx*dx + dy*dy + dz*dz;
			// if(r2==0.f) continue;
			float pij = jp.mass * rsqrtf(r2);
			// phii = float2_accum(phii, pij);
			if(r2 > 0.f) phii = float2_accum(phii, pij);
		}
      }else{
#pragma unroll 8
        for(int jj=0; jj<NTHREAD; jj++){
			// if(j+jj == i) continue;
			Jparticle &jp = jpbuf_f[jj];
			float dx = jp.pos.x - ip.pos.x;
			float dy = jp.pos.y - ip.pos.y;
			float dz = jp.pos.z - ip.pos.z;
			float r2 = dx*dx + dy*dy + dz*dz;
			// if(r2==0.f) continue;
			float pij = jp.mass * rsqrtf(r2);
			// phii = float2_accum(phii, pij);
			if(r2 > 0.f) phii = float2_accum(phii, pij);
		}
      }
      phii = float2_regularize(phii);
	}
#ifdef HIGH_CORR
    int istart = (ni* (jbid ))  / NJBLOCK;
    int iend   = (ni* (jbid+1)) / NJBLOCK;

    for(int j=istart; j<iend; j+= NTHREAD){
      __shared__ Jparticle jpbuf_f[NTHREAD];

      __syncthreads();
      float4 *src = (float4 *)&ipbuf[j];
      float4 *dst = (float4 *)jpbuf_f;
      dst[threadIdx.x]         = src[threadIdx.x];
      dst[NTHREAD+threadIdx.x] = src[NTHREAD+threadIdx.x];
      __syncthreads();

      if(iend-j < NTHREAD) {
#pragma unroll 4        
        for(int jj=0; jj<iend-j; jj++){
			// if(j+jj == i) continue;
            Jparticle &jp = jpbuf_f[jj];
			float dx = jp.pos.x - ip.pos.x;
			float dy = jp.pos.y - ip.pos.y;
			float dz = jp.pos.z - ip.pos.z;
			float r2 = dx*dx + dy*dy + dz*dz;
			// if(r2==0.f) continue;
			float pij = -0.5 * jp.mass * rsqrtf(r2);
			// phii = float2_accum(phii, pij);
			if(r2 > 0.f) phii = float2_accum(phii, pij);
		}
      }else{
#pragma unroll 8
        for(int jj=0; jj<NTHREAD; jj++){
			// if(j+jj == i) continue;
			Jparticle &jp = jpbuf_f[jj];
			float dx = jp.pos.x - ip.pos.x;
			float dy = jp.pos.y - ip.pos.y;
			float dz = jp.pos.z - ip.pos.z;
			float r2 = dx*dx + dy*dy + dz*dz;
			// if(r2==0.f) continue;
			float pij = -0.5 * jp.mass * rsqrtf(r2);
			// phii = float2_accum(phii, pij);
			if(r2 > 0.f) phii = float2_accum(phii, pij);
		}
      }
      phii = float2_regularize(phii);
	}
#endif
	phipart[i][jbid] = phii;
}

__global__ void pot_kernel(
          int ni,                           
          int n,
          Particle *ipbuf,
          Particle *ptcl,
          float2 phipart[][NJBLOCK]){
	int i = NTHREAD * blockIdx.x + threadIdx.x;
    int jbid = blockIdx.y;
    int jstart = (n * (jbid  )) / NJBLOCK;
    int jend   = (n * (jbid+1)) / NJBLOCK;

    Particle ip=Particle();
    if(i<ni) ip = ipbuf[i];
	float2 phii = make_float2(0.f, 0.f);
	for(int j=jstart; j<jend; j+= NTHREAD){
      __shared__ Particle jpbuf[NTHREAD];

      __syncthreads();
      float4 *src = (float4 *)&ptcl[j];
      float4 *dst = (float4 *)jpbuf;
      dst[threadIdx.x]         = src[threadIdx.x];
      dst[NTHREAD+threadIdx.x] = src[NTHREAD+threadIdx.x];
      __syncthreads();

      if(jend-j < NTHREAD) {
#pragma unroll 4        
        for(int jj=0; jj<jend-j; jj++){
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
      }else{
#pragma unroll 8
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
      }
      phii = float2_regularize(phii);
	}
#ifdef HIGH_CORR
    int istart = (ni* (jbid ))  / NJBLOCK;
    int iend   = (ni* (jbid+1)) / NJBLOCK;

    for(int j=istart; j<iend; j+= NTHREAD){
      __shared__ Particle jpbuf_f[NTHREAD];

      __syncthreads();
      float4 *src = (float4 *)&ipbuf[j];
      float4 *dst = (float4 *)jpbuf_f;
      dst[threadIdx.x]         = src[threadIdx.x];
      dst[NTHREAD+threadIdx.x] = src[NTHREAD+threadIdx.x];
      __syncthreads();

      if(iend-j < NTHREAD) {
#pragma unroll 4        
        for(int jj=0; jj<iend-j; jj++){
			// if(j+jj == i) continue;
            Particle &jp = jpbuf_f[jj];
			float dx = (jp.pos[0].x - ip.pos[0].x) + (jp.pos[0].y - ip.pos[0].y);
			float dy = (jp.pos[1].x - ip.pos[1].x) + (jp.pos[1].y - ip.pos[1].y);
			float dz = (jp.pos[2].x - ip.pos[2].x) + (jp.pos[2].y - ip.pos[2].y);
			float r2 = dx*dx + dy*dy + dz*dz;
			// if(r2==0.f) continue;
			float pij = 0.5 * jp.mass * rsqrtf(r2);
			// phii = float2_accum(phii, pij);
			if(r2 > 0.f) phii = float2_accum(phii, pij);
		}
      }else{
#pragma unroll 8
        for(int jj=0; jj<NTHREAD; jj++){
			// if(j+jj == i) continue;
			Particle &jp = jpbuf_f[jj];
			float dx = (jp.pos[0].x - ip.pos[0].x) + (jp.pos[0].y - ip.pos[0].y);
			float dy = (jp.pos[1].x - ip.pos[1].x) + (jp.pos[1].y - ip.pos[1].y);
			float dz = (jp.pos[2].x - ip.pos[2].x) + (jp.pos[2].y - ip.pos[2].y);
			float r2 = dx*dx + dy*dy + dz*dz;
			// if(r2==0.f) continue;
			float pij = 0.5 * jp.mass * rsqrtf(r2);
			// phii = float2_accum(phii, pij);
			if(r2 > 0.f) phii = float2_accum(phii, pij);
		}
      }
      phii = float2_regularize(phii);
	}
#endif
	phipart[i][jbid] = phii;
}

extern "C"  void gpunb_devinit_(int *irank);
//extern "C"  void gpunb_share_jp_(Jparticle *jpoint);
extern cudaPointer <Jparticle> jpbuf;
static cudaPointer <Jparticle> ipbuf;
static cudaPointer <float2> phi;
static cudaPointer <float2[NJBLOCK]> phipart;
static cudaPointer <Particle> ptcl;
static cudaPointer <Particle> ibuf;
static double tsend,tcalc;
static int icall,ini;

void init_ni() {
  phi.allocate(NIMAX);
  phipart.allocate(NIMAX);
}

void close_ni() {
  phi.free();
  phipart.free();
}  

void gpupot_init_float(int *irank){
  gpunb_devinit_(irank);
  ipbuf.allocate(NIMAX);
  init_ni();
  tsend=0.;
  tcalc=0.;
  icall=0;
  ini=0;
}

void gpupot_close_float(){
  close_ni();
  ipbuf.free();
}

void gpupot_init(int *irank, int n){
  int ntg = NTHREAD * (n/NTHREAD + (n%NTHREAD ? 1 : 0));
  gpunb_devinit_(irank);
  ibuf.allocate(NIMAX);
  init_ni();
  ptcl.allocate(ntg);
  tsend=0.;
  tcalc=0.;
  icall=0;
  ini=0;
}

void gpupot_close(){
  ibuf.free();
  ptcl.free();
  close_ni();
}

void gpupot_float(
        int *rank,
        int ni,
		int n,
        int ishift,
        int list[],
        double dm[],
		double pot[]){

	gpunb_devinit_(rank);
    assert(ni<=NIMAX);
    //    Jparticle *jpoint=NULL;
    // gpunb_share_jp_(jpoint);
    //    assert(jpoint);
    assert(jpbuf.size);
    
	tcalc -= get_wtime();
    icall++;
    ini +=ni;
    //    cudaPointer <int> plist;
    int ng = NTHREAD * (ni/NTHREAD + (ni%NTHREAD ? 1 : 0));

    // ipbuf.allocate(ng);
	// phi.allocate(ng);
    // phipart.allocate(ng);
    //    plist.allocate(ni);
    // int *plist;
    // CUDA_SAFE_CALL(cudaMalloc((void**)&plist, ni*sizeof(int)));
    // CUDA_SAFE_CALL(cudaMemcpy(plist,list, ni*sizeof(int), cudaMemcpyHostToDevice));
    // for(int i=0; i<ni; i++) {
    //   plist[i] = list[i] - 1;
    // }
    for(int i=0; i<ni; i++) {
      ipbuf[i] = jpbuf[list[i]-1+ishift];
      ipbuf[i].mass = dm[i];
    }
    for(int i=ni; i<ng; i++) {
      ipbuf[i] = Jparticle(0);
    }

    ipbuf.htod(ng);
    // plist.htod(ni);
	dim3 grid(ng/NTHREAD, NJBLOCK, 1);
	dim3 threads(NTHREAD, 1, 1);
	pot_kernel_float <<<grid, threads>>> (ni, n, ipbuf, jpbuf, phipart);

    const int ni8 = 1 + (ni-1) / NYREDUCE;
    dim3 rgrid (ni8, 1, 1);
    dim3 rthreads(NXREDUCE, NYREDUCE, 1);
    pot_reduce_kernel <<< rgrid, rthreads >>> (ni, phipart, phi);

    phi.dtoh(ni);

    for(int i=0; i<ni; i++){
      pot[i] = (double)phi[i].x + (double)phi[i].y;
	}

	tcalc += get_wtime();
// #ifdef PROFILE
// 	fprintf(stderr, "GPU potential (correction, float) - rank %d; Ni: %d; NTOT: %d; Time: %f sec\n",*rank,ni,n,t1 - t0);
// #endif
}

void gpupot_send(
        int *rank,
		int n,
		double m[],
		double x[][3]){
  	gpunb_devinit_(rank);
    tsend -= get_wtime();
	int ntg = NTHREAD * (n/NTHREAD + (n%NTHREAD ? 1 : 0));
	for(int i=0; i<n; i++){
      ptcl[i] = Particle(x[i], m[i]);
	}
    for(int i=n; i<ntg; i++){
      ptcl[i] = Particle(0);
    }
	ptcl.htod(ntg);
	tsend += get_wtime();
// #ifdef PROFILE
// 	fprintf(stderr, "GPU potential send - rank %d; NTOT: %d; Time: %f sec\n",*rank,n,t1 - t0);
// #endif
}    

void gpupot(
        int *rank,
        int ni,
		int n,
        int ishift,
        int list[],
        double dm[],
		double pot[]){
	gpunb_devinit_(rank);

    //    assert(ni<=NIMAX);
    //DEBUG=============================================//
    //    printf("ni %d n %d list[0] %d m %lf\n",ni,n,list[0],m[0]);
	tcalc -= get_wtime();
    icall++;
    ini +=ni;
	// cudaPointer <float2> phi;
    // cudaPointer <float2[NJBLOCK]> phipart;
	// cudaPointer <Particle> ptcl;
    // cudaPointer <Particle> ipbuf;
    //    cudaPointer <int> plist;
    int ng = NTHREAD * (ni/NTHREAD + (ni%NTHREAD ? 1 : 0));
    //	int ntg = NTHREAD * (n/NTHREAD + (n%NTHREAD ? 1 : 0));

    // ibuf.allocate(ng);
	// phi.allocate(ng);
    // phipart.allocate(ng);
	// ptcl.allocate(ntg);
    //    plist.allocate(ni);

    //  std::cout << n << " " << ng << " "<< ntg << std::endl;
    // for(int i=0; i<ni; i++) {
    //   plist[i] = list[i] - 1;
    // }
    for(int i=0; i<ni; i++) {
      ibuf[i] = ptcl[list[i]-1+ishift];
      ibuf[i].mass = dm[i];
    }
    for(int i=ni; i<ng; i++) {
      ibuf[i] = Particle(0);
    }
    
    //    printf("plist[0] %d\n",plist[0]);
    
    ibuf.htod(ng);
    //    plist.htod(ni);
    //    ptcl.dtoh(ntg);
	dim3 grid(ng/NTHREAD, NJBLOCK, 1);
	dim3 threads(NTHREAD, 1, 1);
    //	int sharedMemSize = NTHREAD * sizeof(Particle);
	pot_kernel <<<grid, threads>>> (ni, n, ibuf, ptcl, phipart);

    //    ptcl.dtoh(ni);
    
    const int ni8 = 1 + (ni-1) / NYREDUCE;
    dim3 rgrid (ni8, 1, 1);
    dim3 rthreads(NXREDUCE, NYREDUCE, 1);
    pot_reduce_kernel <<< rgrid, rthreads >>> (ni, phipart, phi);

    phi.dtoh(ni);

    for(int i=0; i<ni; i++){
      pot[i] = (double)phi[i].x + (double)phi[i].y;
	}

	// phi.free();
    // phipart.free();
	// ptcl.free();
    // ipbuf.free();
	tcalc += get_wtime();
// #ifdef PROFILE
// 	fprintf(stderr, "GPU potential (correction) - rank %d; Ni: %d; NTOT: %d; Time: %f sec\n",*rank,ni,n,t1 - t0);
// #endif
}

void gpupot_profile(int irank) {
#ifdef PROFILE
  if(icall) {
    // R.: rank;
    // Ncall: number of call gpupot and gpupot_float during two checking time(adjust time interval)
    // <Ni>: averaged i particles per call;
    // send: j particle sending time;
    // pot:  potential calculation time;
    fprintf(stderr, "[R.%d GPU Pot.C] Ncall %d  <Ni> %d   send(s) %f (ave) %f  pot(s) %f (ave) %f\n",irank,icall,ini/icall,tsend,tsend/icall,tcalc,tcalc/icall);
  }
  ini = icall = 0;
  tsend = tcalc = 0.;
#else
  return;
#endif
}

extern "C"{
  void gpupot_init_(
             int *irank,
             int *n){
    gpupot_init(irank,*n);
  }
  void gpupot_close_(){
    gpupot_close();
  }
  void gpupot_send_(int *irank, int *n, double m[], double x[][3]) {
    gpupot_send(irank,*n,m,x);
  }
  void gpupot_dm_(
            int *irank,
            int *ni,
			int *n,
            int *ishift,
            int list[],
            double dm[],
			double pot[]){
    gpupot(irank, *ni, *n, *ishift, list, dm, pot);
  }
  void gpupot_init_float_(
           int *irank ){
    gpupot_init_float(irank);
  }
  void gpupot_close_float_(){
    gpupot_close_float();
  }
  void gpupot_float_(
            int *irank,
            int *ni,
			int *n,
            int *ishift,
            int list[],
            double dm[],
			double pot[]){
    gpupot_float(irank, *ni, *n, *ishift, list, dm, pot);
  }
  void gpupot_mdot_profile_(int *irank){
    gpupot_profile(*irank);
  }
}
