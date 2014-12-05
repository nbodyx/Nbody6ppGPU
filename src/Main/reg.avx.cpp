#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "simd_define.h"

#define TMAX 32 // maximum number of threads
#if 1
#include <omp.h>
#else
static inline int omp_get_num_threads(){return 1;}
static inline int omp_get_thread_num() {return 0;}
#endif

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

static double time_send, time_grav;
static long long numInter;
static int icall,ini,isend;

template <class T>
struct myvector{
	size_t num;
	T *val;
	myvector() : num(0), val(NULL) {}
	~myvector(){
		delete [] val;
		val = NULL;
	}
	void clear(){
		num = 0;
	}
	void reserve(size_t count){
		val = new T[count];
	}
	void free(){
		delete [] val;
		val = NULL;
	}
	void push_back(const T &t){
		val[num++] = t;
	}
	size_t size() const {
		return num;
	}
	T &operator[](int i){
		return val[i];
	}
	const T &operator[](int i) const{
		return val[i];
	}
};

// typedef float  v4sf __attribute__ ((vector_size(16)));
// typedef float  v8sf __attribute__ ((vector_size(32)));
// typedef double v2df __attribute__ ((vector_size(16)));
// typedef double v4df __attribute__ ((vector_size(32)));

// #define REP4(x) {x, x, x, x}
// #define REP8(x) {x, x, x, x, x, x, x, x}

static inline v8sf v8sf_rsqrt(const v8sf x){
	v8sf y = __builtin_ia32_rsqrtps256(x);
	return ((v8sf)REP8(-0.5f) * y) * (x*y*y + (v8sf)REP8(-3.0f));
}

static v4sf *jparr1; // {x, y, z, m}
static v4sf *jparr2; // {vx, vy, vz, pad}

static myvector<int> nblist[TMAX][4];
static int nbody, nbodymax;

static void *amalloc64(size_t n){
	void *ptr;
	(void)posix_memalign(&ptr, 64, n);
	assert(ptr);
	return ptr;
}

void GPUNB_open(int nbmax, int irank){
	time_send = time_grav = 0.0;
	numInter = 0;
    icall = isend = ini = 0;
	nbodymax = nbmax;
	jparr1 = (v4sf *)amalloc64(sizeof(v4sf) * (nbmax+3));
	jparr2 = (v4sf *)amalloc64(sizeof(v4sf) * (nbmax+3));
	int numCPU = -1;
#pragma omp parallel
	{
		int nth = omp_get_num_threads();
		assert(nth <= TMAX);
		int tid = omp_get_thread_num();
		for(int v=0; v<4; v++){
			nblist[tid][v].reserve(nbmax);
		}
		if(tid == 0) numCPU = omp_get_num_threads();
	}
#ifdef PROFILE
    fprintf(stderr, "# Open AVX regular force - rank: %d; threads: %d\n", irank, numCPU);
#endif
}

void GPUNB_close(){
	free(jparr1); jparr1 = NULL;
	free(jparr2); jparr2 = NULL;
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		for(int v=0; v<4; v++){
			nblist[tid][v].free();
		}
	}
	nbodymax = 0;

// #ifdef PROFILE
// 	fprintf(stderr, "***********************\n");
// 	fprintf(stderr, "Closed NBODY6/AVX library\n");
// 	fprintf(stderr, "time send : %f sec\n", time_send);
// 	fprintf(stderr, "time grav : %f sec\n", time_grav);
// 	fprintf(stderr, "%f  Gflops (gravity part only)\n", 60.e-9 * numInter / time_grav);
// 	fprintf(stderr, "***********************\n");
// #endif
}

void GPUNB_send(
		int nj,
		double mj[],
		double xj[][3],
		double vj[][3])
{
	time_send -= get_wtime();
	nbody = nj;
    ::isend++;
	assert(nbody <= nbodymax);
#pragma omp parallel for
	for(int j=0; j<nj; j++){
		const v4df jp1 = {xj[j][0], xj[j][1], xj[j][2], mj[j]};
		const v4df jp2 = {vj[j][0], vj[j][1], vj[j][2], 0.0  };
		jparr1[j] = __builtin_ia32_cvtpd2ps256(jp1);
		jparr2[j] = __builtin_ia32_cvtpd2ps256(jp2);
	}
	if(nj%2){ // padding
		jparr1[nj] = (v4sf){255.0f, 255.0f, 255.0f, 0.0f};
		jparr2[nj] = (v4sf){0.0f, 0.0f, 0.0f, 0.0f};
		nbody++;
	}
	time_send += get_wtime();
}

static inline v8sf gen_i_particle(double x, double y, double z, double w){
	const v4df vd = {x, y, z, w};
	const v4sf vs = __builtin_ia32_cvtpd2ps256(vd);
	v8sf ret = REP8(0.0f);
	ret = __builtin_ia32_vinsertf128_ps256(ret, vs, 0);
	ret = __builtin_ia32_vinsertf128_ps256(ret, vs, 1);
	return ret;
}

static inline void reduce_force(const v8sf f8, double &x, double &y, double &z, double &w){
	const v4sf fh = __builtin_ia32_vextractf128_ps256(f8, 0);
	const v4sf fl = __builtin_ia32_vextractf128_ps256(f8, 1);
	const v4df fsum = __builtin_ia32_cvtps2pd256(fh)
	                + __builtin_ia32_cvtps2pd256(fl);
	const v2df xy = __builtin_ia32_vextractf128_pd256(fsum, 0);
	const v2df zw = __builtin_ia32_vextractf128_pd256(fsum, 1);
	x = __builtin_ia32_vec_ext_v2df(xy, 0);
	y = __builtin_ia32_vec_ext_v2df(xy, 1);
	z = __builtin_ia32_vec_ext_v2df(zw, 0);
	w = __builtin_ia32_vec_ext_v2df(zw, 1);

}

void GPUNB_regf(
		int ni,
		double h2d[],
		double dtr[],
		double xid[][3],
		double vid[][3],
		double acc[][3],
		double jrk[][3],
		double pot[],
		int lmax,
		int nbmax,
		int *listbase,
        int m_flag)
{
	time_grav -= get_wtime();
	numInter += ni * nbody;
    ::icall++;
    ::ini +=ni;
#pragma omp parallel for
	for(int i=0; i<ni; i+=4){
		int tid = omp_get_thread_num();
		nblist[tid][0].clear();
		nblist[tid][1].clear();
		nblist[tid][2].clear();
		nblist[tid][3].clear();
		int nii = (ni-i < 4) ? ni-i : 4;

		const v8sf xi  = gen_i_particle(xid[i+0][0], xid[i+1][0], xid[i+2][0], xid[i+3][0]);
		const v8sf yi  = gen_i_particle(xid[i+0][1], xid[i+1][1], xid[i+2][1], xid[i+3][1]);
		const v8sf zi  = gen_i_particle(xid[i+0][2], xid[i+1][2], xid[i+2][2], xid[i+3][2]);
		const v8sf vxi = gen_i_particle(vid[i+0][0], vid[i+1][0], vid[i+2][0], vid[i+3][0]);
		const v8sf vyi = gen_i_particle(vid[i+0][1], vid[i+1][1], vid[i+2][1], vid[i+3][1]);
		const v8sf vzi = gen_i_particle(vid[i+0][2], vid[i+1][2], vid[i+2][2], vid[i+3][2]);
		static const v8sf h2mask[5] = {
			{0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},
			{1.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0},
			{1.0, 1.0, 0.0, 0.0,  1.0, 1.0, 0.0, 0.0},
			{1.0, 1.0, 1.0, 0.0,  1.0, 1.0, 1.0, 0.0},
			{1.0, 1.0, 1.0, 1.0,  1.0, 1.0, 1.0, 1.0},
		};
		const v8sf h2i  = gen_i_particle(h2d[i+0], h2d[i+1], h2d[i+2], h2d[i+3]) * h2mask[nii]; 
		const v8sf dtri = gen_i_particle(dtr[i+0], dtr[i+1], dtr[i+2], dtr[i+3]); 

		v8sf Ax = REP8(0.0f);
		v8sf Ay = REP8(0.0f);
		v8sf Az = REP8(0.0f);
		v8sf Jx = REP8(0.0f);
		v8sf Jy = REP8(0.0f);
		v8sf Jz = REP8(0.0f);
		v8sf poti = REP8(0.0f);

		const v4sf *jpp1 = jparr1;
		const v4sf *jpp2 = jparr2;
		const int nj = ::nbody;
		for(int j=0; j<nj; j+=2){
			const v8sf jp1 = *(v8sf *)(jpp1 + j);
			const v8sf jp2 = *(v8sf *)(jpp2 + j);
			const v8sf  xj  = __builtin_ia32_shufps256(jp1, jp1, 0x00);
			const v8sf  yj  = __builtin_ia32_shufps256(jp1, jp1, 0x55);
			const v8sf  zj  = __builtin_ia32_shufps256(jp1, jp1, 0xaa);
			const v8sf  mj  = __builtin_ia32_shufps256(jp1, jp1, 0xff);
			const v8sf  vxj = __builtin_ia32_shufps256(jp2, jp2, 0x00);
			const v8sf  vyj = __builtin_ia32_shufps256(jp2, jp2, 0x55);
			const v8sf  vzj = __builtin_ia32_shufps256(jp2, jp2, 0xaa);

			const v8sf dx = xj - xi;
			const v8sf dy = yj - yi;
			const v8sf dz = zj - zi;
			const v8sf dvx = vxj - vxi;
			const v8sf dvy = vyj - vyi;
			const v8sf dvz = vzj - vzi;

			const v8sf dxp = dx + dtri * dvx;
			const v8sf dyp = dy + dtri * dvy;
			const v8sf dzp = dz + dtri * dvz;

			const v8sf r2  = dx*dx + dy*dy + dz*dz;
			      v8sf rv  = dx*dvx + dy*dvy + dz*dvz;
			const v8sf r2p = dxp*dxp + dyp*dyp + dzp*dzp;

			const v8sf r2min = __builtin_ia32_minps256(r2, r2p);
            v8sf mask;
            if(m_flag) {
              const v8sf mh2i = mj * h2i;
              mask = __builtin_ia32_cmpps256(r2min, mh2i, 17); // 17 : less-than ordered quiet
            }
            else{
              mask = __builtin_ia32_cmpps256(r2min, h2i, 17); // 17 : less-than ordered quiet
            }

			const int bits = __builtin_ia32_movmskps256(mask);
			if(bits){
			//if(bits & 0x0f){
				if (bits&1) nblist[tid][0].push_back(j);
				if (bits&2) nblist[tid][1].push_back(j);
				if (bits&4) nblist[tid][2].push_back(j);
				if (bits&8) nblist[tid][3].push_back(j);
			//}
			//if(bits & 0xf0){
				if (bits&0x10) nblist[tid][0].push_back(j+1);
				if (bits&0x20) nblist[tid][1].push_back(j+1);
				if (bits&0x40) nblist[tid][2].push_back(j+1);
				if (bits&0x80) nblist[tid][3].push_back(j+1);
			//}
			}

			v8sf rinv1 = v8sf_rsqrt(r2);
			rinv1 = __builtin_ia32_andnps256(mask, rinv1);
			const v8sf rinv2 = rinv1 * rinv1;
			rinv1 *= mj;
			poti += rinv1;

			const v8sf rinv3 = rinv1 * rinv2;
			rv *= (v8sf)REP8(-3.0f) * rinv2;

			Ax += rinv3 * dx;
			Ay += rinv3 * dy;
			Az += rinv3 * dz;
			Jx += rinv3 * (dvx + rv * dx);
			Jy += rinv3 * (dvy + rv * dy);
			Jz += rinv3 * (dvz + rv * dz);
		}

		reduce_force(Ax, acc[i+0][0], acc[i+1][0], acc[i+2][0], acc[i+3][0]);
		reduce_force(Ay, acc[i+0][1], acc[i+1][1], acc[i+2][1], acc[i+3][1]);
		reduce_force(Az, acc[i+0][2], acc[i+1][2], acc[i+2][2], acc[i+3][2]);
		reduce_force(Jx, jrk[i+0][0], jrk[i+1][0], jrk[i+2][0], jrk[i+3][0]);
		reduce_force(Jy, jrk[i+0][1], jrk[i+1][1], jrk[i+2][1], jrk[i+3][1]);
		reduce_force(Jz, jrk[i+0][2], jrk[i+1][2], jrk[i+2][2], jrk[i+3][2]);
		reduce_force(poti, pot[i+0], pot[i+1], pot[i+2], pot[i+3]);

		for(int ii=0; ii<nii; ii++){
			const int nnb = nblist[tid][ii].size();
			int *nnbp = listbase + lmax * (i+ii);
			int *nblistp = nnbp + 1;
			if(nnb > nbmax){
				*nnbp = -nnb;
			}else{
				*nnbp = nnb;
				for(int k=0; k<nnb; k++){
					nblistp[k] = nblist[tid][ii][k];
				}
			}
		}
	}
	time_grav += get_wtime();
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
    // Perf: performance for gpu regular calculation
    fprintf(stderr,"[R.%d AVX Reg.F ] Nsend %d  Ngrav %d  <Ni> %d   send(s) %f grav(s) %f  Perf.(Gflops) %f\n",irank,isend,icall,ini/isend,time_send,time_grav,60.e-9*numInter/time_grav);
  }
  time_send = time_grav = 0.0;
  numInter = 0;
  icall = ini = isend= 0;
#else
  return;
#endif
}

extern "C" {
    void gpunb_open_( int *nbmax, int *irank){
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
            int *m_flag){ // list[][lmax]
         GPUNB_regf(*ni, h2, dtr, xi, vi, acc, jrk, pot, *lmax, *nbmax, list,*m_flag);
	}
    void gpunb_profile_(int *irank){
      GPUNB_profile(*irank);
    }
}
