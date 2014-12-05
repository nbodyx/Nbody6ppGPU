#include <cstdio>
#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
#include "simd_define.h"

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

//#define REP4(x) {x, x, x, x}
//#define REP8(x) {x, x, x, x, x, x, x, x}

//typedef float  v4sf __attribute__ ((vector_size(16)));
//typedef float  v8sf __attribute__ ((vector_size(32)));
//typedef double v4df __attribute__ ((vector_size(32)));

static inline v8sf v8sf_rsqrt(const v8sf x){
	v8sf y = __builtin_ia32_rsqrtps256(x);
	return ((v8sf)REP8(-0.5f) * y) * (x*y*y + (v8sf)REP8(-3.0f));
}

static inline void pot_reduce(const v8sf potH, const v8sf potL, double pot[]){
	const v4df dpot0 = 
		  __builtin_ia32_cvtps2pd256(__builtin_ia32_vextractf128_ps256(potH, 0))
		+ __builtin_ia32_cvtps2pd256(__builtin_ia32_vextractf128_ps256(potL, 0));
	const v4df dpot4 = 
		  __builtin_ia32_cvtps2pd256(__builtin_ia32_vextractf128_ps256(potH, 1))
		+ __builtin_ia32_cvtps2pd256(__builtin_ia32_vextractf128_ps256(potL, 1));
	// *(v4df *)(pot + 0) = dpot0;
	// *(v4df *)(pot + 4) = dpot4;
	__builtin_ia32_storeupd256(pot + 0, dpot0);
	__builtin_ia32_storeupd256(pot + 4, dpot4);
}

struct float2{
	float x, y;
};
static inline float2 float2_split(double x){
	float2 ret;
	x *= (1<<16);
	double xi = (int)x;
	double xf = x - xi;
	ret.x = xi * (1./(1<<16));
	ret.y = xf * (1./(1<<16));
	return ret;
}

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
	Particle(){
		pos[0].x = pos[0].y = pos[1].x = pos[1].y = pos[2].x = pos[2].y = mass = pad = 0.f;
	}
};

void gpupot(
        int rank,
        int istart,
        int ni,
		int n,
		double m[],
		double x[][3],
		double pot[]){
	double t0 = get_wtime();

	Particle *ptcl = new Particle[n+8];
	for(int i=0; i<n; i++){
		ptcl[i] = Particle(x[i], m[i]);
	}

    const int ibegin = istart - 1;
    const int iend = ni + istart - 1;
#pragma omp parallel for
	for(int i=ibegin; i<iend; i+=8){
		v8sf potH = REP8(0.0);
		v8sf potL = REP8(0.0);
		const Particle *p = ptcl + i;
		const v8sf xiH = 
			{p[0].pos[0].x, p[1].pos[0].x, p[2].pos[0].x, p[3].pos[0].x,
			 p[4].pos[0].x, p[5].pos[0].x, p[6].pos[0].x, p[7].pos[0].x};
		const v8sf yiH = 
			{p[0].pos[1].x, p[1].pos[1].x, p[2].pos[1].x, p[3].pos[1].x,
			 p[4].pos[1].x, p[5].pos[1].x, p[6].pos[1].x, p[7].pos[1].x};
		const v8sf ziH = 
			{p[0].pos[2].x, p[1].pos[2].x, p[2].pos[2].x, p[3].pos[2].x,
			 p[4].pos[2].x, p[5].pos[2].x, p[6].pos[2].x, p[7].pos[2].x};
		const v8sf xiL = 
			{p[0].pos[0].y, p[1].pos[0].y, p[2].pos[0].y, p[3].pos[0].y,
			 p[4].pos[0].y, p[5].pos[0].y, p[6].pos[0].y, p[7].pos[0].y};
		const v8sf yiL = 
			{p[0].pos[1].y, p[1].pos[1].y, p[2].pos[1].y, p[3].pos[1].y,
			 p[4].pos[1].y, p[5].pos[1].y, p[6].pos[1].y, p[7].pos[1].y};
		const v8sf ziL = 
			{p[0].pos[2].y, p[1].pos[2].y, p[2].pos[2].y, p[3].pos[2].y,
			 p[4].pos[2].y, p[5].pos[2].y, p[6].pos[2].y, p[7].pos[2].y};
		for(int j=0; j<n; j++){
			const v8sf jp0 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&ptcl[j] + 0);
			const v8sf jp1 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&ptcl[j] + 1);
			const v8sf xjH = __builtin_ia32_shufps256(jp0, jp0, 0x00);
			const v8sf xjL = __builtin_ia32_shufps256(jp0, jp0, 0x55);
			const v8sf yjH = __builtin_ia32_shufps256(jp0, jp0, 0xaa);
			const v8sf yjL = __builtin_ia32_shufps256(jp0, jp0, 0xff);
			const v8sf zjH = __builtin_ia32_shufps256(jp1, jp1, 0x00);
			const v8sf zjL = __builtin_ia32_shufps256(jp1, jp1, 0x55);
			const v8sf mj  = __builtin_ia32_shufps256(jp1, jp1, 0xaa);
			
			const v8sf  dx = (xjH - xiH) + (xjL - xiL);
			const v8sf  dy = (yjH - yiH) + (yjL - yiL);
			const v8sf  dz = (zjH - ziH) + (zjL - ziL);
			const v8sf  r2 = dx*dx + dy*dy + dz*dz;
			const v8sf mask = __builtin_ia32_cmpps256((v8sf)REP8(0.0), r2, 17);
			                                    // 17 : less-than ordered quiet
			v8sf rinv = v8sf_rsqrt(r2);
			rinv = __builtin_ia32_andps256(rinv, mask);
			rinv *= mj;

			const v8sf tmp = potH;
			potH += rinv;
			potL -= (potH - tmp) - rinv;
		}
		pot_reduce(potH, potL, pot+i);
	}

	delete [] ptcl;

	double t1 = get_wtime();
#ifdef PROFILE
	fprintf(stderr, "[R.%d AVX Pot.A] Ni %d  NTOT %d  pot(s) %f\n", rank,ni,n,t1 - t0);
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
      gpupot(*irank,*istart,*ni,*n, m, x, pot);
	}
}
