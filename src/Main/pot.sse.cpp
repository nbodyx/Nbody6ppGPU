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

//typedef float v4sf __attribute__ ((vector_size(16)));
static inline v4sf v4sf_rsqrt(v4sf x){
	v4sf y = __builtin_ia32_rsqrtps(x);
	return ((v4sf){-0.5f, -0.5f, -0.5f, -0.5f} * y) * 
			(x*y*y + (v4sf){-3.f, -3.f, -3.f, -3.f});
}

static inline void pot_reduce(v4sf potH, v4sf potL, double pot[]){
	pot[0] = (double)__builtin_ia32_vec_ext_v4sf(potH, 0)
		   + (double)__builtin_ia32_vec_ext_v4sf(potL, 0);
	pot[1] = (double)__builtin_ia32_vec_ext_v4sf(potH, 1)
		   + (double)__builtin_ia32_vec_ext_v4sf(potL, 1);
	pot[2] = (double)__builtin_ia32_vec_ext_v4sf(potH, 2)
		   + (double)__builtin_ia32_vec_ext_v4sf(potL, 2);
	pot[3] = (double)__builtin_ia32_vec_ext_v4sf(potH, 3)
		   + (double)__builtin_ia32_vec_ext_v4sf(potL, 3);
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

	Particle *ptcl = new Particle[n+4];
	for(int i=0; i<n; i++){
		ptcl[i] = Particle(x[i], m[i]);
	}

    const int ibegin = istart - 1;
    const int iend = ni + istart - 1;

#pragma omp parallel for
	for(int i=ibegin; i<iend; i+=4){
		v4sf potH = {0.f, 0.f, 0.f, 0.f};
		v4sf potL = {0.f, 0.f, 0.f, 0.f};
		Particle *p = ptcl + i;
		v4sf xiH = {p[0].pos[0].x, p[1].pos[0].x, p[2].pos[0].x, p[3].pos[0].x};
		v4sf yiH = {p[0].pos[1].x, p[1].pos[1].x, p[2].pos[1].x, p[3].pos[1].x};
		v4sf ziH = {p[0].pos[2].x, p[1].pos[2].x, p[2].pos[2].x, p[3].pos[2].x};
		v4sf xiL = {p[0].pos[0].y, p[1].pos[0].y, p[2].pos[0].y, p[3].pos[0].y};
		v4sf yiL = {p[0].pos[1].y, p[1].pos[1].y, p[2].pos[1].y, p[3].pos[1].y};
		v4sf ziL = {p[0].pos[2].y, p[1].pos[2].y, p[2].pos[2].y, p[3].pos[2].y};
		for(int j=0; j<n; j++){
			v4sf jp0 = ((v4sf *)&ptcl[j])[0];
			v4sf jp1 = ((v4sf *)&ptcl[j])[1];
			v4sf xjH = __builtin_ia32_shufps(jp0, jp0, 0x00);
			v4sf xjL = __builtin_ia32_shufps(jp0, jp0, 0x55);
			v4sf yjH = __builtin_ia32_shufps(jp0, jp0, 0xaa);
			v4sf yjL = __builtin_ia32_shufps(jp0, jp0, 0xff);
			v4sf zjH = __builtin_ia32_shufps(jp1, jp1, 0x00);
			v4sf zjL = __builtin_ia32_shufps(jp1, jp1, 0x55);
			v4sf mj  = __builtin_ia32_shufps(jp1, jp1, 0xaa);
			
			v4sf dx = (xjH - xiH) + (xjL - xiL);
			v4sf dy = (yjH - yiH) + (yjL - yiL);
			v4sf dz = (zjH - ziH) + (zjL - ziL);
			v4sf r2 = dx*dx + dy*dy + dz*dz;
			v4sf mask = (v4sf)__builtin_ia32_cmpltps((v4sf){0,0,0,0}, r2);
			v4sf rinv = v4sf_rsqrt(r2);
			rinv = __builtin_ia32_andps(rinv, mask);
			rinv *= mj;

			v4sf tmp = potH;
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
