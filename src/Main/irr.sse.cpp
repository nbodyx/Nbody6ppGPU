#define __USE_GNU
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <unistd.h>
#include <omp.h>
#include "vector3.h"
#include "simd_define.h"

#define NNBMAX 600

#define _out_
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

//typedef float  v4sf __attribute__((vector_size(16)));
//typedef double v2df __attribute__((vector_size(16)));

//#define REP4(x) {x,x,x,x}

static inline v4sf rsqrt_NR(const v4sf x){
	const v4sf y = __builtin_ia32_rsqrtps(x);
	const v4sf c1 = REP4(-0.5f);
	const v4sf c2 = REP4(-3.0f);
	return (c1 * y) * (x*y*y + c2);
}

struct Particle{
	v4sf posH; // elemen w is for mass
	v4sf posL;
	v4sf vel;
	v4sf acc2;
	v4sf jrk6;
	v2df time; // 6 xmm words;

	void make_pos(const double *x, const double *m){
		v2df xy = {x[0], x[1]};
		v2df zw = {x[2], m[0]};
		v4sf xyH = __builtin_ia32_cvtpd2ps(xy);
		v4sf xyL = __builtin_ia32_cvtpd2ps(xy - __builtin_ia32_cvtps2pd(xyH));
		v4sf zwH = __builtin_ia32_cvtpd2ps(zw);
		v4sf zwL = __builtin_ia32_cvtpd2ps(zw - __builtin_ia32_cvtps2pd(zwH));

		posH = __builtin_ia32_movlhps(xyH, zwH);
		posL = __builtin_ia32_movlhps(xyL, zwL);
	}

    static v4sf make_v4sf(const double *x){
      v2df xy = {x[0], x[1]};
      v2df zw = {x[2],     };
      return __builtin_ia32_movlhps(
				__builtin_ia32_cvtpd2ps(xy),
				__builtin_ia32_cvtpd2ps(zw));
	}

    Particle(
		const double _pos [3],
		const double _vel [3],
		const double _acc2[3],
		const double _jrk6[3],
		const double _mass,
		const double _time)
	{
		make_pos(_pos, &_mass);
		vel  = make_v4sf(_vel);
		acc2 = make_v4sf(_acc2);
		jrk6 = make_v4sf(_jrk6);
		time = (v2df){_time, _time};
	}

	Particle(int) // constructor for a dummy particle
	{
		posH = (v4sf){255.0f, 255.0f, 255.0f, 0.0f};
		posL = vel = acc2 = jrk6 = (v4sf)REP4(0.0f);
		time = (v2df){0.0, 0.0};
	}

	void prefetch(const int rw=0, const int locality=3) const {
		const char *addr = (char *)&posH;
		__builtin_prefetch( 0+addr, rw , locality);
		__builtin_prefetch(64+addr, rw , locality);
	}
	bool is_aligned() const {
		unsigned long addr = (unsigned long)(&posH);
		return (0 == addr%32);
	}
	void store(Particle *pdst) const {
		const v4sf xm0 = posH;
		const v4sf xm1 = posL;
		const v4sf xm2 = vel;
		const v4sf xm3 = acc2;
		const v4sf xm4 = jrk6;
		const v2df xm5 = time;
		v4sf *dst = (v4sf *)pdst;
		dst[0] = xm0;
		dst[1] = xm1;
		dst[2] = xm2;
		dst[3] = xm3;
		dst[4] = xm4;
		dst[5] = (v4sf)xm5;
	}
	void sstore(Particle *pdst) const {
		const v4sf xm0 = posH;
		const v4sf xm1 = posL;
		const v4sf xm2 = vel;
		const v4sf xm3 = acc2;
		const v4sf xm4 = jrk6;
		const v2df xm5 = time;
		v4sf *dst = (v4sf *)pdst;
		__builtin_ia32_movntps((float *)(dst + 0), xm0);
		__builtin_ia32_movntps((float *)(dst + 1), xm1);
		__builtin_ia32_movntps((float *)(dst + 2), xm2);
		__builtin_ia32_movntps((float *)(dst + 3), xm3);
		__builtin_ia32_movntps((float *)(dst + 4), xm4);
		__builtin_ia32_movntps((float *)(dst + 5), (v4sf)xm5);
	}
};

struct Predictor{
	v4sf posH; // elemen w is for mass
	v4sf posL;
	v4sf vel;  // 3 xmm words

    Predictor(const Particle &p, const v2df ti)
	{
        const v4sf dt = __builtin_ia32_cvtpd2ps(ti - v2df(p.time));
		const v4sf s0 = __builtin_ia32_shufps(dt, dt, 0x00);
		const v4sf s1 = s0 + s0;
		const v4sf s2 = s0 * (v4sf)REP4(1.5f);

		this->posH = p.posH;
		this->posL = v4sf(p.posL) + s0*(v4sf(p.vel) + s0*(v4sf(p.acc2) + s0*(v4sf(p.jrk6))));
		this->vel = v4sf(p.vel) + s1*(v4sf(p.acc2) + s2*(v4sf(p.jrk6)));
	}
};

struct Pred4{
	v4sf xH, yH, zH, mass;
	v4sf xL, yL, zL;
	v4sf vx, vy, vz;

	static v4sf bcast(const float f){
		v4sf v = {f, f, f, f};
		return v;
	}
	static v4sf bcast0(const v4sf v){
		return __builtin_ia32_shufps(v, v, 0);
	}
	static v4sf bcast1(const v4sf v){
		return __builtin_ia32_shufps(v, v, 0x55);
	}
	static v4sf bcast2(const v4sf v){
		return __builtin_ia32_shufps(v, v, 0xaa);
	}
	static v4sf bcast3(const v4sf v){
		return __builtin_ia32_shufps(v, v, 0xff);
	}
	static void transpose(v4sf &v0, v4sf &v1, v4sf &v2, v4sf &v3){
		const v4sf t0 = __builtin_ia32_unpcklps(v0, v2);
		const v4sf t1 = __builtin_ia32_unpckhps(v0, v2);
		const v4sf t2 = __builtin_ia32_unpcklps(v1, v3);
		const v4sf t3 = __builtin_ia32_unpckhps(v1, v3);

		v0 = __builtin_ia32_unpcklps(t0, t2);
		v1 = __builtin_ia32_unpckhps(t0, t2);
		v2 = __builtin_ia32_unpcklps(t1, t3);
		v3 = __builtin_ia32_unpckhps(t1, t3);
	}

    Pred4(){}

	// for the i-particle
	Pred4(const Predictor &pi){
		xH   = bcast0(pi.posH);
		yH   = bcast1(pi.posH);
		zH   = bcast2(pi.posH);
		mass = bcast3(pi.posH);
		xL   = bcast0(pi.posL);
		yL   = bcast1(pi.posL);
		zL   = bcast2(pi.posL);
		vx   = bcast0(pi.vel );
		vy   = bcast1(pi.vel );
		vz   = bcast2(pi.vel );
	}
	// for the j-particle
	Pred4(const Predictor &p0, const Predictor &p1, const Predictor &p2, const Predictor &p3){
		xH   = p0.posH;
		yH   = p1.posH;
		zH   = p2.posH;
		mass = p3.posH;
		transpose(xH, yH, zH, mass);
		v4sf dum0;
		xL   = p0.posL;
		yL   = p1.posL;
		zL   = p2.posL;
		dum0 = p3.posL;
		transpose(xL, yL, zL, dum0);
		v4sf dum1;
		vx   = p0.vel;
		vy   = p1.vel;
		vz   = p2.vel;
		dum1 = p3.vel;
		transpose(vx, vy, vz, dum1);
	}
};

// single precision version
struct Force{
	v4sf ax, ay, az;
	v4sf jx, jy, jz;
    v4sf vnnb;

	void clear(){
		const v4sf zero = REP4(0.0f);
		ax = ay = az = zero;
		jx = jy = jz = zero;

        v2df tmp = {HUGE,HUGE};
        vnnb = (v4sf)tmp;
	}
 
	static double reduce(const v4sf v){
        const v2df lo = __builtin_ia32_cvtps2pd(v);
        const v2df hi =  __builtin_ia32_cvtps2pd(__builtin_ia32_movhlps(v, v));
		const v2df tmpd = lo + hi;
		const v2df sum  = __builtin_ia32_haddpd(tmpd, tmpd);
		return __builtin_ia32_vec_ext_v2df(sum , 0);
	}

    int reduce_nnb() const{
        // min( [i1|r1|i2|r2] , [i2|r2|i2|r2] )
        v2df min = __builtin_ia32_minpd((v2df)vnnb,
                         __builtin_ia32_unpckhpd((v2df)vnnb, (v2df)vnnb));
        // [i1|r1] <-> [i|]
        union{
            v2df v;
            int  i;
        } mem;
        mem.v = min;
        return mem.i;
    }

  void write(double *acc, double *jrk, int &nnbid) const{

		const double a0 = reduce(ax);
		const double a1 = reduce(ay);
		const double a2 = reduce(az);
		const double j0 = reduce(jx);
		const double j1 = reduce(jy);
		const double j2 = reduce(jz);
        const int    ii = reduce_nnb() + 1;
		acc[0] = a0;
		acc[1] = a1;
		acc[2] = a2;
		jrk[0] = j0;
		jrk[1] = j1;
		jrk[2] = j2;
        nnbid  = ii;
	}
  void calc_and_accum(const Pred4 &pi, const Pred4 &pj, const v4sf idx){
		const v4sf dx = (pj.xH - pi.xH) + (pj.xL - pi.xL);
		const v4sf dy = (pj.yH - pi.yH) + (pj.yL - pi.yL);
		const v4sf dz = (pj.zH - pi.zH) + (pj.zL - pi.zL);

		const v4sf dvx = pj.vx - pi.vx;
		const v4sf dvy = pj.vy - pi.vy;
		const v4sf dvz = pj.vz - pi.vz;

		const v4sf r2 = dx*dx  + dy*dy  + dz*dz;
		const v4sf rv = dx*dvx + dy*dvy + dz*dvz;
		const v4sf rinv   = rsqrt_NR(r2);
		const v4sf rinv2  = rinv * rinv;
		const v4sf c1     = REP4(-3.0f);
		const v4sf alpha  = c1 * rinv2 * rv;
		const v4sf mrinv3 = pj.mass * rinv * rinv2;

        // idx     = [i1|i2|i3|i4]; r2      = [r1|r2|r3|r4]
        // r2_idx0 = [i1|r1|i2|r2]; r2_idx1 = [i3|r3|i4|r4]
        const v4sf r2_idx0 = __builtin_ia32_unpcklps(idx, r2);
        const v4sf r2_idx1 = __builtin_ia32_unpckhps(idx, r2);
        // Find minimum distance
        vnnb = (v4sf)__builtin_ia32_minpd((v2df)vnnb, (v2df)r2_idx0);
        vnnb = (v4sf)__builtin_ia32_minpd((v2df)vnnb, (v2df)r2_idx1);

		ax += mrinv3 * dx;
		ay += mrinv3 * dy;
		az += mrinv3 * dz;

		jx += mrinv3 * (dvx + alpha * dx); 
		jy += mrinv3 * (dvy + alpha * dy); 
		jz += mrinv3 * (dvz + alpha * dz); 
	}
};

struct NBlist{
	enum{ NB_MAX = NNBMAX };
	int pad[7];
	int nnb;
	int nb[NB_MAX];

	NBlist(){
		nnb = 0;
		for(int i=0; i<NB_MAX; i++){
			nb[i] = 0;
		}
	}

	void print(const int i, FILE *fp = stdout) const{
		fprintf(fp, "%6d%6d :", i, nnb);
		for(int k=0; k<nnb; k++){
			fprintf(fp, " %d", nb[k]);
		}
		fprintf(fp, "\n");
		fflush(fp);
	}
};

// The irregular force library
static bool       is_open = false;
static Particle  *ptcl = NULL;
static NBlist    *list = NULL;
static int        nmax;
static int        num_threads;
static double     time_grav;
static unsigned long long num_inter, num_fcall, num_steps;
static v2df       vec_tnow;

static void irr_simd_open(
		const int nmax,
		const int lmax,
        const int rank)
{
	if(is_open){
		fprintf(stderr, "irr_simd: it is already open\n");
		return;
	}

	assert(lmax <= 1 + NBlist::NB_MAX);

	fprintf(stderr, "# Opening IRR_simd lib. SSE ver. - rank: %d; nmax: %d, lmax: %d\n",rank,nmax,lmax); 
	assert(0 == sizeof(NBlist)%16);

	void *ptr;
	assert( 0 == posix_memalign(&ptr, 64, (1+nmax) * sizeof(Particle)) );
	memset(ptr, 0xff, (1+nmax) * sizeof(Particle));
	ptcl = (Particle *)ptr;

	// store dummy predictor to the last of the array
	ptcl[nmax] = Particle(0);

	// list.resize(nmax);
	assert( 0 == posix_memalign(&ptr, 64, nmax * sizeof(NBlist)) );
	memset(ptr, 0xff, nmax * sizeof(NBlist));
	list = (NBlist *)ptr;

	::nmax = nmax;
#pragma omp parallel
	{
		num_threads = omp_get_num_threads();
	}
	time_grav = 0.0;
	num_inter = num_fcall = num_steps = 0;

	is_open = true;
}

static void irr_simd_profile(
        const int rank)
{
  if(is_open){
	const double Gflops = 60.0 * double(num_inter) * 1.e-9 / time_grav;
	const double usec_fcall    = 1.e6 * (time_grav / num_fcall);
	const double nnb_avr = double(num_inter) / double(num_steps);
    const int ni_avr = num_inter / num_fcall;

	fprintf(stderr, "[R.%d SSE Irr.F ] Ncall: %llu <NI>: %d <NB>: %f grav: %f s, %f Gflops, %f usec\n",
            rank, num_fcall, ni_avr, nnb_avr, time_grav, Gflops, usec_fcall);

    time_grav = 0.0;
    num_inter = num_fcall = num_steps = 0;
  }
}

static void irr_simd_close(
        const int rank)
{
	if(!is_open){
		fprintf(stderr, "irr_simd: it is already close\n");
		return;
	}

	free(ptcl); ptcl = NULL;
	free(list); list = NULL;

	fprintf(stderr, "Closing IRR_simd lib. CPU ver. - rank: %d\n",rank); 

	is_open = false;
}

static void irr_simd_set_jp(
		const int addr,
		const double pos [3],
		const double vel [3],
		const double acc2[3],
		const double jrk6[3],
		const double mass,
		const double time)
{
	Particle(pos, vel, acc2, jrk6, mass, time).store(&ptcl[addr]);
}

static void irr_simd_set_list(
		const int addr,
		const int nnb,
		const int nblist[])
{
	assert(nnb <= NBlist::NB_MAX);
	list[addr].nnb = nnb;

	const int *src = nblist;
	      int *dst = list[addr].nb;

	// for(int k=0; k<nnb; k+=4){
	// 	const int i0 = src[k+0] - 1;
	// 	const int i1 = src[k+1] - 1;
	// 	const int i2 = src[k+2] - 1;
	// 	const int i3 = src[k+3] - 1;
	// 	dst[k+0] = i0;
	// 	dst[k+1] = i1;
	// 	dst[k+2] = i2;
	// 	dst[k+3] = i3;
	// }
	
    for(int k=0; k<nnb; k+=4){

    // assert((unsigned long)dst %16 == 0);
	 	typedef int       v4si __attribute__((vector_size(16)));
	 	typedef long long v2di __attribute__ ((__vector_size__ (16)));
	 	const v4si one = REP4(1);
	 	const v4si idx0 = (v4si)__builtin_ia32_loaddqu((const char *)(src+k+0));
	 	__builtin_ia32_movntdq((v2di *)(dst+k+0), (v2di)(idx0-one));
	}
    
	// fill dummy
	const int nmax = ::nmax;
	const int kmax = 4 + 4 * (1 + (nnb-1)/4);
	for(int k=nnb; k<kmax; k++){
		dst[k] = nmax;
	}

}

static void irr_simd_firr(
		const int addr,
		_out_ double accout[3],
		_out_ double jrkout[3],
        _out_ int    &nnbid)
{
	const Particle *ptcl = ::ptcl;
	const v2df      tnow = ::vec_tnow;
	const int  nnb  = list[addr].nnb;

	Force force;
	force.clear();

	const Pred4 pri(Predictor(ptcl[addr], tnow));

	for(int k=0; k<nnb; k++){
		const int jaddr = list[addr].nb[k];
		ptcl[jaddr].prefetch();
	}

	for(int k=0; k<nnb; k+=4){
		const int *jptr = &(list[addr].nb[k]);
		const Predictor p0(ptcl[jptr[0]], tnow);
		const Predictor p1(ptcl[jptr[1]], tnow);
		const Predictor p2(ptcl[jptr[2]], tnow);
		const Predictor p3(ptcl[jptr[3]], tnow);

		const Pred4 prj(p0, p1, p2, p3);

        assert(0 == (size_t)jptr % 16);
        v4sf idx = *(v4sf *)jptr;
        force.calc_and_accum(pri, prj, idx);
	}

	force.write(accout, jrkout, nnbid);
}

static void irr_simd_firr_vec(
        const double ti,
		const int ni,
		const int addr[],
		_out_ double accout[][3],
		_out_ double jrkout[][3],
        _out_ int nnbid[])
{
  //  printf("NI %d TIME %f FIRST %d",ni,ti,addr[0]);
	const double t0 = get_wtime();
    ::vec_tnow = (v2df){ti, ti};
	int ninter = 0;

#pragma omp parallel for reduction(+: ninter) schedule(guided)
	for(int i=0; i<ni; i++){
      irr_simd_firr(addr[i]-1, accout[i], jrkout[i], nnbid[i]);
		ninter += list[addr[i]-1].nnb;
	}
	::num_inter += ninter;
	const double t1 = get_wtime();
	::time_grav += t1-t0;
	::num_fcall++;
	::num_steps += ni;
}

// FORTRAN interface
extern "C"{
    void irr_simd_open_(int *nmax, int *lmax, int *rank){
        irr_simd_open(*nmax, *lmax, *rank);
	}
	void irr_simd_close_(int *rank){
		irr_simd_close(*rank);
	}
    void irr_simd_profile_(int *rank){
        irr_simd_profile(*rank);
    }
    
	void irr_simd_set_jp_(
		int    *addr,
		double  pos [3],
		double  vel [3],
		double  acc2[3],
		double  jrk6[3],
		double *mass,
		double *time)
	{
		irr_simd_set_jp((*addr)-1, pos, vel, acc2, jrk6, *mass, *time);
	}
	void irr_simd_set_list_(
		int *addr,
		int *nblist)
	{
		irr_simd_set_list((*addr)-1, *nblist, nblist+1);
	}
	void irr_simd_firr_vec_(
            double *ti,
			int   *ni,
			int    addr[],
			double acc [][3],
			double jrk [][3],
            int nnbid[])
	{
      irr_simd_firr_vec(*ti, *ni, addr, acc, jrk, nnbid);
	}
}
