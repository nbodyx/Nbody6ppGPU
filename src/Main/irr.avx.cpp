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

static inline v8sf rsqrt_NR(const v8sf x){
	const v8sf y = __builtin_ia32_rsqrtps256(x);
	const v8sf c1 = REP8(-0.5f);
	const v8sf c2 = REP8(-3.0f);
	return (c1 * y) * (x*y*y + c2);
}

static inline v8sf pack_2_v4sf(const v4sf a, const v4sf b){
	// v8sf p;
	v8sf p = REP8(0.0f); // just avoid warning
	p = __builtin_ia32_vinsertf128_ps256(p, a, 0);
	p = __builtin_ia32_vinsertf128_ps256(p, b, 1);
	return p;
}

// more efficieent?
static inline v8sf pack_2_v4sf(const v4sf a, const v4sf b, v8sf &p){
	p = __builtin_ia32_vinsertf128_ps256(p, a, 0);
	p = __builtin_ia32_vinsertf128_ps256(p, b, 1);
	return p;
}

static inline v4df pack_2_v2df(const v2df a, const v2df b){
	v4df p = REP4(0.0f); // just avoid warning
	p = __builtin_ia32_vinsertf128_pd256(p, a, 0);
	p = __builtin_ia32_vinsertf128_pd256(p, b, 1);
	return p;
}

struct Particle{
	v4sf posH; // elemen w is for mass
	v4sf posL;
	v4sf vel;
	v4sf acc2;
	v4sf jrk6;
	v2df time; // 6 xmm words;

	void make_pos(const double *x, const double *m){
		const v4df xd = {x[0], x[1], x[2], m[0]};
		const v4sf xs = __builtin_ia32_cvtpd2ps256(xd);
		const v4df yd = xd - __builtin_ia32_cvtps2pd256(xs);
		const v4sf ys = __builtin_ia32_cvtpd2ps256(yd);
		posH = xs;
		posL = ys;
	}
	static v4sf make_v4sf(const double *x){
		const v4df xd = {x[0], x[1], x[2], 0.0};
		const v4sf xs = __builtin_ia32_cvtpd2ps256(xd);
		return xs;
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
		const v8sf ym0 = pack_2_v4sf(posH, posL);
		const v8sf ym1 = pack_2_v4sf(vel,  acc2);
		const v8sf ym2 = pack_2_v4sf(jrk6, (v4sf)time);
		v8sf *dst = (v8sf *)pdst;
		dst[0] = ym0;
		dst[1] = ym1;
		dst[2] = ym2;
	}
	void sstore(Particle *pdst) const {
		const v8sf ym0 = pack_2_v4sf(posH, posL);
		const v8sf ym1 = pack_2_v4sf(vel,  acc2);
		const v8sf ym2 = pack_2_v4sf(jrk6, (v4sf)time);
		v8sf *dst = (v8sf *)pdst;
		__builtin_ia32_movntps256((float *)(dst + 0), ym0);
		__builtin_ia32_movntps256((float *)(dst + 1), ym1);
		__builtin_ia32_movntps256((float *)(dst + 2), ym2);
	}
};

struct Pred2{ // couple of 2 predictors
	v8sf posH; // |m|z|y|x|m|z|y|x|
	v8sf posL; // |*|z|y|x|*|z|y|x|
	v8sf vel;  // |*|w|v|u|*|w|v|u|

	static inline v8sf gen_dt01(const v2df t0, const v2df t1, const v2df tnow){
		const v2df dt01d = tnow - __builtin_ia32_unpcklpd(t0, t1);
		const v4sf dt01s = __builtin_ia32_cvtpd2ps(dt01d);
		const v4sf dt0   = __builtin_ia32_shufps(dt01s, dt01s, 0x00);
		const v4sf dt1   = __builtin_ia32_shufps(dt01s, dt01s, 0x55);
		return pack_2_v4sf(dt0, dt1);
	}

	Pred2(const Particle &p0, const Particle &p1, const v2df tnow){
		const v8sf s0 = gen_dt01(p0.time, p1.time, tnow);
		const v8sf s1 = s0 + s0;
		const v8sf s2 = s0 * (v8sf)REP8(1.5f);

		const v8sf pos8 = pack_2_v4sf(p0.posL, p1.posL);
		const v8sf vel8 = pack_2_v4sf(p0.vel , p1.vel );
		const v8sf acc8 = pack_2_v4sf(p0.acc2, p1.acc2);
		const v8sf jrk8 = pack_2_v4sf(p0.jrk6, p1.jrk6);

		this->posH = pack_2_v4sf(p0.posH, p1.posH);
		// this->posH = pack_2_v4sf(p0.posH, p1.posH, this->posH);
		this->posL = pos8 + s0*(vel8 + s0*(acc8 + s0*(jrk8)));
		this->vel  = vel8 + s1*(acc8 + s2*(jrk8));
	}
};

struct Pred8{ // almost an array of structures
	v8sf xH, yH, zH, m;
	v8sf xL, yL, zL;
	v8sf vx, vy, vz;

	static void transpos(v8sf &v0, v8sf &v1, v8sf &v2, v8sf &v3){
		// transpose from
		// [ |w4|z4|y4|x4||w0|z0|y0|x0|, |w5|z5|y5|x5||w1|z1|y1|x1|, |w6|z6|y6|x6||w2|z2|y2|x2|, |w7|z7|y7|x7||w3|z3|y3|x3| ]
		// to
		// [ |x7|x6|x5|x4||x3|x2|x1|x0|, |y7|y6|y5|y4||y3|y2|y1|y0|, |z7|z6|z5|z4||z3|z2|z1|z0|, |w7|w6|w5|w4||w3|w2|w1|w0| ] 
		const v8sf y2y0x2x0 = __builtin_ia32_unpcklps256(v0, v2);
		const v8sf w2w0z2z0 = __builtin_ia32_unpckhps256(v0, v2);
		const v8sf y3y1x3x1 = __builtin_ia32_unpcklps256(v1, v3);
		const v8sf w3w1z3z1 = __builtin_ia32_unpckhps256(v1, v3);
		const v8sf xxxx = __builtin_ia32_unpcklps256(y2y0x2x0, y3y1x3x1);
		const v8sf yyyy = __builtin_ia32_unpckhps256(y2y0x2x0, y3y1x3x1);
		const v8sf zzzz = __builtin_ia32_unpcklps256(w2w0z2z0, w3w1z3z1);
		const v8sf wwww = __builtin_ia32_unpckhps256(w2w0z2z0, w3w1z3z1);
		v0 = xxxx;
		v1 = yyyy;
		v2 = zzzz;
		v3 = wwww;
	}

	Pred8(){}

	// for j-particle
	Pred8(const Pred2 &p0, const Pred2 &p1, const Pred2 &p2, const Pred2 &p3){
		{
			v8sf a0 = p0.posH;
			v8sf a1 = p1.posH;
			v8sf a2 = p2.posH;
			v8sf a3 = p3.posH;
			transpos(a0, a1, a2, a3);
			xH = a0;
			yH = a1;
			zH = a2;
			m  = a3;
		}
		{
			v8sf b0 = p0.posL;
			v8sf b1 = p1.posL;
			v8sf b2 = p2.posL;
			v8sf b3 = p3.posL;
			transpos(b0, b1, b2, b3);
			xL = b0;
			yL = b1;
			zL = b2;
		}
		{
			v8sf c0 = p0.vel;
			v8sf c1 = p1.vel;
			v8sf c2 = p2.vel;
			v8sf c3 = p3.vel;
			transpos(c0, c1, c2, c3);
			vx = c0;
			vy = c1;
			vz = c2;
		}
	}
	// for j-particle
	Pred8(const Pred2 &pr){
		xH = __builtin_ia32_shufps256(pr.posH, pr.posH, 0x00);
		yH = __builtin_ia32_shufps256(pr.posH, pr.posH, 0x55);
		zH = __builtin_ia32_shufps256(pr.posH, pr.posH, 0xaa);

		xL = __builtin_ia32_shufps256(pr.posL, pr.posL, 0x00);
		yL = __builtin_ia32_shufps256(pr.posL, pr.posL, 0x55);
		zL = __builtin_ia32_shufps256(pr.posL, pr.posL, 0xaa);

		vx = __builtin_ia32_shufps256(pr.vel,  pr.vel , 0x00);
		vy = __builtin_ia32_shufps256(pr.vel,  pr.vel , 0x55);
		vz = __builtin_ia32_shufps256(pr.vel,  pr.vel , 0xaa);
	}
};

// single precision version
struct Force{
	v8sf ax, ay, az;
	v8sf jx, jy, jz;
	v8sf vnnb;

	void clear(){
		const v8sf zero = REP8(0.0f);
		ax = ay = az = zero;
		jx = jy = jz = zero;

		v4df tmp = {HUGE, HUGE, HUGE, HUGE};
		vnnb = (v8sf)tmp;
	}
	static double reduce(const v8sf v){
		const v4df d0 = __builtin_ia32_cvtps2pd256(
		                  __builtin_ia32_vextractf128_ps256(v, 0));
		const v4df d1 = __builtin_ia32_cvtps2pd256(
		                  __builtin_ia32_vextractf128_ps256(v, 1));
		const v4df dd = d0 + d1;
		const v4df dh = __builtin_ia32_haddpd256(dd, dd);
		const v2df sum = __builtin_ia32_vextractf128_pd256(dh, 0)
		               + __builtin_ia32_vextractf128_pd256(dh, 1);
#ifdef __USE_INTEL
        double p;
        _mm_store_sd(&p, sum);
        return p;
#else
		return __builtin_ia32_vec_ext_v2df(sum, 0);
#endif
	}

	int reduce_nnb() const{
		v2df min = __builtin_ia32_minpd(
				__builtin_ia32_vextractf128_pd256((v4df)vnnb, 0),
				__builtin_ia32_vextractf128_pd256((v4df)vnnb, 1));
		min = __builtin_ia32_minpd(min,
				__builtin_ia32_unpckhpd(min, min));
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
    void calc_and_accum(const Pred8 &pi, const Pred8 &pj, const v8sf idx){
		const v8sf dx = (pj.xH - pi.xH) + (pj.xL - pi.xL);
		const v8sf dy = (pj.yH - pi.yH) + (pj.yL - pi.yL);
		const v8sf dz = (pj.zH - pi.zH) + (pj.zL - pi.zL);

		const v8sf dvx = pj.vx - pi.vx;
		const v8sf dvy = pj.vy - pi.vy;
		const v8sf dvz = pj.vz - pi.vz;

		const v8sf r2 = dx*dx  + dy*dy  + dz*dz;
		const v8sf rv = dx*dvx + dy*dvy + dz*dvz;
		const v8sf rinv   = rsqrt_NR(r2);
		const v8sf rinv2  = rinv * rinv;
		const v8sf c1     = REP8(-3.0f);
		const v8sf alpha  = c1 * rinv2 * rv;
		const v8sf mrinv3 = pj.m * rinv * rinv2;

		const v8sf r2_idx0 = __builtin_ia32_unpcklps256(idx, r2);
		const v8sf r2_idx1 = __builtin_ia32_unpckhps256(idx, r2);
		vnnb = (v8sf)__builtin_ia32_minpd256((v4df)vnnb, (v4df)r2_idx0);
		vnnb = (v8sf)__builtin_ia32_minpd256((v4df)vnnb, (v4df)r2_idx1);

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

	fprintf(stderr, "# Opening IRR_simd lib. AVX ver. - rank: %d; nmax: %d, lmax: %d\n",rank,nmax,lmax); 
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

	fprintf(stderr, "[R.%d AVX Irr.F ] Ncall: %llu <NI>: %d <NB>: %f grav: %f s, %f Gflops, %f usec\n",
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

#ifdef __USE_INTEL
     const int one4[4] = REP4(1);
     __m128i one;
     _mm_stream_si128(&one, *((__m128i*) one4));
    //    const __m128i one =  _mm_stream_load_si128 ((__m128i*) one4 );
#endif
	const int *src = nblist;
	      int *dst = list[addr].nb;
	for(int k=0; k<nnb; k+=8){
      // assert((unsigned long)dst %16 == 0);
#ifdef __USE_INTEL
      //        const __m128i one  = 
      //        const __m128i idx0 = _mm_stream_load_si128 ((__m128i*) (src+k+0));
      //        const __m128i idx1 = _mm_stream_load_si128 ((__m128i*) (src+k+4));
        __m128i idx0,idx1;
        _mm_stream_si128(&idx0, *((__m128i*) (src+k+0)));
        _mm_stream_si128(&idx1, *((__m128i*) (src+k+4)));
  
        _mm_stream_si128((__m128i *) (dst+k+0), idx0-one);
        _mm_stream_si128((__m128i *) (dst+k+4), idx1-one);
#else
        typedef long long v2di __attribute__ ((__vector_size__ (16)));
        typedef int       v4si __attribute__((vector_size(16)));
		const v4si one = REP4(1);
		const v4si idx0 = (v4si)__builtin_ia32_loaddqu((const char *)(src+k+0));
		const v4si idx1 = (v4si)__builtin_ia32_loaddqu((const char *)(src+k+4));
		__builtin_ia32_movntdq((v2di *)(dst+k+0), (v2di)(idx0-one));
		__builtin_ia32_movntdq((v2di *)(dst+k+4), (v2di)(idx1-one));
#endif
	}
	// fill dummy
	const int nmax = ::nmax;
	const int kmax = 8 * (1 + (nnb-1)/8);
	for(int k=nnb; k<kmax; k++){
		dst[k] = nmax;
	}

}

static void irr_simd_firr(
		const int addr,
		_out_ double accout[3],
		_out_ double jrkout[3],
        _out_ int  &nnbid)
{
	const Particle *ptcl = ::ptcl;
	const v2df      tnow = ::vec_tnow;
	const int  nnb  = list[addr].nnb;

	Force force;
	force.clear();

	const Pred8 pri(Pred2(ptcl[addr], ptcl[addr], tnow));

	for(int k=0; k<nnb; k++){
		const int jaddr = list[addr].nb[k];
		ptcl[jaddr].prefetch();
	}

	for(int k=0; k<nnb; k+=8){
		const int *jptr = &(list[addr].nb[k]);
		const Pred2 p04(ptcl[jptr[0]], ptcl[jptr[4]], tnow);
		const Pred2 p15(ptcl[jptr[1]], ptcl[jptr[5]], tnow);
		const Pred2 p26(ptcl[jptr[2]], ptcl[jptr[6]], tnow);
		const Pred2 p37(ptcl[jptr[3]], ptcl[jptr[7]], tnow);

		const Pred8 prj(p04, p15, p26, p37);
		// const Pred2 p01(ptcl[jptr[0]], ptcl[jptr[1]], tnow);
		// const Pred2 p23(ptcl[jptr[2]], ptcl[jptr[3]], tnow);
		// const Pred2 p45(ptcl[jptr[4]], ptcl[jptr[5]], tnow);
		// const Pred2 p67(ptcl[jptr[6]], ptcl[jptr[7]], tnow);

		// const Pred8 prj(p01, p23, p45, p67);

		assert(0 == (size_t)jptr % 32);
		v8sf idx = *(v8sf *)jptr;
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
            int    nnbid[])
	{
      irr_simd_firr_vec(*ti, *ni, addr, acc, jrk, nnbid);
	}
}
