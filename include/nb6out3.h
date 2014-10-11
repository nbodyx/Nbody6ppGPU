//The structure for nbody6/6++ OUT3=================//
//I/O for OUT3 =====================================//
/*Structure:
  For nbody6
   header:
     4 members: ntot (all particles including center mass)
                model (index of model defined in input file)
                nrun
                nk (parameter number)
     3 member functions:
                print_all(FILE* fout): print all data with titles to file I/O fout with format
                print_data(FILE* fout): print all data with format
                print_title(FILE* fout): print title (for the first line in the data file)
   particle6:
     4 members: mass (mass in NB unit)
                x[3] (position in NB unit)
                v[3] (velocity in NB unit)
                name (particle index)
     2 member functions:
                print_data(FILE* fout) similar as header.print_data
                print_title(FILE* fout) similar as header.print_title
   pars6:
     18 members: 1 ttot_nb   total time of snapshot in Nbody unit 
                 2 npairs    number of binaries 
                 3 rbar   distance scaling factor: (PC). R(astro. unit)=RBAR * R(Nbody unit) 
                 4 zmbar  mass scaling factor: (M_sun) 
                 5 rtide  tidal radius in Nbody Unit 
                 6 tidal
                 7 rdens[3]
                 8 ttot_tcr    total time per crossing time 
                 9 tscale      time scaling factor
                 10 vstar       velocity scaling factor
                 11 rc          core radius in Nbody unit 
                 12 nc          number of stars inside core radius 
                 13 vc 
                 14 rhom
                 15 cmax
                 16 rscale
                 17 rsmin
                 18 dmin1
     3 member functions: print_all(FILE *fout)
                         print_data(FILE *fout)
                         print_title(FILE *fout)
  For nbody6++
    particle6pp:
      7 members: mass
                 x[3]
                 v[3]
                 rhos
                 xns
                 phi    potential
                 name
      2 member functions: print_data(FILE *fout)
                          print_title(FILE *fout)
    pars6pp:
     18 members: 1 ttot_nb   total time of snapshot in Nbody unit 
                 2 npairs    number of binaries 
                 3 rbar   distance scaling factor: (PC). R(astro. unit)=RBAR * R(Nbody unit) 
                 4 zmbar  mass scaling factor: (M_sun) 
                 5 rtide  tidal radius in Nbody Unit 
                 6 tidal
                 7 rdens[3]
                 8 ttot_tcr    total time per crossing time 
                 9 ttot_astro  integer time with astronomical unit
                 10 nzero
                 11 rc          core radius in Nbody unit 
                 12 nc          number of stars inside core radius 
                 13 vc 
                 14 rhom
                 15 cmax
                 16 rscale
                 17 rsmin
                 18 dmin1
      2 member functions: print_all(FILE *fout);
                          print_data(FILE *fout);
                          print_title(FILE *fout);
*/
#ifndef NB6OUT3_H
#define NB6OUT3_H

#include <cstdio>
#include <cstring>

typedef float float3[3];

struct header{
  int ntot,model,nrun,nk;
  void print_all(FILE *fout);
  void print_data(FILE *fout);
  void print_title(FILE *fout);
};

//For nbody6========================================//
struct particle6{
  float mass,x[3],v[3];
  int name;
  void print_data(FILE *fout);
  void print_title(FILE *fout);
  particle6& operator=(const particle6& a) {
    memcpy(this,&a,sizeof(particle6));
    return *this;
  }
};
  
struct pars6{
  float ttot_nb,npairs,rbar,zmbar,rtide,tidal,rdens[3],ttot_tcr,tscale,vstar,rc,nc,vc,rhom,cmax,rscale,rsmin,dmin1;
  void print_all(FILE *fout);
  void print_data(FILE *fout);
  void print_title(FILE *fout);
};

//For nbody6++======================================//
struct particle6pp{
  float mass,x[3],v[3],rhos,xns,phi;
  int name;
  void print_data(FILE *fout);
  void print_title(FILE *fout);
  particle6pp& operator=(const particle6pp& a) {
    memcpy(this,&a,sizeof(particle6pp));
    return *this;
  }
};

struct pars6pp{
  float ttot_nb,npairs,rbar,zmbar,rtide,tidal,rdens[3],ttot_tcr,tscale,vstar,rc,nc,vc,rhom,cmax,rscale,rsmin,dmin1;
  void print_all(FILE *fout);
  void print_data(FILE *fout);
  void print_title(FILE *fout);
};

//READING ==========================================//
bool readp6(const int &NMAX, particle6 *data, header *h, pars6 *p, FILE* fin, bool fromnb6=true, int offset=1);
// NMAX: maximum particle number can be read.

bool readp6(const int &NMAX, float *m, float3 x[], float3 v[], int *n, header *h, pars6 *p, FILE* fin, int offset=1);

bool readp6pp(const int &NMAX, particle6pp *data, header *h, pars6pp *p, FILE* fin, bool fromnb6pp=true, int offset=1);

bool readp6pp(const int &NMAX, float *m, float3 x[], float3 v[], float *rh, float *xns, float*phi, int *n, header *h, pars6pp *p, FILE* fin, int offset=1);

//WRITING===========================================//
bool writep6(particle6 *din, header *h, pars6 *p, FILE* fout, bool binary=true,bool with_title=false);

bool writep6pp(particle6pp *din, header *h, pars6pp *p, FILE* fout, bool binary=true, bool with_title=false);

#endif
