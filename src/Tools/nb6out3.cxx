#include "nb6out3.h"
#include <cstdio>

//Header============================================//
void header::print_all(FILE *fout) {
  fprintf(fout,"NTOT %d MODEL %d NRUN %d NK %d\n",ntot,model,nrun,nk);
}

void header::print_data(FILE *fout) {
  fprintf(fout,"%d %d %d %d\n",ntot,model,nrun,nk);
}

void header::print_title(FILE *fout) {
  fprintf(fout,"NTOT MODEL NRUN NK\n");
}

//PARTICLE6================================================//
void particle6::print_data(FILE *fout) {
  fprintf(fout,"%g %g %g %g %g %g %g %d\n",mass,x[0],x[1],x[2],v[0],v[1],v[2],name);
}

void particle6::print_title(FILE *fout) {
  fprintf(fout,"mass x[1] x[2] x[3] v[1] v[2] v[3] name\n");
}

//Pars6============================================//
void pars6::print_all(FILE *fout) {
  fprintf(fout,"Ttot:(nb) %g (tcr) %g tscale: %g rbar: %g zmbar: %g vstar: %g rscale: %g rtide: %g rdens: [0] %g [1] %g [2] %g rsmin: %g tidal: %g rc: %g nc: %g vc: %g rhom: %g cmax: %g nparis: %g dmin1: %g\n",ttot_nb,ttot_tcr,tscale,rbar,zmbar,vstar,rscale,rtide,rdens[0],rdens[1],rdens[2],rsmin,tidal,rc,nc,vc,rhom,cmax,npairs,dmin1);
}

void pars6::print_data(FILE *fout) {
  fprintf(fout,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",ttot_nb,ttot_tcr,tscale,rbar,zmbar,vstar,rscale,rtide,rdens[0],rdens[1],rdens[2],rsmin,tidal,rc,nc,vc,rhom,cmax,npairs,dmin1);
}

void pars6::print_title(FILE *fout) {
  fprintf(fout,"Ttot(nb) Ttot(tcr) tscale rbar zmbar vstar rscale rtide rdens[0] rdens[1] rdens[2] rsmin tidal rc nc vc rhom cmax nparis dmin1\n");
}

//PARTICLE6pp==============================================//
void particle6pp::print_data(FILE *fout) {
  fprintf(fout,"%g %g %g %g %g %g %g %g %g %g %d\n",mass,x[0],x[1],x[2],v[0],v[1],v[2],rhos,xns,phi,name);
}

void particle6pp::print_title(FILE *fout) {
  fprintf(fout,"mass x[1] x[2] x[3] v[1] v[2] v[3] rhos xns phi name\n");
}

//Pars6pp========================================//
void pars6pp::print_all(FILE *fout) {
  fprintf(fout,"Ttot:(nb) %g (tcr) %g tscale: %g rbar: %g zmbar: %g vstar: %g rscale: %g rtide: %g rdens: [0] %g [1] %g [2] %g rsmin: %g tidel: %g rc: %g nc: %g vc: %g rhom: %g cmax: %g nparis: %g dmin1: %g\n",ttot_nb,ttot_tcr,tscale,rbar,zmbar,vstar,rscale,rtide,rdens[0],rdens[1],rdens[2],rsmin,tidal,rc,nc,vc,rhom,cmax,npairs,dmin1);
}

void pars6pp::print_data(FILE *fout) {
  fprintf(fout,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",ttot_nb,ttot_tcr,tscale,rbar,zmbar,vstar,rscale,rtide,rdens[0],rdens[1],rdens[2],rsmin,tidal,rc,nc,vc,rhom,cmax,npairs,dmin1);
}

void pars6pp::print_title(FILE *fout) {
  fprintf(fout,"Ttot(nb) Ttot(tcr) tscale rbar zmbar vstar rscale rtide rdens[0] rdens[1] rdens[2] rsmin tidel rc nc vc rhom cmax nparis dmin1\n");
}

bool readp6(const int &NMAX, particle6 *data, header *h, pars6 *p, FILE* stream, bool fromnb6, int off){
  if(fromnb6) {
    float *tmp=new float[off?off*4:1];
    if(off) fread(tmp,sizeof(int),off,stream);
    fread(h,sizeof(header),1,stream);
    if(h->ntot>NMAX) return false;
    if(off) fread(tmp,sizeof(int),2*off,stream);
    fread(p,sizeof(pars6),1,stream);
    float *m=new float[h->ntot];
    float *x=new float[3*h->ntot];
    float *v=new float[3*h->ntot];
    int *n=new int[h->ntot];
    fread(m,sizeof(float),h->ntot,stream);
    fread(x,sizeof(float),3*h->ntot,stream);
    fread(v,sizeof(float),3*h->ntot,stream);
    fread(n,sizeof(int),h->ntot,stream);
    if(off) fread(tmp,sizeof(int),off,stream);
    for (int i=0;i<h->ntot;i++) {
      data[i].mass=m[i];
      data[i].x[0]=x[3*i];
      data[i].x[1]=x[3*i+1];
      data[i].x[2]=x[3*i+2];
      data[i].v[0]=v[3*i];
      data[i].v[1]=v[3*i+1];
      data[i].v[2]=v[3*i+2];
      data[i].name=n[i];
    }
  }
  else {
    fread(h,sizeof(header),1,stream);
    if(h->ntot>NMAX) return false;
    fread(p,sizeof(pars6),1,stream);
    fread(data,sizeof(particle6),h->ntot,stream);
  }
#ifdef DEBUG
  h->print_all(stdout);
  p->print_all(stdout);
#endif
  return true;
}

bool readp6(const int &NMAX, float *m, float3 *x, float3 *v, int *n, header *h, pars6 *p, FILE* stream, int off){
  float *tmp=new float[off?off*4:1];
  if(off) fread(tmp,sizeof(int),off,stream);
  fread(h,sizeof(header),1,stream);
  if(h->ntot>NMAX) return false;
  if(off) fread(tmp,sizeof(int),off*2,stream);
  fread(p,sizeof(pars6),1,stream);
  fread(m,sizeof(float),h->ntot,stream);
  fread(x,sizeof(float),3*h->ntot,stream);
  fread(v,sizeof(float),3*h->ntot,stream);
  fread(n,sizeof(int),h->ntot,stream);
  if(off) fread(tmp,sizeof(int),off,stream);
#ifdef DEBUG
  h->print_all(stdout);
  p->print_all(stdout);
#endif
  return true;
}

bool readp6pp(const int &NMAX, particle6pp *data, header *h, pars6pp *p, FILE* stream, bool fromnb6pp, int off){
  if(fromnb6pp) {
    float *tmp=new float[off?off*4:1];
    if(off) fread(tmp,sizeof(int),off,stream);
    fread(h,sizeof(header),1,stream);
    if(h->ntot>NMAX) return false;
    if(feof(stream)) {
      perror("Error: reach end !\n");
      return false;
    }
    if(off) fread(tmp,sizeof(int),off*2,stream);
    fread(p,sizeof(pars6pp),1,stream);
    float *m=new float[h->ntot];
    float *rh=new float[h->ntot];
    float *xns=new float[h->ntot];    
    float *x=new float[3*h->ntot];
    float *v=new float[3*h->ntot];
    float *phi=new float[h->ntot];        
    int *n=new int[h->ntot];
    fread(m,sizeof(float),h->ntot,stream);
    fread(rh,sizeof(float),h->ntot,stream);
    fread(xns,sizeof(float),h->ntot,stream);
    fread(x,sizeof(float),3*h->ntot,stream);
    fread(v,sizeof(float),3*h->ntot,stream);
    fread(phi,sizeof(float),h->ntot,stream);
    fread(n,sizeof(int),h->ntot,stream);
    if(off) fread(tmp,sizeof(int),off,stream);
    for (int i=0;i<h->ntot;i++) {
      data[i].mass=m[i];
      data[i].x[0]=x[3*i];
      data[i].x[1]=x[3*i+1];
      data[i].x[2]=x[3*i+2];
      data[i].v[0]=v[3*i];
      data[i].v[1]=v[3*i+1];
      data[i].v[2]=v[3*i+2];
      data[i].rhos=rh[i];
      data[i].xns=xns[i];
      data[i].phi=phi[i];
      data[i].name=n[i];
    }
  }
  else {
    fread(h,sizeof(header),1,stream);
    if(h->ntot>NMAX) return false;    
    fread(p,sizeof(pars6pp),1,stream);
    fread(data,sizeof(particle6pp),h->ntot,stream);
  }
#ifdef DEBUG
  h->print_all(stdout);
  p->print_all(stdout);
#endif
  return true;
}

bool readp6pp(const int &NMAX, float *m, float3 *x, float3 *v, float *rh, float *xns, float *phi, int *n, header *h, pars6pp *p, FILE* stream, int off){
  float *tmp=new float[off?off*4:1];
  if(off) fread(tmp,sizeof(int),off,stream);
  fread(h,sizeof(header),1,stream);
  if(h->ntot>NMAX) return false;
  if(off) fread(tmp,sizeof(int),off*2,stream);
  fread(p,sizeof(pars6pp),1,stream);
  fread(m,sizeof(float),h->ntot,stream);
  fread(rh,sizeof(float),h->ntot,stream);
  fread(xns,sizeof(float),h->ntot,stream);
  fread(x,sizeof(float3),h->ntot,stream);
  fread(v,sizeof(float3),h->ntot,stream);
  fread(phi,sizeof(float),h->ntot,stream);
  fread(n,sizeof(int),h->ntot,stream);
  if(off) fread(tmp,sizeof(int),off,stream);
#ifdef DEBUG
  h->print_all(stdout);
  p->print_all(stdout);
#endif
  return true;
}

bool writep6(particle6 *din, header *h, pars6 *p, FILE* fout, bool binary, bool with_title){
  if(!fout||!din||!h||!p) return false;
  if (binary) {
    fwrite(h,sizeof(header),1,fout);
    fwrite(p,sizeof(pars6pp),1,fout);
    fwrite(din,sizeof(particle6),h->ntot,fout);
  }
  else {
    if(with_title) h->print_title(fout);
    h->print_data(fout);
    if(with_title) p->print_title(fout);
    p->print_data(fout);
    if(with_title) din->print_title(fout);
    for (int i=0;i<h->ntot;i++) din[i].print_data(fout);
  }
  return true;
}

bool writep6pp(particle6pp *din, header *h, pars6pp *p, FILE* fout, bool binary, bool with_title){
  if(!fout||!din||!h||!p) return false;
  if (binary) {
    fwrite(h,sizeof(header),1,fout);
    fwrite(p,sizeof(pars6pp),1,fout);
    fwrite(din,sizeof(particle6pp),h->ntot,fout);
  }
  else {
    if(with_title) h->print_title(fout);    
    h->print_data(fout);
    if(with_title) p->print_title(fout);    
    p->print_data(fout);
    if(with_title) din->print_title(fout);
    for (int i=0;i<h->ntot;i++) din[i].print_data(fout);
  }
  return true;
}
