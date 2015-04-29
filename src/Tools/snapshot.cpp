//transform bdat.9, bwdat.19, sev.83, bew.82 and conf.3 to specific format//
//
// Output option 1:
// x1(1:3)[pc], v1(1:3)[km/s], x2(1:3)[pc], v2(1:3)[km/s], m1[M_sun], m2[M_sun], log10[L1[L_sun]), log10(L2[L_sun]), rs1[R_sun], rs2[R_sun], K1, K2, a[R_sun], eccentricity, xcm(1:3)[pc], vcm(1:3)[km/s]
// Here L1/2 is luminosity, rs1/2 is stellar radius, a is semi-major axis.
//   In binary case, x1/2,v1/2 are the two components positions and velocities, xcm, vcm are the center-of-mass positions and velocities.
//   In single case, x1 and xcm are same, v1 and vcm are same. the data for *2 are zero
//
// Output option 2:
//   single star data: mass[M_sun], x(1:3)[NB],v(1:3)[NB]
//   binary star data: eccentricity, log10(a[R_sun]), mass1[M_sun], mass2[M_sun], xcm(1:3)[NB], vcm(1:3)[NB]

#include <stdio.h>
#include <cstring>
#include <initial.h>
#include <nb6out3.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <uftools.h>
#include <string>
#include <algorithm>

//new data structure================================//
struct ndata{
  double x1[3],v1[3],x2[3],v2[3],m1,m2,L1,L2,rs1,rs2;  //x/v(1/2): binary component(1/2) or single star(1)[NB], m1/2: mass[M_sun], L1/2: log10 luminocity[L_sun], rs1/2: star radius[R_sun]
  int k1,k2; //k1/2: stellar type
  double a,e,x[3],v[3]; //a: semi-major[R_sun], e: eccentricity, x/v: binary c.m. or single star[NB]
  ndata() {}
  ndata(float a[3],float b[3],float c[3],float d[3],double p1,double p2,double p3,double p4,double p5,double p6,int p7,int p8,double p9,double p10,double p11[3],double p12[3]): m1(p1),m2(p2),L1(p3),L2(p4),rs1(p5),rs2(p6),k1(p7),k2(p8),a(p9),e(p10) {
    x1[0]=a[0];    x1[1]=a[1];    x1[2]=a[2];    v1[0]=b[0];    v1[1]=b[1];    v1[2]=b[2];
    x2[0]=c[0];    x2[1]=c[1];    x2[2]=c[2];    v2[0]=d[0];    v2[1]=d[1];    v2[2]=d[2];
    x[0]=p11[0];   x[1]=p11[1];   x[2]=p11[2];   v[0]=p12[0];   v[1]=p12[1];   v[2]=p12[2];
  }
  ndata(float a[3],float b[3],float c[3],float d[3],double p1,double p2,double p3,double p4,double p5,double p6,int p7,int p8,double p9,double p10,float p11[3],float p12[3]): m1(p1),m2(p2),L1(p3),L2(p4),rs1(p5),rs2(p6),k1(p7),k2(p8),a(p9),e(p10) {
    x1[0]=a[0];    x1[1]=a[1];    x1[2]=a[2];    v1[0]=b[0];    v1[1]=b[1];    v1[2]=b[2];
    x2[0]=c[0];    x2[1]=c[1];    x2[2]=c[2];    v2[0]=d[0];    v2[1]=d[1];    v2[2]=d[2];
    x[0]=p11[0];   x[1]=p11[1];   x[2]=p11[2];   v[0]=p12[0];   v[1]=p12[1];   v[2]=p12[2];
  }
};

//sev data structure================================//
struct sev{
  double time;
  int i,name,k;
  double ri,m,logL,logR,Temp;
  sev& operator=(const sev& a) {
    memcpy(this,&a,sizeof(sev));
    return *this;
  }
};

//bev data structure================================//  
struct bev{
  double time;
  int i1,i2,name,name2,k1,k2,kcm;
  double ri,ecc,P,a,m1,m2,logL1,logL2,logR1,logR2,Temp1,Temp2;
  bev& operator=(const bev& a) {
    memcpy(this,&a,sizeof(bev));
    return *this;
  }
};

//Get module square of vector==========================//
template <typename T>
T dot (T x[3]) {
  return x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
}

//Sorting function==================================//
bool ismaller_r(ndata a, ndata b) {
  return dot(a.x)<dot(b.x);
}

//Lagrangian radii fraction==========================//
const float fraction[18]={0.001,0.003,0.005,0.01,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,1.0};

int main(int argc, char *argv[]){
    
  //initial arguments=================================//
  pars_initial init(".snapshot_config");
  init.add("time","time of snapshot file","0.0");
  init.add("out","output file name","data");
  init.add("off","offset for first line",(int)1);
  init.add("nmax","maximum particle number",(int)1000000);
  init.add("semi","semi-major axis[AU] upper limit",(double)50);
  init.add("format","1: Snapshot for projection, generate lagrangian radii, 2: Snapshot for MOCCA",(int)1);
  init.add("bin","0: without binaries, 1: with binaries",(int)1);
  init.add("datapath","nbody data path","./");
  init.initial(argc,argv);

  std::string time=init.gets("time");
  const bool bflag=init.geti("bin");
  std::string datapath=init.gets("datapath");

  //m,x,v data========================================//
  FILE *sin;
  if ( (sin = fopen((datapath+"/conf.3_"+time).c_str(),"r")) == NULL) {
    fprintf(stderr,"Error: Cannot open input file %s.\n",("conf.3_"+time).c_str());
    return 0;
  }

  //Maximum particle number===========================//
  const int NMAX=init.geti("nmax");

  //conf.3 data ======================================//
  particle6pp *datain=new particle6pp[NMAX];
  particle6pp *data=new particle6pp[NMAX];
  memset(data,0,NMAX*sizeof(particle6pp));
  header *h=new header;
  pars6pp *par=new pars6pp;
  int off=init.geti("off");

  //reading conf.3====================================//
  printf("Reading m,x,v data (conf.3_%s)\n",time.c_str());
  if(!readp6pp(NMAX,datain,h,par,sin,true,off)) {
    perror("Error: NMAX too small!\n");
    return 0;
  }

  int n=(int)h->ntot;            //Ntot singles + 3*binaries
  int nt=n-(int)par->npairs;     //N singles + 2*binaries
  int ifirst=(int)par->npairs*2; //first single star position
  int ns=nt-ifirst+(int)par->npairs;       //N singles + binaries

  //locate particle data by name as array index=======//
  printf("Re-aligning particle data\n");
  int smin=NMAX,smax=0;
  for (int i=0;i<n;i++) {
    int name=datain[i].name;
    if(name>=NMAX) {
      fprintf(stderr,"Error!: name %d > NMAX.\n",name);
      return 0;
    }
    else if(name<=0) {
      fprintf(stdout,"Ingore special name: %d\n",name);
      continue;
    }
    if(i>=nt&&!bflag) break;
    if(smin>name) smin = name;
    if(smax<name) smax = name;
    data[name]=datain[i];
  }
  fclose(sin);

  //Single star evolution data========================//
  FILE *sein;
  if ( (sein = fopen((datapath+"/sev.83_"+time).c_str(),"r")) == NULL) {
    fprintf(stderr,"Error: Cannot open input file %s.\n.",("sev.83_"+time).c_str());
    return 0;
  }

  //reading sev.83====================================//
  //escape first line=================================//
  printf("Reading SSE. data (sev.83_%s)\n",time.c_str());
  for (int i=0;i<1;i++) {
    char c;
    do {
      c = fgetc(sein);
      //      printf("%c",c);
    }
    while (c != '\n');
  } 
  sev *sedat=new sev[smax+10];
  memset(sedat,0,(smax+10)*sizeof(sev));
  int sei=0;
  while(1) {
    sev tmp;
    fscanf(sein,"%lg %d %d %d %lg %lg %lg %lg %lg",&tmp.time,&tmp.i,&tmp.name,&tmp.k,&tmp.ri,&tmp.m,&tmp.logL,&tmp.logR,&tmp.Temp);
    if (feof(sein)) break;
    if(tmp.name>=smax+10) {
      fprintf(stderr,"Error!: name %d > Name(MAX) %d.\n",tmp.name,smax+10);
      return 0;
    }
    else if(tmp.name<=0) {
      fprintf(stdout,"Ingore special name: %d\n",tmp.name);
      continue;
    }
    sedat[tmp.name]=tmp;
    sei++;
  }
  fclose(sein);
  //  printf("SEV0: T %lg I %d\n",sedat[0].time,sedat[0].i);
  //  std::sort(sedat,&sedat[sei],ismaller_name_sev);

  //filter for binary-true:binary,false:single========//
  bool *mask=new bool[smax+1];
  for(int i=0;i<smax+1;i++) mask[i]=false;

  //output data=======================================//
  ndata *dat=new ndata[ns];

  //counter===========================================//
  int bcount=0,bwcount=0,scount=0; //N
  double mstot=0,mbtot=0;          //mass

  if(bflag) {
    //KS binary star evolution data=====================//
    FILE *bein;
    if ( (bein = fopen((datapath+"/bev.82_"+time).c_str(),"r")) == NULL) {
      fprintf(stderr,"Error: Cannot open input file %s.\n.",("bev.82_"+time).c_str());
      return 0;
    }

    //read bev.82=======================================//
    //escape first line=================================//
    printf("Reading BSE. data (bev.82_%s)\n",time.c_str());
    for (int i=0;i<1;i++) {
      char c;
      do {
        c = fgetc(bein);
        //      printf("%c",c);
      }
      while (c != '\n');
    } 
    bev *bedat=new bev[smax+10];
    memset(bedat,0,smax+10*sizeof(bev));
    int bei=0;
    while(1) {
      bev tmp;
      fscanf(bein,"%lg %d %d %d %d %d %d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",&tmp.time,&tmp.i1,&tmp.i2,&tmp.name,&tmp.name2,&tmp.k1,&tmp.k2,&tmp.kcm,&tmp.ri,&tmp.ecc,&tmp.P,&tmp.a,&tmp.m1,&tmp.m2,&tmp.logL1,&tmp.logL2,&tmp.logR1,&tmp.logR2,&tmp.Temp1,&tmp.Temp2);
      if (feof(bein)) break;
      if(tmp.name>=smax+10) {
        fprintf(stderr,"Error!: name %d > Name(MAX) %d.\n",tmp.name,smax+10);
        return 0;
      }
      else if(tmp.name<=0) {
        fprintf(stdout,"Ingore special name: %d\n",tmp.name);
        continue;
      }
      bedat[tmp.name]=tmp;
      bei++;
    }
    fclose(bein);
    //  printf("BEV0: T %lg I %d\n",bedat[0].time,bedat[0].i1);

    //binary data=======================================//
    int name[2],k1,k2,ncm,kcm;
    double m[2],E,ecc,p,a,rcm,vcm,zn,rp,step1,ecm;

    //KS binary data====================================//
    FILE *bin;
    if ( (bin = fopen((datapath+"/bdat.9_"+time).c_str(),"r")) == NULL) {
      fprintf(stderr,"Error: Cannot open input file %s.\n.",("bdat.9_"+time).c_str());
      return 0;
    }
    //escape first 4 lines==============================//
    printf("Scanning KS binary data (bdat.9_%s)\n",time.c_str());
    for (int i=0;i<4;i++) {
      char c;
      do {
        c = fgetc(bin);
        //      printf("%c",c);
      }
      while (c != '\n');
    }
    while(!feof(bin)) {
      fscanf(bin,"%d %d %lg %lg %lg %lg %lg %lg %lg %lg %d %d %lg %lg %lg %d %lg %d",&name[0],&name[1],&m[0],&m[1],&E,&ecc,&p,&a,&rcm,&vcm,&k1,&k2,&zn,&rp,&step1,&ncm,&ecm,&kcm);
      if (feof(bin)) break;
      if (bedat[name[0]].name!=name[0]||bedat[name[0]].name2!=name[1]) {
        fprintf(stderr,"Warning! bev name[0] = %d, name[1] = %d have empty data %d %d\n",name[0],name[1],bedat[name[0]].name,bedat[name[1]].name);
        //        return 0;
        continue;
      }
      if (data[name[0]].name!=name[0]||data[name[1]].name!=name[1]||data[ncm].name!=ncm) {
        fprintf(stderr,"Error, conf.3 name[0] = %d, name[1] = %d, ncm= %d have empty data %d %d %d\n",name[0],name[1],ncm,data[name[0]].name,data[name[1]].name,data[ncm].name);
        return 0;
      }
      mask[ncm]=true;
      mask[name[0]]=true;
      mask[name[1]]=true;
      dat[bcount] =
        ndata(data[name[0]].x, data[name[0]].v, data[name[1]].x, data[name[1]].v, //x1,v1,x2,v2 
              m[0], m[1], bedat[name[0]].logL1, bedat[name[0]].logL2,     //mass1,mass2,logL1,logL2
              pow(10,bedat[name[0]].logR1), pow(10,bedat[name[0]].logR2), //Rs1, Rs2
              k1, k2, a*215.095, ecc, data[ncm].x, data[ncm].v);          //K1,K2,a,e,xcm,vcm      
      bcount++;
      mbtot +=m[0]+m[1];
    }
    fclose(bin);

    //Diagnostic========================================//
    // for (int i=nt;i<n;i++) if(!mask[i]) fprintf(stderr,"Warning!: missing cm index %d name %d\n",i,data[i].name);

    //Wide binary data==================================//
    FILE *bwin;
    if ( (bwin = fopen((datapath+"/bwdat.19_"+time).c_str(),"r")) == NULL) {
      fprintf(stderr,"Error: Cannot open input file %s.\n.",("bwdat.19_"+time).c_str());
      return 0;
    }
    //reading bwdat.19==================================//
    //escape first two line=============================//
    printf("Scanning wide binary data (bwdat.19_%s)\n",time.c_str());
    for (int i=0;i<2;i++) {
      char c;
      do {
        c = fgetc(bwin);
        //      printf("%c",c);
      }
      while (c != '\n');
    } 
    while(!feof(bwin)) {
      fscanf(bwin,"%d %d %lg %lg %lg %lg %lg %lg %lg %lg %d %d",&name[0],&name[1],&m[0],&m[1],&E,&ecc,&p,&a,&rcm,&vcm,&k1,&k2);
      if (feof(bwin)) break;
      if (a>init.getd("semi")) continue;
      if (data[name[0]].name!=name[0]||data[name[1]].name!=name[1]) {
        fprintf(stderr,"Error, conf.3 name[0] = %d, name[1] = %d, ncm= %d have empty data %d %d %d\n",name[0],name[1],ncm,data[name[0]].name,data[name[1]].name,data[ncm].name);
        return 0;
      }
      if (sedat[name[0]].name!=name[0]||sedat[name[1]].name!=name[1]) {
        fprintf(stderr,"Ignore sev name[0] = %d, name[1] = %d\n",name[0],name[1]);
        continue;
      }
      mask[name[0]]=true;
      mask[name[1]]=true;
      double x[3],v[3];
      double m1=data[name[0]].mass;
      double m2=data[name[1]].mass;
      mbtot +=m[0]+m[1];
      for (int i=0;i<3;i++) {
        x[i]=(m1*data[name[0]].x[i] + m2*data[name[1]].x[i])/(m1+m2);
        v[i]=(m1*data[name[0]].v[i] + m2*data[name[1]].v[i])/(m1+m2);
      }
      dat[bcount+bwcount] =
        ndata(data[name[0]].x, data[name[0]].v, data[name[1]].x, data[name[1]].v,//x1,v1,x2,v2
              m[0], m[1], sedat[name[0]].logL, sedat[name[1]].logL,     //mass1,mass2,logL1,logL2 
              pow(10,sedat[name[0]].logR), pow(10,sedat[name[1]].logR), //Rs1, Rs2
              k1, k2, a*215.095, ecc, x, v);                            //K1,K2,a,e,xcm,vcm
      bwcount++;                                                          
    }
    fclose(bwin);
  }
  
  //Single star data==================================//
  printf("Scanning single star data. Name(min) = %d, Name(max) = %d\n",smin,smax);
  for(int i=smin;i<=smax;i++) {
    if(!mask[i]&&data[i].name==i){
      if(sedat[i].name!=i) {
        fprintf(stderr,"Warning!: sev data with name = %d empty\n",i);
        continue;
      }
      mstot +=data[i].mass*par->zmbar;
      float xzero[3]={};
      float vzero[3]={};
      dat[bcount+bwcount+scount] =
        ndata(data[i].x, data[i].v, xzero, vzero, //x1,v1,x2,v2
              data[i].mass*par->zmbar, 0.,  //mass1 mass2
              sedat[i].logL, 0.,           //logL1 logL2
              pow(10,sedat[i].logR), 0.,   //Rs1, Rs2
              sedat[i].k, 0,  0.,  0.,     //K1,K2,a,e
              data[i].x, data[i].v);        //xcm,vcm
      scount++;
    }
  }

  //counter===========================================//
  int ntt=bcount+bwcount+scount;
  int nb=bcount+bwcount;

  if(init.geti("format")==1) {
    //output IO=========================================//
    FILE *out;
    FILE *lout;
    if ( (out = fopen((init.gets("out")+"_"+time).c_str(),"w")) == NULL) {
      fprintf(stderr,"Error: Cannot open input file %s.\n",(init.gets("out")+"_"+time).c_str());
      return 0;
    }
    if ( (lout = fopen((init.gets("out")+"_lagr").c_str(),"a")) == NULL) {
      fprintf(stdout,"No %s file found, try to create new\n",(init.gets("out")+"_lagr").c_str());
      if ( (lout = fopen((init.gets("out")+"_lagr").c_str(),"w")) == NULL) {
        fprintf(stderr,"Error: Cannot open input file %s.\n",(init.gets("out")+"_lagr").c_str());
        return 0;
      }
      fprintf(lout,"Time[NB] M_tot[M*] M_stot[M*] M_btot[M*]");
      for (int i=0;i<18;i++)
        fprintf(lout," %g",fraction[i]);
      fprintf(lout,"\n");
    }

    //Sort data by distance to cluster center===========//
    printf("Sorting data by distance to cluster center\n");
    std::sort(dat,&dat[ntt],ismaller_r);

    printf("Output data\n");
    double countmass=0.;
    double cbmass=0.;
    int bicount=0;
    int lagri=0;
    double mttot=mstot+mbtot;
    double nbfrac[18]={};
    double mbfrac[18]={};
    double rlagr[18]={};
    int lcount[18]={};
    for(int i=0;i<ntt;i++) {
      fprintf(out,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d %d %lg %lg %lg %lg %lg %lg %lg %lg\n",
              dat[i].x1[0]*par->rbar, dat[i].x1[1]*par->rbar, dat[i].x1[2]*par->rbar,
              dat[i].v1[0]*par->vstar, dat[i].v1[1]*par->vstar, dat[i].v1[2]*par->vstar,
              dat[i].x2[0]*par->rbar, dat[i].x2[1]*par->rbar, dat[i].x2[2]*par->rbar,
              dat[i].v2[0]*par->vstar, dat[i].v2[1]*par->vstar, dat[i].v2[2]*par->vstar,
              dat[i].m1, dat[i].m2, dat[i].L1, dat[i].L2, dat[i].rs1, dat[i].rs2,
              dat[i].k1, dat[i].k2, dat[i].a, dat[i].e,
              dat[i].x[0]*par->rbar, dat[i].x[1]*par->rbar, dat[i].x[2]*par->rbar,
              dat[i].v[0]*par->vstar, dat[i].v[1]*par->vstar, dat[i].v[2]*par->vstar);
      countmass +=dat[i].m1+dat[i].m2;
      if(dat[i].m2>0.){
        cbmass +=dat[i].m1+dat[i].m2;
        bicount++;
      }
      if(countmass>fraction[lagri]*mttot) {
        nbfrac[lagri] = (float)bicount/(float)i;
        mbfrac[lagri] = cbmass/countmass;
        rlagr[lagri] = std::sqrt(dot(dat[i].x));
        lcount[lagri] = i;
        lagri++;
        if(lagri>18) printf("Warning!: Lagrangian radii counter larger than 18!, i=%d, mttot=%lg, countmass=%lg\n",i,mttot,countmass);
      }
    }
    fprintf(lout,"%s %lg %lg %lg",time.c_str(),mttot,mstot,mbtot);
    for(int i=0;i<18;i++) {
      fprintf(lout," %lg",rlagr[i]);
      fprintf(lout," %d",lcount[i]);
      fprintf(lout," %lg",nbfrac[i]);
      fprintf(lout," %lg",mbfrac[i]);
    }
    fprintf(lout,"\n");
    fclose(out);
    fclose(lout);
  }
  else if(init.geti("format")==2) {
    if(bflag) {
      //for binary========================================//
      FILE *bout;
      if ( (bout = fopen((init.gets("out")+"_binary.dat").c_str(),"w")) == NULL) {
        fprintf(stderr,"Error: Cannot open input file %s.\n",(init.gets("out")+"_binary.dat").c_str());
        return 0;
      }

      for(int i=0;i<nb;i++) 
        fprintf(bout,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",dat[i].e,log10(dat[i].a),dat[i].m1,dat[i].m2,dat[i].x[0],dat[i].x[1],dat[i].x[2],dat[i].v[0],dat[i].v[1],dat[i].v[2]);

      fclose(bout);
    }

    //for single========================================//    
    FILE *sout;
    if ( (sout = fopen((init.gets("out")+"_single.dat").c_str(),"w")) == NULL) {
      fprintf(stderr,"Error: Cannot open input file %s.\n",(init.gets("out")+"_single.dat").c_str());
      return 0;
    }

    for(int i=nb;i<ntt;i++)
      fprintf(sout,"%lg %lg %lg %lg %lg %lg %lg\n",dat[i].m1,dat[i].x[0],dat[i].x[1],dat[i].x[2],dat[i].v[0],dat[i].v[1],dat[i].v[2]);

    fclose(sout);
  }
  
  fprintf(stdout,"Finished. Mass single: %lg Mass binary: %lg Total objects: %d KS binaries: %d Wide binaries %d\n",mstot,mbtot,ntt,bcount,bwcount);
  return 0;
}
