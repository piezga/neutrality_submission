#include <time.h>
#include "rndknuth.h"


int myrndstart(int seme){ 
#ifndef FIXEDSEED
  seme=time(NULL);
#endif

#ifdef KNUTH  
  fprintf(stderr, "##start PRNG KNUTH with seed=%d\n",seme);
  ranf_start(seme);
#else
  fprintf(stderr, "##start PRNG srand48 with seed=%d\n",seme);
  srand48(seme);
#endif
  return seme;
}



double myrnd(){ 
#ifdef KNUTH
return ranf_arr_next();
#else
return drand48();
#endif
}


double exprnd(double a){
  return -log(myrnd())/a;
}



double gauss(void){
  double x,y,r,fac=0;
  while(fac==0){
    x=2.0*myrnd()-1.0;
    y=2.0*myrnd()-1.0;
    r=x*x+y*y;
    if(r<1.0 && r!=0.0){
      fac = sqrt(-2.0*log(r)/r);
      return x*fac;
    }
  }
}



void gauss2(double G[2],double sigma){
  double x,y,r,fac=0;
  while(fac==0){
    x=2.0*myrnd()-1.0;
    y=2.0*myrnd()-1.0;
    r=x*x+y*y;
    if(r<1.0 && r!=0.0){
      fac = sqrt(-2.0*log(r)/r);
      G[0]= x*fac*sigma;
      G[1]= y*fac*sigma;
    }
  }
}

