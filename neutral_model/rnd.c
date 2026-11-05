#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functionsdef.h"

#define KK 100                 
#define LL  37                 
#define mod_sum(x,y) (((x)+(y))-(int)((x)+(y))) 
#define TT  70
#define is_odd(s) ((s)&1)

double ran_u[KK];  

void ranf_array(double *aa,int n)
{
  register int i,j;
  for (j=0;j<KK;j++) aa[j]=ran_u[j];
  for (;j<n;j++) aa[j]=mod_sum(aa[j-KK],aa[j-LL]);
  for (i=0;i<LL;i++,j++) ran_u[i]=mod_sum(aa[j-KK],aa[j-LL]);
  for (;i<KK;i++,j++) ran_u[i]=mod_sum(aa[j-KK],ran_u[i-LL]);
}


void ranf_start(long seed) 
{
  register int t,s,j;
  double u[KK+KK-1],ul[KK+KK-1];
  double ulp=(1.0/(1L<<30))/(1L<<22); 
  double ss=2.0*ulp*((seed&0x3fffffff)+2);

  for (j=0;j<KK;j++) {
    u[j]=ss; ul[j]=0.0;               
    ss+=ss; if (ss>=1.0) ss-=1.0-2*ulp;
  }
  for (;j<KK+KK-1;j++) u[j]=ul[j]=0.0;
  u[1]+=ulp;ul[1]=ulp;  
  s=seed&0x3fffffff;
  t=TT-1; while (t) {
    for (j=KK-1;j>0;j--) ul[j+j]=ul[j],u[j+j]=u[j];  
    for (j=KK+KK-2;j>KK-LL;j-=2)
        ul[KK+KK-1-j]=0.0,u[KK+KK-1-j]=u[j]-ul[j];
    for (j=KK+KK-2;j>=KK;j--) if(ul[j]) {
      ul[j-(KK-LL)]=ulp-ul[j-(KK-LL)],
        u[j-(KK-LL)]=mod_sum(u[j-(KK-LL)],u[j]);
      ul[j-KK]=ulp-ul[j-KK],u[j-KK]=mod_sum(u[j-KK],u[j]);
    }
    if (is_odd(s)) {                           
      for (j=KK;j>0;j--)  ul[j]=ul[j-1],u[j]=u[j-1];
      ul[0]=ul[KK],u[0]=u[KK];      
      if (ul[KK]) ul[LL]=ulp-ul[LL],u[LL]=mod_sum(u[LL],u[KK]);
    }
    if (s) s>>=1; else t--;
  }
  for (j=0;j<LL;j++) ran_u[j+KK-LL]=u[j];
  for (;j<KK;j++) ran_u[j-LL]=u[j];
}



#define QUALITY 1009
double ranf_arr_buf[QUALITY];
double ranf_arr_sentinel=-1.0;
double *ranf_arr_ptr=&ranf_arr_sentinel;

#define ranf_arr_next() (*ranf_arr_ptr>=0? *ranf_arr_ptr++: ranf_arr_cycle())
double ranf_arr_cycle()
{
  ranf_array(ranf_arr_buf,QUALITY);
  ranf_arr_buf[100]=-1;
  ranf_arr_ptr=ranf_arr_buf+1;
  return ranf_arr_buf[0];
}





void myrndstart(int seme){ 
#ifdef DEBUG
  ranf_start(seme);
#else
  /* V2D_SEED (environment) fixes the seed, for reproducible runs; otherwise the clock */
  char *env = getenv("V2D_SEED");
  if (env) ranf_start(atol(env) + seme);
  else ranf_start(time(NULL));
#endif
}



double myrnd(void){ 
  return ranf_arr_next();
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





