/* 
Voter Model in due dimensioni 

versione NX x NY (per fare domini non quadrati)

*/   

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"
#include "gvar.h"
#include "functionsdef.h"


/* >>>>>>>>>> M A I N <<<<<<<<<<<<<< */



int main(int argc, char *argv[]){
int NREL,irel,iq;
int j,k,kx,ky;
int nummosse;
int **speciesmap;
int numscale;
double sepscale;
double **sar;
double *NS_vs_NI,mediaspecie;
char nomefi[128];
FILE *fout;

 if(argc!=6){
   fprintf(stderr, "usage %s <NX> <NY> <NU> <NUMREL> <kernel>\n",argv[0]);
   exit(9);
 }
 sscanf(argv[1],"%d",&NX);
 sscanf(argv[2],"%d",&NY);
 sscanf(argv[3],"%lf",&NU);
 sscanf(argv[4],"%d",&NREL);
 sscanf(argv[5],"%d",&kernsize);
 NTOT=NX*NY;

 

 fprintf(stderr,"## Voter model \n");
 fprintf(stderr,"## Size_X=%d Size_Y=%d \n",NX,NY);
 fprintf(stderr,"## Speciaction Rate=%e \n",NU);




 /* initialize some stuff */
 myrndstart(2111);
#ifdef GAUSSIAN
 if(kernsize<=1){
   fprintf(stderr, "Kernsize should be positive for Gaussian\n");
   exit(1);
 }
 fprintf(stderr,"#Gaussian Kernel with sigma=%d\n",kernsize);
#else
 fprintf(stderr, "I'm using square Kernel\n");
 fprintf(stderr, "if kernsize=-1 use NN kernel otherwise square with halfsize=kernsize\n");
 movedef(&nummosse);
#endif
 //fprintf(stderr, "numero di mosse %d\n",nummosse);

 sepscale=sqrt(sqrt(2.0));
 numscale=generascale(sepscale);



 sar = (double **) malloc(sizeof(double *) * NQ);
 for(iq=0;iq<NQ;iq++)
   sar[iq] = (double *) malloc(sizeof(double) * numscale);

 NS_vs_NI = (double *) malloc(sizeof(double) * NTOT);

 for(iq=0;iq<NQ;iq++)
   for(j=0;j<numscale;j++)
     sar[iq][j]=0.0;


 for(j=0;j<NTOT;j++)
   NS_vs_NI[j]=0.0;





 /* ++++++++++++ MAIN LOOP +++++++++++++++++++ */
 mediaspecie=0.0;
 for(irel=1;irel<=NREL;irel++){

   initsystem();
   dynamics(nummosse);
   mediaspecie+= (double)NS;
   fprintf(stderr, "## Realization=%d Nspecies=%d\n",irel,NS);


   if ((speciesmap = (int **) malloc(sizeof(int *) * NX) ) == NULL) {
     fprintf(stderr,"ERROR ALLOCATING speciesmap\n");
     exit(9);
   }   
   for (kx=0;kx<NX;kx++){
     if ((speciesmap[kx] = (int *) malloc(sizeof(int) * NY) ) == NULL) {
       fprintf(stderr,"ERROR ALLOCATING speciesmap\n");
       exit(9);
     } 
   }

   
   placespecies(speciesmap);

#ifdef DUMP
/* dump configuration if you not specify BINARY it will be printed in ASCII */
   printconf(irel,speciesmap);
#endif

/* PUT HERE CALLS TO FUNCTIONS FOR STATISTICAL ANALYSIS which are contained in statistics.c */
//   computeSAR(speciesmap,sar,numscale, irel,NREL);
   computeABUN(speciesmap,irel);
   

   for (kx=0;kx<NX;kx++)
     free(speciesmap[kx]);   
   free(speciesmap);
 }
 fprintf(stderr,"FINE\n");

 /* RAW data per fare il Preston plot (mediati su tutte le realizzazioni */
 /*
 sprintf(nomefi,"NS_vs_NI_NX%d_NY%d_NU%e_NREL%d.ker%d",NX,NY,NU,NREL,kernsize);
 fout=fopen(nomefi,"w");
 for (k=0;k<NTOT;k++)
   if(NS_vs_NI[k]>0)fprintf(fout,"%d %g\n",k,NS_vs_NI[k]/(double) NREL);
 fclose(fout);

 fprintf(stderr, "in media il sistema contiene=%g specie\n",mediaspecie/(double) NREL);
 mediaspecie= prestonplot(NS_vs_NI,NREL);
 
 fprintf(stderr, "Somma da preston plot=%g\n",mediaspecie);

 free(NS_vs_NI);


 for(iq=0;iq<NQ;iq++)
   free(sar[j]);
 free(sar);

 */
#ifndef GAUSSIAN
  moveundef();
#endif  

 exit(1);
}






