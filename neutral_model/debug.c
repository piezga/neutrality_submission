#include <stdlib.h>
#include <stdio.h>
#ifndef MACOSX
#include <malloc.h>
#endif
#include "common.h"
#include "functionsdef.h"


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

******************** TESTING FUNCTIONS **********************

<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/* *************** PRINT INFO CONTAINED IN TABLE ****************** */
void printtable(void){
int k,j;
 for(k=0;k<NTOT;k++){
   fprintf(stderr,"#table[%d]  numdata=%d ",k,(table+k)->numdata);
   fprintf(stderr,"size=%d\t",(table+k)->size);
   for(j=0;j<(table+k)->numdata;j++)
     fprintf(stderr,"[%d]=%d->%d\t",
	     j,(table+k)->list[j],(walker+(table+k)->list[j])->name);
   fprintf(stderr,"\n");
 }
}

/* *************** CHECK CONSISTENCY OF TABLE ****************** */

void checktable(void){
int k;
int sumdata,maxdata,maxsize;

sumdata=0;
maxdata=0;
maxsize=0;
for(k=0;k<NTOT;k++){
  sumdata +=(table+k)->numdata;
  if((table+k)->numdata>maxdata)maxdata=(table+k)->numdata;
  if((table+k)->size>maxsize)maxsize=(table+k)->size;  
 }
 if(sumdata!=NW) fprintf(stderr,"NW=%d, sumdata=%d maxdata=%d, maxsize=%d\n", NW,sumdata,maxdata,maxsize);
}







