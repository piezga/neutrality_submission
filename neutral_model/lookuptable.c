#include <stdio.h>
#include <stdlib.h>
#ifndef MACOSX
#include <malloc.h>
#endif
#include "common.h"
#include "functionsdef.h"
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

****************  FUNCTIONS FOR TABLE MANIPULATION ****************

<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/* ************** LOCATE WALKER POSITION IN TABLE *********** */

 int posinside(int ix, int iy){
  int kx,ky;  

  if(ix>=0) kx=ix-(ix/NX)*NX;
  else kx= ix+(( (-ix-1)/NX) +1)*NX;

  if(iy>=0) ky=iy-(iy/NY)*NY;
  else ky= iy+(( (-iy-1)/NY) +1)*NY;

  return kx*NY+ky;
}


/* ****************  PUT WALKER POSITION IN TABLE ************** */

void put_in_table(int toput){
int pt;
#ifdef VERBOSE
 int k;
#endif
 pt=posinside((walker+toput)->ix,(walker+toput)->iy);
 if((table+pt)->numdata+1<=(table+pt)->size)
   (table+pt)->list[((table+pt)->numdata)++]=toput;
 else{
#ifdef VERBOSE
   fprintf(stderr,"reallocation for table %d\n", pt);
#endif
   (table+pt)->size *=2;
   (table+pt)->list =(int *) realloc((table+pt)->list, sizeof(int)*(table+pt)->size);
   (table+pt)->list[((table+pt)->numdata)++]=toput;

#ifdef VERBOSE
   fprintf(stderr,"table[%d] size=%d numdata=%d\n with elements \n", 
	   pt,(table+pt)->size,(table+pt)->numdata);
   for(k=0;k<(table+pt)->numdata;k++){
     fprintf(stderr,"[%d]=%d->%d (posx=%d,posy=%d)\t",
	     k,(table+pt)->list[k],(walker+(table+pt)->list[k])->name,
	     (walker+(table+pt)->list[k])->ix,(walker+(table+pt)->list[k])->iy);
   }
   fprintf(stderr,"\n");
#endif

 }  
}

/* *************** REMOVE WALKER POSITION FROM TABLE *************** */

void remove_from_table(int toremove){
  int pt,k;
  pt=posinside((walker+toremove)->ix,(walker+toremove)->iy);  

  k=0;
  while((table+pt)->list[k]!=toremove) k++;
  
  (table+pt)->list[k]=(table+pt)->list[--(table+pt)->numdata];
}

