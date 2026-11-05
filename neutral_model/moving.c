#include <stdlib.h>
#include <stdio.h>
#ifndef MACOSX
#include <malloc.h>
#endif
#include "common.h"
#include "functionsdef.h"



/* ***************** INITIALIZE KERNEL FOR MOVEMENTS **************** */
#ifdef GAUSSIAN
void moveg(int i, int kx,int ky){
  (walker+i)->ix += kx;     
  (walker+i)->iy += ky;     
}

#else
void movedef(int *nummosse){
  int ix,iy,nm;

  if(kernsize==-1){
    *nummosse=4;
    mx = (int *) malloc(sizeof(int)*(*nummosse));
    my = (int *) malloc(sizeof(int)*(*nummosse));
    nm=0;
    for(ix=-1;ix<=1;ix++){
      for(iy=-1;iy<=1;iy++){
	if(ix!=iy && (ix+iy)!=0){
	  mx[nm]=ix;
	  my[nm]=iy;
	  nm++;
	}
      }
    }
  }else{
    *nummosse=(2*kernsize+1)*(2*kernsize+1)-1;
    mx = (int *) malloc(sizeof(int)*(*nummosse));
    my = (int *) malloc(sizeof(int)*(*nummosse));
    nm=0;
    for(ix=-kernsize;ix<=kernsize;ix++){
      for(iy=-kernsize;iy<=kernsize;iy++){
	if(ix!=0 || iy!=0) {
	  mx[nm]=ix;
	  my[nm]=iy;
	  //	  fprintf(stderr,"%d (%d,%d)\n",nm,ix,iy);
	  nm++;
	}
      }
    }
  }
    if(nm!=(*nummosse)){
      fprintf(stderr,"ERROR in movedef nm=%d nummosse=%d\n",nm,*nummosse);
    exit(9);
  }
}


/* ***************** RULE  FOR MOVEMENTS **************** */
void move(int i, int mossa){
  (walker+i)->ix += mx[mossa];     
  (walker+i)->iy += my[mossa];     
}



void moveundef(void){
  free(mx);
  free(my);
}
#endif
