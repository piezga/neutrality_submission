#include <stdlib.h>
#include <stdio.h>
#ifndef MACOSX
#include <malloc.h>
#endif
#include "common.h"
#include "functionsdef.h"


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

**********  FUNCTIONS FOR INITIALIZATION  *******

<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/* **********  INITIALIZE WALKERS AND TABLE  ******* */

void initsystem(void){
  int kx,ky,pos;  

  if((walker = (walker_type *) malloc(sizeof(walker_type) * NTOT) ) == NULL) {
    fprintf(stderr,"ERROR ALLOCATING walker\n");
    exit(9);
  }
 
  if((table = (table_type *) malloc(sizeof(table_type) * NTOT) ) == NULL) {
    fprintf(stderr,"ERROR ALLOCATING table\n");
    exit(9);
  }
  
  NW=0;
  for(kx=0;kx<NX;kx++){
    for(ky=0;ky<NY;ky++){
      pos=kx*NY+ky;
      //      fprintf(stderr,"kx=%d ky=%d pos=%d NW=%d\n", kx,ky,pos,NW);
      (walker+NW)->ix=kx;
      (walker+NW)->iy=ky;
      (walker+NW)->name= pos;
      (walker+NW)->nameptr= -1;
      (walker+NW)->specie= -1; 
      (table+pos)->size=2;
      (table+pos)->numdata=0 ;
      (table+pos)->list =(int *) malloc(sizeof(int)*(table+pos)->size);
      (table+pos)->list[(table+pos)->numdata++]=NW;
      NW++;
    }
  }
#ifdef VERBOSE
  fprintf(stderr,"first check\n");
  checktable();
#endif

  // useless check
  if(NW!=NTOT)
    fprintf(stderr,"Some problem nw=%d is different from NTOT=%d\n", NW,NTOT); 
}












/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

****************  FUNCTIONS FOR THE DYNAMICS **********************

<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/* ************* SWAP WALKER i AND j **************************** */

void swap(int i,int j){
  walker_type walker_tmp; 
  walker_tmp  = *(walker+i); 
  *(walker+i) = *(walker+j); 
  *(walker+j) = walker_tmp;
}

/*  *************  KILL A WALKER AND CREATE NEW SPECIES ********* */

int Kill(int pick){

 if(pick!=NW-1){
   remove_from_table(pick);
   remove_from_table(NW-1);
   swap(pick,NW-1);
   put_in_table(pick);
 }
 else
   remove_from_table(pick);
 (walker+(NW-1))->specie= NS; //assign species   
 NS++;
 NW--;

 return 0;
}

/*  *************  SEARCH PARTNER OF i FOR COALESCING ********* */

int search(int pick){
int k,j,pt;  

 pt=posinside((walker+pick)->ix,(walker+pick)->iy);
 j=-1;
 for(k=0;k<(table+pt)->numdata;k++)
   if((walker+pick)->ix==(walker+(table+pt)->list[k])->ix)j=(table+pt)->list[k];
 return j;
}

/*  *************  COALESCENCE BETWEEN WALKER i AND j ********* */

int Coalesce (int i, int j){

#ifdef VERBOSE
 fprintf(stderr,"coalescence of (%d,%d)=",i,j);
 fprintf(stderr,"(name=%d,name=%d)\n",(walker+i)->name,(walker+j)->name);
#endif


 (walker+j)->nameptr= (walker+i)->name; 

 if(j!=NW-1 && i!=NW-1){
   remove_from_table(NW-1);
   remove_from_table(j);
   swap(j,NW-1);
   put_in_table(j);
   put_in_table(i);
 } else if(j==NW-1){
   remove_from_table(j);
   put_in_table(i);
 } else if(i==NW-1){
   remove_from_table(j);
   swap(j,NW-1);
   put_in_table(j);
 }
 NW--;
 return 0;
}



/* *****  EVOLVE WALKERS TILL ONLY ONE IS SURVIVED (CORE FUNCTION FOR EVOLUTION) ***** */

int dynamics(int nummosse){
int  mossa,pick,j;
#ifdef GAUSSIAN
    int  jx,jy,imove;
    double G[2];
#endif

NS=0;  
while(NW>1){

#ifdef VERBOSE
  checktable(); //printtable();
  //  fprintf(stderr, "\n\n");
#endif

  pick=(int) (NW*myrnd()); 
  if(myrnd()<NU){ 
#ifdef VERBOSE
    fprintf(stderr, "Kill %d->%d in pos=%d\n", pick,(walker+pick)->name,posinside((walker+pick)->ix,(walker+pick)->iy));
#endif
    Kill(pick);
  }else{            
    remove_from_table(pick);          
    
#ifdef GAUSSIAN
    imove=1;
    while(imove){
      gauss2(G,kernsize);
      jx=(int) G[0];
      jy=(int) G[1];
      if(jx!=0 || jy!=0)imove=0;
    }
    moveg(pick,jx,jy);
#else    
    mossa=(int) (nummosse*myrnd());  
    move(pick,mossa); 
#endif
    j=search(pick);
    if(j!=-1)Coalesce(pick,j);
    else put_in_table(pick);     
  }
 }

 Kill(0);
#ifdef VERBOSE
 pick=0;
 fprintf(stderr, "Kill %d->%d in pos=%d\n", pick,(walker+pick)->name,posinside((walker+pick)->ix,(walker+pick)->iy));
 // printtable();
 fprintf(stderr,"exit dynamics\n");
#endif
 return NS;
}

/*  **************** RECONSTRUCT SPECIES POSITIONS ************** */

void placespecies(int **mappa){
  int k,kx,ky,j,pos,specie,next;
  species_type *habitat;




  /* deallocate lookup table */
  for(k=0;k<NTOT;k++)
    free((table+k)->list);
  free(table);

  /* allocate temporary array for species location */
  if ((habitat = (species_type *) malloc(sizeof(species_type) * NTOT) ) == NULL) {
    fprintf(stderr,"ERROR ALLOCATING habitat \n");
    exit(9);
  } 

  /* Reorganize walker in the initial configuration taking track of species and pointer to colliding partner */
  for (k=0;k<NTOT;k++){
    j=(walker+k)->name;
    (habitat+j)->specie=(walker+k)->specie;
    (habitat+j)->nameptr=(walker+k)->nameptr;
  }
 
  /* deallocate walkers */
  free(walker);   



  /* define occupation by individuals of the various species */
  for(kx=0;kx<NX;kx++){  
    for(ky=0;ky<NY;ky++){  
      pos=kx*NY+ky;
      specie=(habitat+pos)->specie;
      if(specie==-1){
	next=(habitat+pos)->nameptr;
	while(specie==-1){
	  specie=(habitat+next)->specie;
	  next=(habitat+next)->nameptr;
	}
	(habitat+pos)->specie=specie;
      }
      mappa[kx][ky]=specie;
    }
  }
  
  /* deallocate temporary array for species location */
  free(habitat);

}


