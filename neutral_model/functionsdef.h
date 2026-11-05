//#define GAUSSIAN
/* functions in rnd.c */
extern void ranf_array(double *,int);
extern void ranf_start(long ); 
extern double ranf_arr_cycle(void);
extern void myrndstart(int );
extern double myrnd(void); 
extern void gauss2(double *,double );
/* functions in debug.c */
void printtable(void);
void checktable(void);
/* functions in moving.c */
#ifdef GAUSSIAN
    void moveg(int pick,int jx,int jy);
#else
extern void movedef(int *);
extern void move(int , int );
extern void moveundef(void);
#endif

/* functions in statistics.c */
extern int generascale(double);
extern void computeSAR(int **, double **,int, int ,int);
extern void computeABUN(int **,int);
extern double prestonplot(double *,int);
extern void  printconf(int,int **);


/* functions in lookuptable.c */
extern int posinside(int,int);
extern void put_in_table(int);
extern void remove_from_table(int);

/* functions in evolution.c */
extern void initsystem(void);
extern int dynamics(int);
extern void placespecies(int **);
extern void swap(int,int);
extern int Kill(int);
extern int search(int);
extern int Coalesce (int, int);

