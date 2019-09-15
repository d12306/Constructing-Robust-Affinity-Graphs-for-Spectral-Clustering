#ifndef CLASSTREE_H
#define CLASSTREE_H
//----------------------
#ifdef MATLAB
#define Rprintf  mexPrintf
#include "mex.h"
#endif

/*#ifndef MATLAB
#define Rprintf printf
#include "mex_stub.h"
#include "stdio.h"
#endif*/

typedef unsigned long uint32;
extern void seedMT(uint32 seed);
extern uint32 reloadMT(void);
extern uint32 randomMT(void);
extern double unif_rand();
extern void R_qsort_I(double *v, int *I, int i, int j);

extern "C"{
    #ifdef WIN64
	void _catmax_(double *parentDen, double *tclasscat,
                      double *tclasspop, int *nclass, int *lcat,
                      int *ncatsp, double *critmax, int *nhit,
                      int *maxcat, int *ncmax, int *ncsplit);
    #endif
    
    #ifndef WIN64
	void catmax_(double *parentDen, double *tclasscat,
                      double *tclasspop, int *nclass, int *lcat,
                      int *ncatsp, double *critmax, int *nhit,
                      int *maxcat, int *ncmax, int *ncsplit);
    #endif                
}
extern "C"{
    #ifdef WIN64
	void _catmaxb_(double *totalWt, double *tclasscat, double *classCount,
                       int *nclass, int *nCat, int *nbest, double *critmax,
                       int *nhit, double *catCount) ;
    #endif
            
    #ifndef WIN64
	void catmaxb_(double *totalWt, double *tclasscat, double *classCount,
                       int *nclass, int *nCat, int *nbest, double *critmax,
                       int *nhit, double *catCount) ;
    #endif
}
//----------------------

#ifdef WIN64
void F77_NAME(_catmax)
#endif

#ifndef WIN64
void F77_NAME(catmax)
#endif
(double *parentDen, double *tclasscat,
                      double *tclasspop, int *nclass, int *lcat,
                      int *ncatsp, double *critmax, int *nhit,
                      int *maxcat, int *ncmax, int *ncsplit) ;
#ifdef WIN64
void F77_NAME(_catmaxb)
#endif
#ifndef WIN64
void F77_NAME(catmaxb)
#endif
(double *totalWt, double *tclasscat, double *classCount,
                       int *nclass, int *nCat, int *nbest, double *critmax,
                       int *nhit, double *catCount) ;
void predictClassTree(double *x, int n, int mdim, int *treemap,
		      int *nodestatus, double *xbestsplit,
		      int *bestvar, int *nodeclass,
		      int treeSize, int *cat, int nclass,
		      int *jts, int *nodex, int maxcat) ;


#endif