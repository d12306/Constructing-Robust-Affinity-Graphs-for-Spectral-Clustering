#ifndef RFUTILS_H
#define RFUTILS_H

#define MAX_UINT_COKUS 4294967295  //basically 2^32-1

typedef unsigned long uint32;
extern void seedMT(uint32 seed);
extern uint32 reloadMT(void);
extern uint32 randomMT(void);
extern double unif_rand();

void zeroInt(int *x, int length) ;
void zeroDouble(double *x, int length);

void createClass(double *x, int realN, int totalN, int mdim, int num_nv_src_discrete, int* Y_nv_src_discrete, int num_nv_src_continuous, double* Y_nv_src_continuous);
void normClassWt(int *cl, const int nsample, const int nclass, 
                 const int useWt, double *classwt, int *classFreq);
void makeA(double *x, const int mdim, const int nsample, int *cat, int *a, int *b);
void modA(int *a, int *nuse, const int nsample, const int mdim, int *cat, const int maxcat, int *ncase, int *jin);
void Xtranslate(double *x, int mdim, int nrnodes, int nsample, 
		int *bestvar, int *bestsplit, int *bestsplitnext,
		double *xbestsplit, int *nodestatus, int *cat, int treeSize) ;
void permuteOOB(int m, double *x, int *in, int nsample, int mdim) ;
void computeProximity(double *prox, int oobprox, int *node, int *inbag, 
                      int *oobpair, int n);
int pack(int nBits, int *bits);
void unpack(unsigned int pack, int *bits) ;



#endif