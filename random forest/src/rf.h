/**************************************************************
 * mex interface to Andy Liaw et al.'s C code (used in R package randomForest)
 * Added by Abhishek Jaiantilal ( abhishek.jaiantilal@colorado.edu )
 * License: GPLv2
 * Version: 0.02
 *
 * other than adding the macros for F77_* and adding this message
 * nothing changed .
 *************************************************************/

/*******************************************************************
   Copyright (C) 2001-7 Leo Breiman, Adele Cutler and Merck & Co., Inc.
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.                            
*******************************************************************/
#ifndef RF_H
#define RF_H

/* test if the bit at position pos is turned on */
#define isBitOn(x,pos) (((x) & (1 << (pos))) > 0)
/* swap two integers */
#define swapInt(a, b) ((a ^= b), (b ^= a), (a ^= b))

#define F77_CALL(x) x ## _
#define F77_NAME(x) F77_CALL(x)
#define F77_SUB(x) F77_CALL(x)


void normClassWt(int *cl, const int nsample, const int nclass, 
                 const int useWt, double *classwt, int *classFreq);

void classForest(int *mdim, int *ntest, int *nclass, int *maxcat, 
				int *nrnodes, int *jbt, double *xts, double *xbestsplit, 
				double *pid, double *cutoff, double *countts, int *treemap, 
				int *nodestatus, int *cat, int *nodeclass, int *jts, 
				int *jet, int *bestvar, int *nodexts, int *ndbigtree, 
				int *keepPred, int *prox, double *proxmatrix, int *nodes,
				int *path_stat, int path_length);

void regTree(double *x, double *y, int mdim, int nsample, 
	     int *lDaughter, int *rDaughter, double *upper, double *avnode, 
             int *nodestatus, int nrnodes, int *treeSize, int nthsize, 
             int mtry, int *mbest, int *cat, double *tgini, int *varUsed);

void findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample, 
		   int ndstart, int ndend, int *msplit, double *decsplit, 
		   double *ubest, int *ndendl, int *jstat, int mtry,
		   double sumnode, int nodecnt, int *cat);

void predictRegTree(double *x, int nsample, int mdim, 
		    int *lDaughter, int *rDaughter, int *nodestatus, 
                    double *ypred, double *split, double *nodepred, 
                    int *splitVar, int treeSize, int *cat, int maxcat,
                    int *nodex);

void predictClassTree(double *x, int n, int mdim, int *treemap,
		      int *nodestatus, double *xbestsplit,
		      int *bestvar, int *nodeclass,
		      int ndbigtree, int *cat, int nclass,
		      int *jts, int *nodex, int *path_stat, int maxcat);

			  
int pack(int l, int *icat);
void unpack(unsigned int npack, int *icat);

void zeroInt(int *x, int length);
void zeroDouble(double *x, int length);
void createClass(double *x, int realN, int totalN, int mdim, int num_nv_src_discrete, int* Y_nv_src_discrete, int num_nv_src_continuous, double* Y_nv_src_continuous);
void prepare(int *cl, const int nsample, const int nclass, const int ipi, 
	     double *pi, double *pid, int *nc, double *wtt);
//void makeA(double *x, const int mdim, const int nsample, int *cat, int *a, int *b);
void modA(int *a, int *nuse, const int nsample, const int mdim, int *cat, 
          const int maxcat, int *ncase, int *jin);
void Xtranslate(double *x, int mdim, int nrnodes, int nsample, 
		int *bestvar, int *bestsplit, int *bestsplitnext,
		double *xbestsplit, int *nodestatus, int *cat, int treeSize);
void permuteOOB(int m, double *x, int *in, int nsample, int mdim);
void computeProximity(double *prox, int oobprox, int *node, int *inbag, 
                      int *oobpair, int n);

/* Node status */
#define NODE_TERMINAL -1
#define NODE_TOSPLIT  -2
#define NODE_INTERIOR -3

#endif /* RF_H */
