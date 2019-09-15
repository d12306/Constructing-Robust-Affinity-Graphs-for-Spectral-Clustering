/**************************************************************
 * mex interface to Andy Liaw et al.'s C code (used in R package randomForest)
 * Added by Abhishek Jaiantilal ( abhishek.jaiantilal@colorado.edu )
 * License: GPLv2
 * Version: 0.02
 *
 * File: contains all the supporting code for a standalone C or mex for
 *       Classification RF. 
 * Copied all the code from the randomForest 4.5-28 or was it -29?
 * 
 * important changes (other than the many commented out printf's)
 * 1. realized that instead of changing individual S_allocs to callocs
 *    a better way is to emulate them
 * 2. found some places where memory is not freed in classRF via valgrind so 
 *    added frees
 * 3. made sure that C can now interface with brieman's fortran code so added 
 *    externs "C"'s and the F77_* macros
 * 4. added cokus's mersenne twister.
 *
 *************************************************************/

/*****************************************************************
 * Copyright (C) 2001-7 Leo Breiman, Adele Cutler and Merck & Co., Inc.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * C driver for Breiman & Cutler's random forest code.
 * Re-written from the original main program in Fortran.
 * Andy Liaw Feb. 7, 2002.
 * Modifications to get the forest out Matt Wiener Feb. 26, 2002.
 *****************************************************************/

#include "stdlib.h"
#include "memory.h"
#include "stdio.h"
#include "math.h"
#include "time.h"
#include <iostream>
#include <fstream>
#include "rf.h"
#include "rfutils.h"
#include "buildtree.h"
#include "classRF.h"

#ifndef MATLAB
#include "mex_stub.h"
#endif

#ifdef MATLAB
#include "mex.h"
#define Rprintf mexPrintf
#endif

#define F77_CALL(x) x ## _
#define F77_NAME(x) F77_CALL(x)
#define F77_SUB(x) F77_CALL(x)


#define MAX_UINT_COKUS 4294967295  //basically 2^32-1

typedef unsigned long uint32;
extern void seedMT(uint32 seed);
extern uint32 reloadMT(void);
extern uint32 randomMT(void);

double unif_rand(){
    return (((double)randomMT())/((double)MAX_UINT_COKUS));
}

void GetRNGstate(){};
void PutRNGstate(){};

void oob(int nsample, int nclass, int *jin, int *cl, int *jtr, int *jerr,
        int *counttr, int *out, double *errtr, int *jest, double *cutoff);

void TestSetError(double *countts, int *jts, int *clts, int *jet, int ntest,
        int nclass, int nvote, double *errts,
        int labelts, int *nclts, double *cutoff);


/* Define the R RNG for use from Fortran. */
#ifdef WIN64
void _rrand_(double *r) { *r = unif_rand(); }
#endif

#ifndef WIN64
void rrand_(double *r) { *r = unif_rand(); }
#endif

/*
	Varialbes:
	- x: data points, a matrix with $mdim$(row) x $nsample$(column)
	
	The papameters w.r.t. non-visual data sources:
	- num_nv_src_discrete: the number of discrete non-visual sources;
	- Y_nv_src_discrete: the response of discrete non-visual data, a matrix with $nsample$(row) x $num_nv_src_discrete$(column)
	- alpha_nv_src_discrete: the weight of discrete non-visual sources;
	- max_num_class_nv_src_discrete: the maximam class number of discrete non-visual data;
	- num_nv_src_continuous: the number of contunious non-visual sources;
	- Y_nv_src_continuous: the response of continuous non-visual data, a matrix with $nsample$(row) x $num_nv_src_continuous$(column);
	- alpha_nv_src_continuous: the weight of continuous non-visual sources;
*/

void classRF(double *x, int *dimx, int *cl, int *ncl, int *cat, int *maxcat,
			int *sampsize, int *strata, int *Options, int *ntree, int *nvar,
			int *ipi, double *classwt, double *cut, int *nodesize,
			int *outcl, int *counttr, double *prox,
			double *imprt, double *impsd, double *impmat, int *nrnodes,
			int *ndbigtree, int *nodestatus, int *bestvar, int *treemap,
			int *nodeclass, double *xbestsplit, double *errtr,
			int *testdat, double *xts, int *clts, int *nts, double *countts,
			int *outclts, int labelts, double *proxts, double *errts, 
			int *inbag, int print_verbose_tree_progression,
			int num_nv_src_discrete, int *Y_nv_src_discrete, 
			double *alpha_nv_src_discrete, int max_num_class_nv_src_discrete,
			int num_nv_src_continuous, double *Y_nv_src_continuous, 
			double *alpha_nv_src_continuous, double *meanSurrVarAssoc,
			double *visual_info_gain, double *num_flat_region, double *num_split,
			int *node_size_table, double *pseudo_X) 
{
    /******************************************************************
     *  C wrapper for random forests:  get input from R and drive
     *  the Fortran routines.
     *
     *  Input:
     *
     *  x:        matrix of predictors (transposed!)
     *  dimx:     two integers: number of variables and number of cases
     *  cl:       class labels of the data
     *  ncl:      number of classes in the responsema
     *  cat:      integer vector of number of classes in the predictor;
     *            1=continuous
     * maxcat:    maximum of cat
     * Options:   7 integers: (0=no, 1=yes)
     *     add a second class (for unsupervised RF)?
     *         1: sampling from product of marginals
     *         2: sampling from product of uniforms
     *     assess variable importance?
     *     calculate proximity?
     *     calculate proximity based on OOB predictions?
     *     calculate outlying measure?
     *     how often to print output?
     *     keep the forest for future prediction?
     *  ntree:    number of trees
     *  nvar:     number of predictors to use for each split
     *  ipi:      0=use class proportion as prob.; 1=use supplied priors
     *  pi:       double vector of class priors
     *  nodesize: minimum node size: no node with fewer than ndsize
     *            cases will be split
     *
     *  Output:
     *
     *  outcl:    class predicted by RF
     *  counttr:  matrix of votes (transposed!)
     *  imprt:    matrix of variable importance measures
     *  impmat:   matrix of local variable importance measures
     *  prox:     matrix of proximity (if iprox=1)
     ******************************************************************/ 
    
	
    int nsample0, mdim, nclass, addClass, mtry, ntest, nsample, ndsize,
            mimp, nimp, near, nuse, noutall, nrightall, nrightimpall,
            keepInbag, nstrata;
    int j, n, m, k, idxByNnode, idxByNsample, imp, localImp, iprox,
            oobprox, keepf, replace, stratify, trace, *nright,
            *nrightimp, *nout, *nclts = 0, Ntree;
    
    int *out, *bestsplitnext, *bestsplit, *nodepop, *jin, *nodex, *path_stat,
		*nodexts, *path_stat_xts, *nodestart, *ta, *ncase, *jerr, *varUsed,
		*jtr, *classFreq, *idmove, *jvr,
		*at, *a, *b, *mind, *nind = 0, *jts, *oobpair = 0;
    int **strata_idx = 0, *strata_size = 0, last, ktmp, anyEmpty, ntry;
    
    double av = 0.0;
    
    double *tgini, *tx, *wl, *classpop, *tclasscat, *tclasspop, *win, *tp, *wr;
	
	   
    srand(time(NULL));
    //Do initialization for COKUS's Random generator
    seedMT(2*rand()+1);  //works well with odd number so why don't use that
    
    addClass = Options[0];
    imp      = Options[1];
    localImp = Options[2];
    iprox    = Options[3];
    oobprox  = Options[4];
    trace    = Options[5];
    keepf    = Options[6];
    replace  = Options[7];
    stratify = Options[8];
    keepInbag = Options[9];
	int keep_pseudo_X = Options[10];
    mdim     = dimx[0];
    nsample0 = dimx[1];
    nclass   = (*ncl==1) ? 2 : *ncl;
    ndsize   = *nodesize;
    Ntree    = *ntree;
    mtry     = *nvar;
    ntest    = *nts;
    nsample = addClass ? (nsample0 + nsample0) : nsample0;
    mimp = imp ? mdim : 1;
    nimp = imp ? nsample : 1;
    near = iprox ? nsample0 : 1;
    if (trace == 0) trace = Ntree + 1;
	
	/* The variables for all non-visual sources*/
	// -- the weight of samples, the more times one sample is selected, the larger weight it has
	double *win_nv_src = new double[nsample];
	
	
    /* The variables w.r.t. discrete non-visual data sources */
	// -- The variables for computing the information criteria of a split
	// -- The information critierion on the splitting node, used in each splitting process
	double* parent_crit_nv_src_discrete = new double[num_nv_src_discrete];
	// ---- The numerator of the information criteria formula for the splitting node
	double* pno_nv_src_discrete = new double[num_nv_src_discrete];
	// ---- The denominator of the information criteria formula for the splitting node
	double* pdo_nv_src_discrete = new double[num_nv_src_discrete];
	// ---- The numerator of the information criteria formula for the left child
	double* rln_nv_src_discrete = new double[num_nv_src_discrete];
	// ---- The denominator of the information criteria formula for the left child
	double* rld_nv_src_discrete = new double[num_nv_src_discrete];
	// ---- The numerator of the information criteria formula for the right child
	double* rrn_nv_src_discrete = new double[num_nv_src_discrete];
	// ---- The denominator of the information criteria formula for the right child
	double* rrd_nv_src_discrete = new double[num_nv_src_discrete];
	
	
	// -- The aggregated weight of discrete non-visual data for the left child
	double* wl_nv_src_discrete = new double[max_num_class_nv_src_discrete * num_nv_src_discrete];
	// -- The aggregated weight of discrete non-visual data for the right child
	double* wr_nv_src_discrete = new double[max_num_class_nv_src_discrete * num_nv_src_discrete];
	
	// -- The category histogram of discrete non-visual data, a matrix used in each splitting process
	double* tclasspop_nv_src_discrete = new double[max_num_class_nv_src_discrete * num_nv_src_discrete];
	// -- The category histogram of discrete non-visual data for all nodes of a tree
	double* classpop_nv_src_discrete_all_nodes = new double[num_nv_src_discrete * max_num_class_nv_src_discrete * (*nrnodes)];
	
	
	/* The variables w.r.t. continuous non-visual data sources */
	// -- The information critierion on the splitting node, used in each splitting process
	double* parent_crit_nv_src_continuous = new double[num_nv_src_continuous];
	// -- The sum of non-visual response for all nodes in a tree
	double* sum_nv_src_continuous_all_nodes = new double[num_nv_src_continuous * (*nrnodes)];
	// -- The sum of non-visual response on the splitting node
	double* sum_nv_src_continuous = new double[num_nv_src_continuous];
	// -- The sum of non-visual response on the left child
	double* suml_nv_src_continuous = new double[num_nv_src_continuous];
	// -- The sum of non-visual response on the right child
	double* sumr_nv_src_continuous = new double[num_nv_src_continuous];
	
	
	
	//-----------------------------------------------
    /*printf("\nmdim %d, nclass %d, nrnodes %d, nsample %d, ntest %d\n", mdim, nclass, *nrnodes, nsample, ntest);
    printf("\noobprox %d, mdim %d, nsample0 %d, Ntree %d, mtry %d, mimp %d", oobprox, mdim, nsample0, Ntree, mtry, mimp);
    printf("\nstratify %d, replace %d",stratify,replace);
    printf("\n");*/
    tgini =      (double *) mxCalloc(mdim, sizeof(double));
    wl =         (double *) mxCalloc(nclass, sizeof(double));
    wr =         (double *) mxCalloc(nclass, sizeof(double));
    classpop =   (double *) mxCalloc(nclass* *nrnodes, sizeof(double));
    tclasscat =  (double *) mxCalloc(nclass*32, sizeof(double));
    tclasspop =  (double *) mxCalloc(nclass, sizeof(double));
    tx =         (double *) mxCalloc(nsample, sizeof(double));
    win =        (double *) mxCalloc(nsample, sizeof(double));
    tp =         (double *) mxCalloc(nsample, sizeof(double));
    out =           (int *) mxCalloc(nsample, sizeof(int));
    bestsplitnext = (int *) mxCalloc(*nrnodes, sizeof(int));
    bestsplit =     (int *) mxCalloc(*nrnodes, sizeof(int));
    nodepop =       (int *) mxCalloc(*nrnodes, sizeof(int));
    nodestart =     (int *) mxCalloc(*nrnodes, sizeof(int));
    jin =           (int *) mxCalloc(nsample, sizeof(int));
    nodex =         (int *) mxCalloc(nsample, sizeof(int));
	path_stat =		(int *) mxCalloc(nsample, sizeof(int));
    nodexts =       (int *) mxCalloc(ntest, sizeof(int));
	path_stat_xts = (int *) mxCalloc(ntest, sizeof(int));
    ta =            (int *) mxCalloc(nsample, sizeof(int));
    ncase =         (int *) mxCalloc(nsample, sizeof(int));
    jerr =          (int *) mxCalloc(nsample, sizeof(int));
    varUsed =       (int *) mxCalloc(mdim, sizeof(int));
    jtr =           (int *) mxCalloc(nsample, sizeof(int));
    jvr =           (int *) mxCalloc(nsample, sizeof(int));
    classFreq =     (int *) mxCalloc(nclass, sizeof(int));
    jts =           (int *) mxCalloc(ntest, sizeof(int));
    idmove =        (int *) mxCalloc(nsample, sizeof(int));
    at =            (int *) mxCalloc(mdim*nsample, sizeof(int));
    a =             (int *) mxCalloc(mdim*nsample, sizeof(int));
    b =             (int *) mxCalloc(mdim*nsample, sizeof(int));
    mind =          (int *) mxCalloc(mdim, sizeof(int));
    nright =        (int *) mxCalloc(nclass, sizeof(int));
    nrightimp =     (int *) mxCalloc(nclass, sizeof(int));
    nout =          (int *) mxCalloc(nclass, sizeof(int));
    if (oobprox) 
	{
        oobpair = (int *) mxCalloc(near*near, sizeof(int));
    }
	
    /* Count number of cases in each class. */
    zeroInt(classFreq, nclass);
    for (n = 0; n < nsample; ++n) classFreq[cl[n] - 1] ++;
    /* Normalize class weights. */
    normClassWt(cl, nsample, nclass, *ipi, classwt, classFreq);
   
    if (stratify) 
	{
        /* Count number of strata and frequency of each stratum. */
        nstrata = 0;
        for (n = 0; n < nsample0; ++n)
            if (strata[n] > nstrata) nstrata = strata[n];
        /* Create the array of pointers, each pointing to a vector
         * of indices of where data of each stratum is. */
        strata_size = (int *) mxCalloc(nstrata, sizeof(int));
        for (n = 0; n < nsample0; ++n) 
		{
            strata_size[strata[n] - 1] ++;
        }
        strata_idx = (int **) mxCalloc(nstrata, sizeof(int *));
        for (n = 0; n < nstrata; ++n) {
            strata_idx[n] = (int *) mxCalloc(strata_size[n], sizeof(int));
        }
        zeroInt(strata_size, nstrata);
        for (n = 0; n < nsample0; ++n) {
            strata_size[strata[n] - 1] ++;
            strata_idx[strata[n] - 1][strata_size[strata[n] - 1] - 1] = n;
        }
    } else {
        nind = replace ? NULL : (int *) mxCalloc(nsample, sizeof(int));
    }
    
    /* INITIALIZE FOR RUN */
    if (*testdat) zeroDouble(countts, ntest * nclass);
    zeroInt(counttr, nclass * nsample);
    zeroInt(out, nsample);
    zeroDouble(tgini, mdim);
    zeroDouble(errtr, (nclass + 1) * Ntree);
    
    if (labelts) {
        nclts  = (int *) mxCalloc(nclass, sizeof(int));
        for (n = 0; n < ntest; ++n) nclts[clts[n]-1]++;
        zeroDouble(errts, (nclass + 1) * Ntree);
    }
    //printf("labelts %d\n",labelts);fflush(stdout);
    if (imp) {
        zeroDouble(imprt, (nclass+2) * mdim);
        zeroDouble(impsd, (nclass+1) * mdim);
        if (localImp) zeroDouble(impmat, nsample * mdim);
    }
    if (iprox) {
        zeroDouble(prox, nsample0 * nsample0);
        if (*testdat) zeroDouble(proxts, ntest * (ntest + nsample0));
    }
	
	/* ====== Aggregating original dataset by synthesizing non-existing data ====== */
	createClass(x, nsample0, nsample, mdim, num_nv_src_discrete, Y_nv_src_discrete, num_nv_src_continuous, Y_nv_src_continuous);
	
	if(keep_pseudo_X)
	{
		int num_el = nsample0 * mdim;
		for (int idx = 0; idx < num_el; idx++)
		{
			pseudo_X[idx] = x[idx + num_el];
		}
	}	
	else
	{
		pseudo_X[0] = -1;
	}
	
	/* Order the feature values for each variable, giving $a$ and $b$ */
    makeA(x, mdim, nsample, cat, at, b);

    //R_CheckUserInterrupt();

    /* Starting the main loop over number of trees. */
    GetRNGstate();
    if (trace <= Ntree) 
	{
        /* Print header for running output. */
        Rprintf("ntree      OOB");
        for (n = 1; n <= nclass; ++n) Rprintf("%7i", n);
        if (labelts) 
		{
            Rprintf("|    Test");
            for (n = 1; n <= nclass; ++n) Rprintf("%7i", n);
        }
        Rprintf("\n");
    }
    idxByNnode = 0;
    idxByNsample = 0;
    
	time_t curr_time;
    
    /* ====== The associasion between visual and non-visual data sources ========= */
    int num_nv_src = num_nv_src_discrete + num_nv_src_continuous;
    double* surr_var_assoc = new double[mdim * (mdim + num_nv_src)]; // used in grwoing each tree for measuring the association between predictors/variables
	int* num_measured_assoc = new int[mdim * (mdim + num_nv_src)]; // store the number of measuring association for each pair in a tree
    int* num_measured_assoc_trees = new int[mdim * (mdim + num_nv_src)]; // store the number of measured association over all trees
    int* num_split_var = new int[mdim]; // The number of splits on each variable
    zeroDouble(meanSurrVarAssoc, mdim * (mdim + num_nv_src));
    zeroInt(num_measured_assoc_trees, mdim * (mdim + num_nv_src));
	int* indicator_LR_child = new int[nsample];
	zeroInt(indicator_LR_child, nsample);
	
	/* ====== To cope with missing nv-source labels ========= */
	double num_original_samples; // Only consider originial samples, since newly created ones have no nv-src labels 
	double* ratio_missing_nv_src = new double[num_nv_src]; // The ratio of samples with missing label to all samples
    int* num_missing_nv_src = new int[num_nv_src]; // The number of samples with missing label
	
	/* ====== Complexitity Analysis====== */
	*visual_info_gain = 0;
	*num_flat_region = 0;
	*num_split = 0;
	
	double visual_info_gain_tree = 0;
	double num_flat_region_tree = 0;
	double num_split_tree = 0;
	
    for(int idx_tree = 0; idx_tree < Ntree; idx_tree++) 
	{
        /* Do we need to simulate data for the second class? */
		/* Originally the following if block is uncommented, which is regarded as a potential bug */		
        //if (addClass) 
		//{
		//	createClass(x, nsample0, nsample, mdim, Y_nv_src, *num_nv_src);
		//}
		
        do // train a tree
		{
            zeroInt(nodestatus + idxByNnode, *nrnodes);
            zeroInt(treemap + 2*idxByNnode, 2 * *nrnodes);
            zeroDouble(xbestsplit + idxByNnode, *nrnodes);
            zeroInt(nodeclass + idxByNnode, *nrnodes);
			zeroInt(node_size_table + idxByNnode, *nrnodes);
            zeroInt(varUsed, mdim);
			
			/* ====== To cope with missing nv-source labels ========= */
			zeroDouble(ratio_missing_nv_src, num_nv_src);
			zeroInt(num_missing_nv_src, num_nv_src);
			num_original_samples = 0;
			
            /* TODO: Put all sampling code into a function. */
            /* drawSample(sampsize, nsample, ); */
            if (stratify) /* stratified sampling */
			{  
                zeroInt(jin, nsample);
                zeroDouble(tclasspop, nclass);
                zeroDouble(win, nsample);
                if (replace) {  /* with replacement */
                    for (n = 0; n < nstrata; ++n) {
                        for (j = 0; j < sampsize[n]; ++j) {
                            ktmp = (int) (unif_rand() * strata_size[n]);
                            k = strata_idx[n][ktmp];
                            tclasspop[cl[k] - 1] += classwt[cl[k] - 1];
                            win[k] += classwt[cl[k] - 1];
                            jin[k] = 1;
                        }
                    }
                } else { /* stratified sampling w/o replacement */
                    /* re-initialize the index array */
                    zeroInt(strata_size, nstrata);
                    for (j = 0; j < nsample; ++j) {
                        strata_size[strata[j] - 1] ++;
                        strata_idx[strata[j] - 1][strata_size[strata[j] - 1] - 1] = j;
                    }
                    /* sampling without replacement */
                    for (n = 0; n < nstrata; ++n) {
                        last = strata_size[n] - 1;
                        for (j = 0; j < sampsize[n]; ++j) {
                            ktmp = (int) (unif_rand() * (last+1));
                            k = strata_idx[n][ktmp];
                            swapInt(strata_idx[n][last], strata_idx[n][ktmp]);
                            last--;
                            tclasspop[cl[k] - 1] += classwt[cl[k]-1];
                            win[k] += classwt[cl[k]-1];
                            jin[k] = 1;
                        }
                    }
                }
            } 
			else /* unstratified sampling */
			{  
                anyEmpty = 0;
                ntry = 0;
                do {
                    zeroInt(jin, nsample);
                    zeroDouble(tclasspop, nclass);
                    zeroDouble(win, nsample);
					
					zeroDouble(sum_nv_src_continuous, num_nv_src_continuous);
					zeroDouble(tclasspop_nv_src_discrete, max_num_class_nv_src_discrete*num_nv_src_discrete);
					zeroDouble(win_nv_src, nsample);
					
                    if (replace) 
					{
                        for (n = 0; n < *sampsize; ++n) 
						{
                            k = unif_rand() * nsample; // should be [0 nsample-1]
                            tclasspop[cl[k] - 1] += classwt[cl[k]-1];
                            win[k] += classwt[cl[k]-1];
                            jin[k] = 1;
							
							/*Note that we only use non-visual data that are attached with the original data points */
							if(k < nsample0)
							{  
								num_original_samples += 1;
								
								/* Count the number of all discrete non-visual data. */
								for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
								{
									int y_nv = Y_nv_src_discrete[idx_nv * nsample + k];
									if(y_nv > 0)
									{
										tclasspop_nv_src_discrete[idx_nv * max_num_class_nv_src_discrete + (y_nv-1)] += 1;
									}
									else // the label is missing
									{
										num_missing_nv_src[idx_nv] += 1;
										//printf("The %d-th sample missing %d-th discrete nv-src data\n", k+1, idx_nv);
									}
								}
							
								/* Aggregate the response of all continuous non-visual data. */
								for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
								{
									double y_nv = Y_nv_src_continuous[idx_nv * nsample + k];
									if(y_nv != y_nv) // the label is missing
									{
										num_missing_nv_src[num_nv_src_discrete + idx_nv] += 1;
										//printf("The %d-th sample missing %d(%d)-th (continuous) nv-src data\n", k+1, num_nv_src_discrete + idx_nv, idx_nv);
									}
									else
									{
										sum_nv_src_continuous[idx_nv] += y_nv;
									}			
								}
								
								win_nv_src[k] += 1;
							}
                        }
                    }
					else 
					{
                        for (n = 0; n < nsample; ++n) nind[n] = n;
                        last = nsample - 1;
                        for (n = 0; n < *sampsize; ++n) 
						{
                            ktmp = (int) (unif_rand() * (last+1));
                            k = nind[ktmp];
                            swapInt(nind[ktmp], nind[last]);
                            last--;
                            tclasspop[cl[k] - 1] += classwt[cl[k]-1];
                            win[k] += classwt[cl[k]-1];
                            jin[k] = 1;
							
							/*Note that we only use non-visual data that are attached with the original data points */
							if(k < nsample0)
							{  
								num_original_samples += 1;
								
								/* Count the number of all discrete non-visual data. */
								for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
								{
									int y_nv = Y_nv_src_discrete[idx_nv * nsample + k];
									if(y_nv > 0)
									{
										tclasspop_nv_src_discrete[idx_nv * max_num_class_nv_src_discrete + (y_nv-1)] += 1;
									}
									else // the label is missing
									{
										num_missing_nv_src[idx_nv] += 1;
										//printf("The %d-th sample missing %d-th discrete nv-src data\n", k+1, idx_nv);
									}
								}
							
								/* Aggregate the response of all continuous non-visual data. */
								for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
								{
									double y_nv = Y_nv_src_continuous[idx_nv * nsample + k];
									if(y_nv != y_nv) // the label is missing
									{
										num_missing_nv_src[num_nv_src_discrete + idx_nv] += 1;
										//printf("The %d-th sample missing %d(%d)-th (continuous) nv-src data\n", k+1, num_nv_src_discrete + idx_nv, idx_nv);
									}
									else
									{
										sum_nv_src_continuous[idx_nv] += y_nv;
									}			
								}
								
								win_nv_src[k] += 1;
							}
                        }
                    }
					
                    /* Check if any class is missing in the sampled subset */
                    for (n = 0; n < nclass; ++n) 
					{
                        if (tclasspop[n] == 0) anyEmpty = 1;
                    }
                    ntry++;
                } while (anyEmpty && ntry <= 10);
            }
            
            /* If need to keep indices of inbag data, do that here. */
            if (keepInbag) 
			{
                for (n = 0; n < nsample0; ++n) 
				{
                    inbag[n + idxByNsample] = jin[n];
                }
            }
            
			/* Copy the original $a$ matrix back. */
            memcpy(a, at, sizeof(int) * mdim * nsample);
            modA(a, &nuse, nsample, mdim, cat, *maxcat, ncase, jin);
				
            // printf("$ Tree-#%d, nsample=%d, dimx1=%d, addClass=%d, num_boost_node=%d \n", idx_tree, nsample, dimx[1], addClass, nuse);
 		
            /* Initialize surr_var_assoc before using */
            zeroDouble(surr_var_assoc, mdim * (mdim + num_nv_src));
            zeroInt(num_measured_assoc, mdim * (mdim + num_nv_src));
			
			/* ====== To cope with missing nv-source labels ========= */
			if(num_original_samples == 0)
			{
				//printf("*************** No original samples are selected for growing the %d-th tree +++++++++++++++++++++++++", idx_tree);
				zeroDouble(ratio_missing_nv_src, num_nv_src);
			}
			else
			{
				//printf("--00-- num_original_samples = %f\n", num_original_samples);
				
				for(int idx_nv = 0; idx_nv < num_nv_src; idx_nv++)
				{
					int num_missing_cur_nv = num_missing_nv_src[idx_nv];
					ratio_missing_nv_src[idx_nv] = num_missing_cur_nv / num_original_samples;
					//printf("----- ratio_missing_nv_src[%d] = %f\n", idx_nv, ratio_missing_nv_src[idx_nv]);
				}
			}
						
			/*for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
			{
				printf("sum_nv_src_continuous: %f \n", sum_nv_src_continuous[idx_nv]);		
			}
			
			for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
			{
				for (int idx_cl = 0; idx_cl < max_num_class_nv_src_discrete; idx_cl++)
					printf("tclasspop_nv_src_discrete: %f\n", tclasspop_nv_src_discrete[idx_nv * max_num_class_nv_src_discrete + (idx_cl)]);
			}*/
			
            buildtree(a, b, cl, cat, maxcat, mdim, nsample, nclass,
                    treemap + 2*idxByNnode, bestvar + idxByNnode,
                    bestsplit, bestsplitnext, tgini,
                    nodestatus + idxByNnode, nodepop,
                    nodestart, classpop, tclasspop, tclasscat,
                    ta, *nrnodes, idmove, ndsize, ncase,
                    mtry, varUsed, nodeclass + idxByNnode,
                    ndbigtree + idx_tree, win, wr, wl, mdim,
                    nuse, mind, dimx[1],
					win_nv_src,
					num_nv_src_discrete,
					Y_nv_src_discrete,
					alpha_nv_src_discrete,
					max_num_class_nv_src_discrete,
					tclasspop_nv_src_discrete, classpop_nv_src_discrete_all_nodes,
					parent_crit_nv_src_discrete,
					pno_nv_src_discrete, pdo_nv_src_discrete, 
					rln_nv_src_discrete, rld_nv_src_discrete,
					rrn_nv_src_discrete, rrd_nv_src_discrete,
					wl_nv_src_discrete, wr_nv_src_discrete,
					num_nv_src_continuous,
					Y_nv_src_continuous,
					alpha_nv_src_continuous,
					parent_crit_nv_src_continuous,
					sum_nv_src_continuous_all_nodes, 
					sum_nv_src_continuous, 
					suml_nv_src_continuous, sumr_nv_src_continuous,
                    surr_var_assoc, num_measured_assoc, indicator_LR_child, 
					ratio_missing_nv_src, 
					&visual_info_gain_tree, &num_flat_region_tree, &num_split_tree,
					node_size_table + idxByNnode);	
            
            //printf("--- To update the association between sources\n");
            for(int idx = 0; idx < mdim * (mdim + num_nv_src); idx++)
            {
                // Divided by the times when var j is an split var when i is the optimal
                int num_measured = num_measured_assoc[idx];
                if(num_measured > 0)
                {
                    // To average over measuring times
                    meanSurrVarAssoc[idx] += surr_var_assoc[idx] / num_measured;
                    num_measured_assoc_trees[idx] += 1;
                }
            }
            
        /* if the "tree" has only the root node, start over */
        } while (ndbigtree[idx_tree] < 1);
		
		/*
			Translate split position back x-value
			(see rfutils.cpp)
		*/
        Xtranslate(x, mdim, *nrnodes, nsample, bestvar + idxByNnode,
                bestsplit, bestsplitnext, xbestsplit + idxByNnode,
                nodestatus + idxByNnode, cat, ndbigtree[idx_tree]);
        
        /* Aggregate information gain */
		//printf("====== visual_info_gain_tree = %f\n", visual_info_gain_tree);
		//printf("====== num_flat_region_tree = %f\n", num_flat_region_tree);
		//printf("====== num_split_tree = %f\n", num_split_tree);
		//printf("=================================== splits tried = %f\n", num_split_tree + num_flat_region_tree);
		*visual_info_gain += visual_info_gain_tree;
		*num_flat_region += num_flat_region_tree/Ntree;
		*num_split += num_split_tree/Ntree;
		
		/* Get test set error */
        if (*testdat) 
		{
            predictClassTree(xts, ntest, mdim, treemap + 2*idxByNnode,
							nodestatus + idxByNnode, xbestsplit + idxByNnode,
							bestvar + idxByNnode,
							nodeclass + idxByNnode, ndbigtree[idx_tree],
							cat, nclass, jts, nodexts, path_stat_xts, *maxcat);
							
			TestSetError(countts, jts, clts, outclts, ntest, nclass, idx_tree+1,
						errts + idx_tree*(nclass+1), labelts, nclts, cut);
        }
        
        /* Get out-of-bag predictions and errors. */
        predictClassTree(x, nsample, mdim, treemap + 2*idxByNnode,
						nodestatus + idxByNnode, xbestsplit + idxByNnode,
						bestvar + idxByNnode,
						nodeclass + idxByNnode, ndbigtree[idx_tree],
						cat, nclass, jtr, nodex, path_stat, *maxcat);
        
        zeroInt(nout, nclass);
        noutall = 0;
        for (n = 0; n < nsample; ++n) 
		{
            if (jin[n] == 0) 
			{
                /* increment the OOB votes */
                counttr[n*nclass + jtr[n] - 1] ++;
                /* count number of times a case is OOB */
                out[n]++;
                /* count number of OOB cases in the current iteration.
                 * nout[n] is the number of OOB cases for the n-th class.
                 * noutall is the number of OOB cases overall. */
                nout[cl[n] - 1]++;
                noutall++;
            }
        }
        
        /* Compute out-of-bag error rate. */
        oob(nsample, nclass, jin, cl, jtr, jerr, counttr, out,
                errtr + idx_tree*(nclass+1), outcl, cut);
        
        if ((idx_tree+1) % trace == 0) 
		{
            Rprintf("%5i: %6.2f%%", idx_tree+1, 100.0*errtr[idx_tree * (nclass+1)]);
            for (n = 1; n <= nclass; ++n) 
			{
                Rprintf("%6.2f%%", 100.0 * errtr[n + idx_tree * (nclass+1)]);
            }
            if (labelts) 
			{
                Rprintf("| ");
                for (n = 0; n <= nclass; ++n) 
				{
                    Rprintf("%6.2f%%", 100.0 * errts[n + idx_tree * (nclass+1)]);
                }
            }
            Rprintf("\n");
            
            //R_CheckUserInterrupt();
        }
        
        /*  DO VARIABLE IMPORTANCE  */
        if (imp) 
		{
            nrightall = 0;
            /* Count the number of correct prediction by the current tree
             * among the OOB samples, by class. */
            zeroInt(nright, nclass);
            for (n = 0; n < nsample; ++n) 
			{
                /* out-of-bag and predicted correctly: */
                if (jin[n] == 0 && jtr[n] == cl[n]) 
				{
                    nright[cl[n] - 1]++;
                    nrightall++;
                }
            }
            for (m = 0; m < mdim; ++m) 
			{
                if (varUsed[m]) 
				{
                    nrightimpall = 0;
                    zeroInt(nrightimp, nclass);
                    for (n = 0; n < nsample; ++n) 
						tx[n] = x[m + n*mdim];
                    /* Permute the m-th variable. */
                    permuteOOB(m, x, jin, nsample, mdim);
                    /* Predict the modified data using the current tree. */
                    predictClassTree(x, nsample, mdim, treemap + 2*idxByNnode,
                            nodestatus + idxByNnode,
                            xbestsplit + idxByNnode,
                            bestvar + idxByNnode,
                            nodeclass + idxByNnode, ndbigtree[idx_tree],
                            cat, nclass, jvr, nodex, path_stat, *maxcat);
                    /* Count how often correct predictions are made with
                     * the modified data. */
                    for (n = 0; n < nsample; n++) 
					{
                        if (jin[n] == 0) 
						{
                            if (jvr[n] == cl[n]) 
							{
                                nrightimp[cl[n] - 1]++;
                                nrightimpall++;
                            }
                            if (localImp && jvr[n] != jtr[n]) 
							{
                                if (cl[n] == jvr[n]) 
								{
                                    impmat[m + n*mdim] -= 1.0;
                                } 
								else 
								{
                                    impmat[m + n*mdim] += 1.0;
                                }
                            }
                        }
                        /* Restore the original data for that variable. */
                        x[m + n*mdim] = tx[n];
                    }
                    /* Accumulate decrease in proportions of correct
                     * predictions. */
                    for (n = 0; n < nclass; ++n) 
					{
                        if (nout[n] > 0) 
						{
                            imprt[m + n*mdim] += ((double) (nright[n] - nrightimp[n])) / nout[n];
                            impsd[m + n*mdim] += ((double) (nright[n] - nrightimp[n]) * (nright[n] - nrightimp[n])) / nout[n];
                        }
                    }
                    if (noutall > 0) 
					{
                        imprt[m + nclass*mdim] += ((double)(nrightall - nrightimpall)) / noutall;
                        impsd[m + nclass*mdim] += ((double) (nrightall - nrightimpall) * (nrightall - nrightimpall)) / noutall;
                    }
                }
            }
        }
        
        /*  DO PROXIMITIES */
        if (iprox) 
		{
            computeProximity(prox, oobprox, nodex, jin, oobpair, near);
            /* proximity for test data */
            if (*testdat) 
			{
				computeProximity(proxts, 0, nodexts, jin, oobpair, ntest);
				/* Compute proximity between testset and training set. */
				for (n = 0; n < ntest; ++n) 
				{
					for (k = 0; k < near; ++k) 
					{
						if (nodexts[n] == nodex[k])
							proxts[n + ntest * (k+ntest)] += 1.0;
					}
				}
            }
        }
        
        if (keepf) idxByNnode += *nrnodes;
        if (keepInbag) idxByNsample += nsample0;

		if(print_verbose_tree_progression)
		{
		#ifdef MATLAB
			time(&curr_time);
            mexPrintf("tree num %d created at %s", idx_tree, ctime(&curr_time));mexEvalString("drawnow;");
        #endif
		}
    
	}
	
    *visual_info_gain = *visual_info_gain / Ntree;
	
	
    // Average association for predictors/variables over all trees
    for(int idx = 0; idx < mdim * (mdim + num_nv_src); idx++)
    {
        // divided by the times wheh var j is an split var when i is the optimal
        int num_measured = num_measured_assoc_trees[idx];
        if(num_measured > 1)
        {
            meanSurrVarAssoc[idx] /= num_measured;
        }
    }
    
    PutRNGstate();
   
    
    /*  Final processing of variable importance. */
    for (m = 0; m < mdim; m++) tgini[m] /= Ntree;
      
    if (imp) 
	{
        for (m = 0; m < mdim; ++m) 
		{
            if (localImp) 
			{ /* casewise measures */
                for (n = 0; n < nsample; ++n) impmat[m + n*mdim] /= out[n];
            }
            /* class-specific measures */
            for (k = 0; k < nclass; ++k) 
			{
                av = imprt[m + k*mdim] / Ntree;
                impsd[m + k*mdim] = sqrt(((impsd[m + k*mdim] / Ntree) - av*av) / Ntree);
                imprt[m + k*mdim] = av;
                /* imprt[m + k*mdim] = (se <= 0.0) ? -1000.0 - av : av / se; */
            }
            /* overall measures */
            av = imprt[m + nclass*mdim] / Ntree;
            impsd[m + nclass*mdim] = sqrt(((impsd[m + nclass*mdim] / Ntree) - av*av) / Ntree);
            imprt[m + nclass*mdim] = av;
            imprt[m + (nclass+1)*mdim] = tgini[m];
        }
    } 
	else 
	{
        for (m = 0; m < mdim; ++m) imprt[m] = tgini[m];
    }
   
    /*  PROXIMITY DATA ++++++++++++++++++++++++++++++++*/
    if (iprox) 
	{
        for (n = 0; n < near; ++n) 
		{
            for (k = n + 1; k < near; ++k) 
			{
                prox[near*k + n] /= oobprox ?
                    (oobpair[near*k + n] > 0 ? oobpair[near*k + n] : 1) :
                        Ntree;
                        prox[near*n + k] = prox[near*k + n];
            }
            prox[near*n + n] = 1.0;
        }
        if (*testdat) 
		{
            for (n = 0; n < ntest; ++n)
			{
                for (k = 0; k < ntest + nsample; ++k)
                    proxts[ntest*k + n] /= Ntree;
                proxts[ntest*n + n]=1.0;
            }
        }
    }
    if (trace <= Ntree)
	{
        printf("\nmdim %d, nclass %d, nrnodes %d, nsample %d, ntest %d\n", mdim, nclass, *nrnodes, nsample, ntest);
        printf("\noobprox %d, mdim %d, nsample0 %d, Ntree %d, mtry %d, mimp %d", oobprox, mdim, nsample0, Ntree, mtry, mimp);
        printf("\nstratify %d, replace %d",stratify,replace);
        printf("\n");
    }
	
	/////////////////////////////////////////
	delete [] win_nv_src;
	
	delete [] parent_crit_nv_src_discrete;
	delete [] tclasspop_nv_src_discrete;
	delete [] classpop_nv_src_discrete_all_nodes;
	
	delete [] pno_nv_src_discrete;
	delete [] pdo_nv_src_discrete;
	
	delete [] rln_nv_src_discrete;
	delete [] rld_nv_src_discrete;
	delete [] rrn_nv_src_discrete;
	delete [] rrd_nv_src_discrete;
	delete [] wl_nv_src_discrete;
	delete [] wr_nv_src_discrete;
	
	
	delete [] sum_nv_src_continuous;
	delete [] suml_nv_src_continuous;
	delete [] sumr_nv_src_continuous;
	delete [] parent_crit_nv_src_continuous;
	delete [] sum_nv_src_continuous_all_nodes;
    
    //printf("��� To delete surr_var_assoc \n");
    delete [] surr_var_assoc;
    delete [] num_split_var;
    delete [] num_measured_assoc;
    delete [] num_measured_assoc_trees;
	delete [] indicator_LR_child;
	
	delete [] ratio_missing_nv_src; // The ratio of samples with missing label to all samples
    delete [] num_missing_nv_src; // The number of samples with missing label
    //printf("��� Delete surr_var_assoc is DONE\n");
	
	/////////////////////////////////////////
    
    //frees up the memory
    mxFree(tgini);
    mxFree(wl);
    mxFree(wr);
    mxFree(classpop);
    mxFree(tclasscat);
    mxFree(tclasspop);
    mxFree(tx);
    mxFree(win);
    mxFree(tp);
    mxFree(out);
    
    mxFree(bestsplitnext);
    mxFree(bestsplit);
    mxFree(nodepop);
    mxFree(nodestart);
    mxFree(jin);
    mxFree(nodex);
	mxFree(path_stat);
    mxFree(nodexts);
	mxFree(path_stat_xts);
    mxFree(ta);
    mxFree(ncase);
    mxFree(jerr);
    
    mxFree(varUsed);
    mxFree(jtr);
    mxFree(jvr);
    mxFree(classFreq);
    mxFree(jts);
    mxFree(idmove);
    mxFree(at);
    mxFree(a);
    mxFree(b);
    mxFree(mind);
    
    mxFree(nright);
    mxFree(nrightimp);
    mxFree(nout);
    
    if (oobprox) 
	{
        mxFree(oobpair);
    }
    
    if (stratify) 
	{
        mxFree(strata_size);
        for (n = 0; n < nstrata; ++n) 
		{
            mxFree(strata_idx[n]);
        }
        mxFree(strata_idx);        
    } 
	else 
	{
        if (replace)
            mxFree(nind);
    }
    //printf("labelts %d\n",labelts);fflush(stdout);
    fflush(stdout);
    if (labelts) 
	{
        mxFree(nclts);        
    }
    //printf("stratify %d",stratify);fflush(stdout);
}


/* predict the category of a new sample */
void classForest(int *mdim, int *ntest, int *nclass, int *maxcat,
		int *nrnodes, int *ntree, double *x, double *xbestsplit,
		double *pid, double *cutoff, double *countts, int *treemap,
		int *nodestatus, int *cat, int *nodeclass, int *jts,
		int *jet, int *bestvar, int *node, int *treeSize,
		int *keepPred, int *prox, double *proxMat, int *nodes,
		int *path_stat, int path_length) 
{
    int j, n, n1, n2, idxNodes, offset1, offset2, *junk, ntie;
    double crit, cmax;
    
    zeroDouble(countts, (*nclass) * (*ntest));
    idxNodes = 0;
    offset1 = 0;
    offset2 = 0;
    junk = NULL;
	int offset_path_stat = 0;

    for (j = 0; j < *ntree; ++j) 
	{
        /* predict by the j-th tree, call the function in file classTree.cpp */
        predictClassTree(x, *ntest, *mdim, treemap + 2*idxNodes,
                nodestatus + idxNodes, xbestsplit + idxNodes,
                bestvar + idxNodes, nodeclass + idxNodes,
                treeSize[j], cat, *nclass,
                jts + offset1, node + offset2, path_stat + offset_path_stat, *maxcat);
        
        /* accumulate votes: */
        for (n = 0; n < *ntest; ++n) 
		{
            countts[jts[n + offset1] - 1 + n * *nclass] += 1.0;
        }
        
        /* if desired, do proximities for this round */
        if (*prox) computeProximity(proxMat, 0, node + offset2, junk, junk,
                *ntest);
        idxNodes += *nrnodes;
        if (*keepPred) offset1 += *ntest;
        if (*nodes)    offset2 += *ntest;
		
		if(path_length)
		{
			offset_path_stat += *ntest;
		}
		
    }
    
    //Rprintf("ntest %d\n", *ntest);
    /* Aggregated prediction is the class with the maximum votes/cutoff */
    for (n = 0; n < *ntest; ++n) 
	{
        //Rprintf("Ap: ntest %d\n", *ntest);
        cmax = 0.0;
        ntie = 1;
        for (j = 0; j < *nclass; ++j) 
		{
            crit = (countts[j + n * *nclass] / *ntree) / cutoff[j];
            if (crit > cmax) 
			{
                jet[n] = j + 1;
                cmax = crit;
            }
            /* Break ties at random: */
            if (crit == cmax) 
			{
                ntie++;
                if (unif_rand() > 1.0 / ntie) jet[n] = j + 1;
            }
        }
    }
    
    //Rprintf("ntest %d\n", *ntest);
    /* if proximities requested, do the final adjustment
     * (division by number of trees) */
    
    //Rprintf("prox %d",*prox);
    if (*prox) 
	{
        //Rprintf("prox: ntest %d\n", *ntest);
        for (n1 = 0; n1 < *ntest; ++n1) 
		{
            for (n2 = n1 + 1; n2 < *ntest; ++n2) 
			{
                proxMat[n1 + n2 * *ntest] /= *ntree;
                proxMat[n2 + n1 * *ntest] = proxMat[n1 + n2 * *ntest];
            }
            proxMat[n1 + n1 * *ntest] = 1.0;
        }
    }
}

/*
 * Modified by A. Liaw 1/10/2003 (Deal with cutoff)
 * Re-written in C by A. Liaw 3/08/2004
 */
void oob(int nsample, int nclass, int *jin, int *cl, int *jtr, int *jerr,
        int *counttr, int *out, double *errtr, int *jest,
        double *cutoff) 
{
    int j, n, noob, *noobcl, ntie;
    double qq, smax, smaxtr;
    
    noobcl  = (int *) mxCalloc(nclass, sizeof(int));
    zeroInt(jerr, nsample);
    zeroDouble(errtr, nclass+1);
    
    noob = 0;
    for (n = 0; n < nsample; ++n) 
	{
        if (out[n]) 
		{
            noob++;
            noobcl[cl[n]-1]++;
            smax = 0.0;
            smaxtr = 0.0;
            ntie = 1;
            for (j = 0; j < nclass; ++j) 
			{
                qq = (((double) counttr[j + n*nclass]) / out[n]) / cutoff[j];
                if (j+1 != cl[n]) smax = (qq > smax) ? qq : smax;
                /* if vote / cutoff is larger than current max, re-set max and
                 * change predicted class to the current class */
                if (qq > smaxtr) 
				{
                    smaxtr = qq;
                    jest[n] = j+1;
                }
                /* break tie at random */
                if (qq == smaxtr) 
				{
                    ntie++;
                    if (unif_rand() > 1.0 / ntie) 
					{
                        smaxtr = qq;
                        jest[n] = j+1;
                    }
                }
            }
            if (jest[n] != cl[n]) 
			{
                errtr[cl[n]] += 1.0;
                errtr[0] += 1.0;
                jerr[n] = 1;
            }
        }
    }
    errtr[0] /= noob;
    for (n = 1; n <= nclass; ++n) errtr[n] /= noobcl[n-1];
    mxFree(noobcl);
}


void TestSetError(double *countts, int *jts, int *clts, int *jet, int ntest,
        int nclass, int nvote, double *errts,
        int labelts, int *nclts, double *cutoff) 
{
    int j, n, ntie;
    double cmax, crit;
    
    for (n = 0; n < ntest; ++n) 
		countts[jts[n]-1 + n*nclass] += 1.0;
    
    /*  Prediction is the class with the maximum votes */
    for (n = 0; n < ntest; ++n) 
	{
        cmax=0.0;
        ntie = 1;
        for (j = 0; j < nclass; ++j) 
		{
            crit = (countts[j + n*nclass] / nvote) / cutoff[j];
            if (crit > cmax) 
			{
                jet[n] = j+1;
                cmax = crit;
            }
            /*  Break ties at random: */
            if (crit == cmax) 
			{
                ntie++;
                if (unif_rand() > 1.0 / ntie) 
				{
                    jet[n] = j+1;
                    cmax = crit;
                }
            }
        }
    }
    if (labelts) 
	{
        zeroDouble(errts, nclass + 1);
        for (n = 0; n < ntest; ++n) 
		{
            if (jet[n] != clts[n]) 
			{
                errts[0] += 1.0;
                errts[clts[n]] += 1.0;
            }
        }
        errts[0] /= ntest;
        for (n = 1; n <= nclass; ++n) errts[n] /= nclts[n-1];
    }
}
