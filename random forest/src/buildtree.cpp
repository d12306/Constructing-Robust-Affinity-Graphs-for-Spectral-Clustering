#include <stdlib.h>
#include "math.h"
#include "fstream"
#include "iostream"
#include "qsort.h"
#include "buildtree.h"
#include "rfutils.h"

#ifndef MATLAB
#include "mex_stub.h"
#endif

#ifdef MATLAB
#include "mex.h"
#define Rprintf mexPrintf
#endif

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
/*     
     Copyright (C) 2001-7  Leo Breiman and Adele Cutler and Merck & Co, Inc.
     This program is free software; you can redistribute it and/or
     modify it under the terms of the GNU General Public License
     as published by the Free Software Foundation; either version 2
     of the License, or (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
     
     Modified by Andy Liaw and Matt Wiener:
     The main program is re-written as a C function to be called from R.
     All calls to the uniform RNG is replaced with R's RNG.  Some subroutines
     not called are excluded.  Variables and arrays declared as double as 
     needed.  Unused variables are deleted. (cons means constraints -- by Eddy)
     
     SUBROUTINE BUILDTREE
 */      

void buildtree(int* a, int* b, int* cl, int* cat, int* maxcat, int mdim, int nsample,                                
    int nclass, int* treemap, int* bestvar, int* bestsplit, int* bestsplitnext, double* tgini,                
    int *nodestatus, int *nodepop, int *nodestart, double *classpop, double *tclasspop,                        
    double* tclasscat, int* ta, int nrnodes, int *idmove, int ndsize, int *ncase, int mtry, int *iv,        
    int *nodeclass, int* ndbigtree, double *win, double *wr, double *wl, int mred, int nuse, int *mind,       
    int num_ori_samples,                                                                 
    double* win_nv_src,
    int num_nv_src_discrete,
    int* Y_nv_src_discrete,
    double* alpha_nv_src_discrete,
    int max_num_class_nv_src_discrete,			
    double* tclasspop_nv_src_discrete, double* classpop_nv_src_discrete_all_nodes,
    double* parent_crit_nv_src_discrete,
    double* pno_nv_src_discrete, double* pdo_nv_src_discrete, 
    double* rln_nv_src_discrete, double* rld_nv_src_discrete,
    double* rrn_nv_src_discrete, double* rrd_nv_src_discrete,
    double* wl_nv_src_discrete, double* wr_nv_src_discrete,
    int num_nv_src_continuous,
    double* Y_nv_src_continuous, 
    double* alpha_nv_src_continuous,
    double* parent_crit_nv_src_continuous,
    double* sum_nv_src_continuous_all_nodes, 
    double* sum_nv_src_continuous, 
    double* suml_nv_src_continuous, double* sumr_nv_src_continuous,
    double* surr_var_assoc, int* num_measured_assoc, int* indicator_LR_child,
	double* ratio_missing_nv_src, 
	double* visual_info_gain, double* num_flat_region, double* num_split,
	int* node_size_table) 		
{
	/* 
	Buildtree consists of repeated calls to two subroutines, Findbestsplit
	and Movedata.  Findbestsplit does just that--it finds the best split of
	the current node.  Movedata moves the data( in the split node right and
	left so that the data corresponding to each child node is contiguous.
	The buildtree bookkeeping is different from that in Friedman's original
	CART program.  ncur is the total number of nodes to date.
	nodestatus(k)=1 if the kth node has been split.  nodestatus(k)=2 if the
	node exists but has not yet been split, and =-1 of the node (is  end terminal.
	A node is terminal if its size is below a threshold value, or if it is
	all one cl  endass, or if all the x-values are equal.  If the curmaxcatrent node k
	is split, then its children are numbered ncur(+1 (left), and
	ncur+2(right), ncur increases to ncur+2 and bestvarthe next node to be split is
	numbered k+1.  When no more nodes can buildtreebe split, buildtree returns to the
	main program. */
	
	//clock_t startTime = clock();
	
	double xrand=0.0, decsplit=0.0;   
	int msplit=0, ntie=0, nbest=0;

	zeroInt(nodestatus, nrnodes);
	zeroInt(nodestart, nrnodes);
	zeroInt(nodepop, nrnodes);
	zeroDouble2(classpop, nclass, nrnodes);
	
	/* ====== Complexity Analysis ======= */
	*visual_info_gain = 0.0;
	*num_flat_region = 0.0;
	*num_split = 0;
	
	/* === Compute the weights of visual and non-visual sources === */
	double* alpha_nv_src_discrete_final = new double[num_nv_src_discrete];
	double* alpha_nv_src_continuous_final = new double[num_nv_src_continuous];
	
	double alpha_sum_nv_src = 0;
	for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
	{
		alpha_nv_src_discrete_final[idx_nv] = alpha_nv_src_discrete[idx_nv];
		alpha_sum_nv_src += alpha_nv_src_discrete[idx_nv];
		//printf("/////// %% alpha_nv_src_discrete_final[%d] = %f\n", idx_nv, alpha_nv_src_discrete_final[idx_nv]);
	}
	for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
	{
		alpha_nv_src_continuous_final[idx_nv] = alpha_nv_src_continuous[idx_nv];
		alpha_sum_nv_src += alpha_nv_src_continuous[idx_nv];
		//printf("/////// ** alpha_nv_src_continuous_final[%d] = %f\n", idx_nv, alpha_nv_src_continuous_final[idx_nv]);
	}
	double alpha_visual = 1 - alpha_sum_nv_src;
	
	if(alpha_visual < 0)
	{
		alpha_visual = 0;
		// Normalise the weights of non-visual data so that the sum of all data sources equals to 1
		for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
		{
			alpha_nv_src_discrete_final[idx_nv] /= alpha_sum_nv_src;
		}
		for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
		{
			alpha_nv_src_continuous_final[idx_nv] /= alpha_sum_nv_src;
		}
	}
	//printf("/////// alpha_visual = %f\n", alpha_visual);
	
	// To adjust source weights based on $ratio_missing_nv_src$, the ratio of missing nv-src data
	double weight_penalty_total = 0;
	int num_unused_nv_src = 0;
	for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
	{
		double ratio_missing = ratio_missing_nv_src[idx_nv];
		double weight_penalty = ratio_missing * alpha_nv_src_discrete_final[idx_nv];
		alpha_nv_src_discrete_final[idx_nv] -= weight_penalty;
		if(alpha_nv_src_discrete_final[idx_nv] == 0)
		{
			num_unused_nv_src += 1;
		}
		weight_penalty_total += weight_penalty;
	}
	for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
	{
		double ratio_missing = ratio_missing_nv_src[num_nv_src_discrete+idx_nv];
		double weight_penalty = ratio_missing * alpha_nv_src_continuous_final[idx_nv];
		alpha_nv_src_continuous_final[idx_nv] -= weight_penalty;
		if(alpha_nv_src_continuous_final[idx_nv] == 0)
		{
			num_unused_nv_src += 1;
		}
		weight_penalty_total += weight_penalty;
	}
	
	// Let all sources with weight > 0 share $weight_penalty_total$ equally
	if(weight_penalty_total > 0)
	{
		//printf("============ Need to adjust source weight");
		double shared_weight_piece = weight_penalty_total / (1 + num_nv_src_discrete + num_nv_src_continuous - num_unused_nv_src);
		
		alpha_visual += shared_weight_piece;
		//printf("/////// alpha_visual = %f\n", alpha_visual);
		for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
		{
			if(alpha_nv_src_discrete_final[idx_nv] > 0)
			{
				alpha_nv_src_discrete_final[idx_nv] += shared_weight_piece;
			}
			//printf("/////// %% alpha_nv_src_discrete_final[%d] = %f\n", idx_nv, alpha_nv_src_discrete_final[idx_nv]);
		}	
		for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
		{
			if(alpha_nv_src_continuous_final[idx_nv] > 0)
			{
				alpha_nv_src_continuous_final[idx_nv] += shared_weight_piece;
			}
			//printf("/////// ** alpha_nv_src_continuous_final[%d] = %f\n", idx_nv, alpha_nv_src_continuous_final[idx_nv]);
		}
	}
	
	int has_nv_src_mode = 1;
	if(num_unused_nv_src == num_nv_src_discrete + num_nv_src_continuous)
	{
		has_nv_src_mode = 0;
		//printf("------------- No non-visual source mode, and the traditional mode is used\n");
	}
	
	/* ==== Initialise the state of the root: 
			visual class population, discrete non-visual population, and continuous non-vsiual sum === */
	// -- For visual data, compute the caterogory histogram
	for(int j = 1; j <= nclass; j++)
	{
		classpop[0*nclass+(j-1)] = tclasspop[j-1];
	}
	
	
	if(has_nv_src_mode == 1)
	{
		// -- For discrete non-visual data, compute the caterogory histogram
		zeroDouble3(classpop_nv_src_discrete_all_nodes, max_num_class_nv_src_discrete, num_nv_src_discrete, nrnodes);
		for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
		{
			// Deal with each discrete source
			for (int idx_cl = 0; idx_cl < max_num_class_nv_src_discrete; idx_cl++)
			{
				classpop_nv_src_discrete_all_nodes[0*max_num_class_nv_src_discrete*num_nv_src_discrete + idx_nv*max_num_class_nv_src_discrete + idx_cl] = tclasspop_nv_src_discrete[idx_nv * max_num_class_nv_src_discrete + idx_cl];
			}
		}
	
		// -- For continuous non-visual data, regerssion
		zeroDouble2(sum_nv_src_continuous_all_nodes, num_nv_src_continuous, nrnodes);
		for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
		{
			sum_nv_src_continuous_all_nodes[0 * num_nv_src_continuous + idx_nv] = sum_nv_src_continuous[idx_nv];
		}
	}
	
	
	/* ==== Calculate the intrinsic Gini impurity (unweighted) for visual and nonvisal data at the root === */
	// -- For visual data (classification)
	double pno = 0;
	double pdo = 0;
	for(int j = 1; j <= nclass; j++)
	{
		pno = pno + tclasspop[j-1] * tclasspop[j-1];
		pdo = pdo + tclasspop[j-1];
	}
	double root_crit_visual = pno / pdo;
	//printf("The initial Gini of visual data: %f\n", root_crit_visual);
	
	double* root_crit_nv_src_discrete = new double[num_nv_src_discrete];
	double* root_crit_nv_src_continuous = new double[num_nv_src_continuous];
	if(has_nv_src_mode == 1)
	{
		// -- For discrete non-visual sources (classification)	
		zeroDouble(root_crit_nv_src_discrete, num_nv_src_discrete);
		for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
		{
			double pno_nv = 0;
			double pdo_nv = 0;
			for(int idx_cl = 0; idx_cl < max_num_class_nv_src_discrete; idx_cl++)
			{
				double pop_cl = tclasspop_nv_src_discrete[idx_nv * max_num_class_nv_src_discrete + idx_cl];
				pno_nv += pop_cl * pop_cl;
				pdo_nv += pop_cl;
			}
			root_crit_nv_src_discrete[idx_nv] = pno_nv / pdo_nv;
			//printf("The initial Gini of %d-th discrete nonvisual data: %f\n", idx_nv, root_crit_nv_src_discrete[idx_nv]);
		}
		
		// -- For continuous non-visual sources (Regression)
		zeroDouble(root_crit_nv_src_continuous, num_nv_src_continuous);
		for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
		{
			double sum_nv = sum_nv_src_continuous[idx_nv];
			root_crit_nv_src_continuous[idx_nv] = (sum_nv * sum_nv) / nuse;
			//printf("The initial Gini of %d-th continuous nonvisual data: %f\n", idx_nv, root_crit_nv_src_continuous[idx_nv]);
		}
	}
	
	/* ====== The associate between visual and non-visual data ========= */
    int num_nv_src = num_nv_src_discrete + num_nv_src_continuous;
	zeroDouble(surr_var_assoc, mdim*(mdim + num_nv_src)); //
    zeroInt(num_measured_assoc, mdim*(mdim + num_nv_src));
	double* purity_increase_nv_src = new double[num_nv_src];
    int* rank_nv_src = new int[num_nv_src];
    int* surr_best_splits = new int[mtry*2]; // (idx_var, split_position)
	
	/* === Start the iterative splitting process === */
	int ncur = 1;
	nodestart[0] = 1;
	nodepop[0] = nuse; // the number of samples used for training current tree 
	node_size_table[0] = nuse;
	nodestatus[0] = 2;
	
	for(int kbuild = 1; kbuild <= nrnodes; kbuild++)
	{        
		//printf("------- Build tree kbuild = %d \n", kbuild);
		if (kbuild > ncur) 
		{
			//printf("-------^^^^^^ Build tree --- kbuild > ncur = %d \n\n", ncur);
			break;
		}
		
		if (nodestatus[kbuild-1] != 2) 
		{
			continue;
		}
		
		/* ====== Initialize parameters before splitting points by calling findbestsplit ====== */
		// -- Deal with visual data: class histogram
		int ndstart = nodestart[kbuild-1];
		int ndend = ndstart + nodepop[kbuild-1] - 1;
		
		for(int j = 1; j <= nclass; j++)
		{
			tclasspop[j-1] = classpop[(kbuild-1)*nclass+(j-1)];
		}
		
		if(has_nv_src_mode == 1)
		{
			// -- Deal with discrete non-visual data: class histogram
			for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
			{
				for(int idx_cl = 0; idx_cl < max_num_class_nv_src_discrete; idx_cl++)
				{
					tclasspop_nv_src_discrete[idx_nv * max_num_class_nv_src_discrete + idx_cl] = classpop_nv_src_discrete_all_nodes[(kbuild-1)*max_num_class_nv_src_discrete*num_nv_src_discrete + idx_nv*max_num_class_nv_src_discrete + idx_cl];
				}
			}
			
			// -- Deal with continuous non-visual data: the sum of non-visual value
			for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
			{
				sum_nv_src_continuous[idx_nv] = sum_nv_src_continuous_all_nodes[(kbuild-1) * num_nv_src_continuous + idx_nv];
			}
		}
		
		/* ====== To find the optimal split ====== */
        int jstat = 0; // the state of node, splitting node or leaf
        findbestsplit(a, b, cl, mdim, nsample, nclass, cat, maxcat,          			\
					ndstart, ndend, tclasspop, tclasscat, &msplit, &decsplit,   		\
					&nbest, ncase, &jstat, mtry, win, wr, wl, mred, mind,              	\
					num_ori_samples,
					root_crit_visual, root_crit_nv_src_discrete, root_crit_nv_src_continuous,
					has_nv_src_mode,
					alpha_visual,
					win_nv_src, 
					num_nv_src_discrete,
					Y_nv_src_discrete,
					alpha_nv_src_discrete_final,
					max_num_class_nv_src_discrete,
					tclasspop_nv_src_discrete, classpop_nv_src_discrete_all_nodes,
					parent_crit_nv_src_discrete,
					pno_nv_src_discrete, pdo_nv_src_discrete, 
					rln_nv_src_discrete, rld_nv_src_discrete,
					rrn_nv_src_discrete, rrd_nv_src_discrete,
					wl_nv_src_discrete, wr_nv_src_discrete,
					num_nv_src_continuous,
					Y_nv_src_continuous,
					alpha_nv_src_continuous_final,
					parent_crit_nv_src_continuous,
					sum_nv_src_continuous_all_nodes, 
					sum_nv_src_continuous, 
					suml_nv_src_continuous, sumr_nv_src_continuous,
                    surr_var_assoc, num_measured_assoc, indicator_LR_child,
                    purity_increase_nv_src, rank_nv_src, surr_best_splits,
					visual_info_gain, num_flat_region, num_split);			
		
		
		
		if (jstat == -1) // if the node is a leaf
		{
			nodestatus[kbuild-1] = -1;
			continue;
		}
		else
		{
			bestvar[kbuild-1] = msplit;
			iv[msplit-1] = 1;
			if (decsplit < 0.0) decsplit = 0.0;
			tgini[msplit-1] = tgini[msplit-1] + decsplit; // Gini importance of one feature is estimated by aggregating this
			if (cat[msplit-1] == 1) 
			{
				bestsplit[kbuild-1] = a[(nbest-1)*mdim+(msplit-1)];
				bestsplitnext[kbuild-1] = a[((nbest-1)+1)*mdim+(msplit-1)];
			}
			else
			{
				bestsplit[kbuild-1] = nbest;
				bestsplitnext[kbuild-1] = 0;
			}
		}
        
		int ndendl = 0;
		/* 
			- ncase: a copy of sample index list on the $msplit$ variable
		*/
        movedata(a, ta, mdim, nsample, ndstart, ndend, idmove, ncase, msplit, cat, nbest, &ndendl);
         
		// leftnode no. = ncur+1, rightnode no. = ncur+2.
		nodepop[(ncur-1)+1] = ndendl - ndstart + 1;
		nodepop[(ncur-1)+2] = ndend - ndendl;
		node_size_table[(ncur-1)+1] = nodepop[(ncur-1)+1];
		node_size_table[(ncur-1)+2] = nodepop[(ncur-1)+2];
		nodestart[(ncur-1)+1] = ndstart;
		nodestart[(ncur-1)+2] = ndendl + 1;

		/* ======= find the response distribution for the left and right child =======*/
		for(int n = ndstart; n <= ndendl; n++) // the left child
		{
			// For visual data - class population
			int nc = ncase[n-1]; // the index of current sample/case
			int j = cl[nc-1];
			classpop[((ncur-1)+1)*nclass+(j-1)] = classpop[((ncur-1)+1)*nclass+(j-1)] + win[nc-1];
			
			/* Note: only use the non-visual data attached with the original data points */
			if (nc <= num_ori_samples)
			{	
				// For discrete non-visual data - class population
				for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
				{
					int y_nv = Y_nv_src_discrete[idx_nv * nsample + (nc-1)];
					if(y_nv > 0)
					{
						// weighted count of non-visual data label
						classpop_nv_src_discrete_all_nodes[((ncur-1)+1)*max_num_class_nv_src_discrete*num_nv_src_discrete + idx_nv*max_num_class_nv_src_discrete + (y_nv-1)] += 1 * win_nv_src[nc-1];
					}
					else // the label is missing
					{
						//printf("The %d-th sample missing %d-th discrete nv-src data\n", nc, idx_nv);
					}
				}
			
				// For continuous non-visual data - regression
				for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
				{
					double y_nv = Y_nv_src_continuous[idx_nv*nsample+(nc-1)];
					if(y_nv != y_nv) // the label is missing
					{
						//printf("The %d-th sample missing %d(%d)-th (continuous) nv-src data\n", nc, num_nv_src_discrete + idx_nv, idx_nv);
					}
					else
					{
						// weighted sum of non-visual data
						double y_nv_weighted = y_nv * win_nv_src[nc-1];
						//printf("y_nv_weighted = %f, y_nv = %f, weight = %f\n", y_nv_weighted, y_nv, win_nv_src[nc-1]);
						sum_nv_src_continuous_all_nodes[((ncur-1)+1) * num_nv_src_continuous + idx_nv] += y_nv_weighted;
					}
				}
			}
		}
		
		for(int n = ndendl+1; n <= ndend; n++) // the right child
		{
			// For visual data
			int nc = ncase[n-1];
			int j = cl[nc-1];
			classpop[((ncur-1)+2)*nclass+(j-1)] = classpop[((ncur-1)+2)*nclass+(j-1)] + win[nc-1];
			
			/* Note: only use the non-visual data attached with the original data points */
			if (nc <= num_ori_samples)
			{	
				// For discrete non-visual data - class populations
				for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
				{
					//int nv_y = Y_nv_src_discrete[idx_nv * nsample + (nc-1)];
					//classpop_nv_src_discrete_all_nodes[((ncur-1)+2)*max_num_class_nv_src_discrete*num_nv_src_discrete + idx_nv*max_num_class_nv_src_discrete + (nv_y-1)] += 1 * win_nv_src[nc-1];
					
					int y_nv = Y_nv_src_discrete[idx_nv * nsample + (nc-1)];
					if(y_nv > 0)
					{
						// weighted count of non-visual data label
						classpop_nv_src_discrete_all_nodes[((ncur-1)+2)*max_num_class_nv_src_discrete*num_nv_src_discrete + idx_nv*max_num_class_nv_src_discrete + (y_nv-1)] += 1 * win_nv_src[nc-1];
					}
					else // the label is missing
					{
						//printf("The %d-th sample missing %d-th discrete nv-src data\n", nc, idx_nv);
					}
				}
			
				// For continuous non-visual data - regression
				for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
				{
					//double nv_y = Y_nv_src_continuous[idx_nv*nsample+(nc-1)];
					//sum_nv_src_continuous_all_nodes[((ncur-1)+2) * num_nv_src_continuous + idx_nv] += nv_y * win_nv_src[nc-1];
					
					double y_nv = Y_nv_src_continuous[idx_nv*nsample+(nc-1)];
					if(y_nv != y_nv) // the label is missing
					{
						//printf("The %d-th sample missing %d(%d)-th (continuous) nv-src data\n", nc, num_nv_src_discrete + idx_nv, idx_nv);
					}
					else
					{
						// weighted sum of non-visual data
						double y_nv_weighted = y_nv * win_nv_src[nc-1];
						//printf("y_nv_weighted = %f, y_nv = %f, weight = %f\n", y_nv_weighted, y_nv, win_nv_src[nc-1]);
						sum_nv_src_continuous_all_nodes[((ncur-1)+2) * num_nv_src_continuous + idx_nv] += y_nv_weighted;
					}
				}
			}
		}
         
		// Check on nodestatust classpop
		nodestatus[(ncur-1)+1] = 2;
		nodestatus[(ncur-1)+2] = 2;
		
		// 1. -- One stop criterion: number of samples <= ndsize 
		if(nodepop[(ncur-1)+1] <= ndsize)
		{
			nodestatus[(ncur-1)+1] = -1;
		}
		if(nodepop[(ncur-1)+2] <= ndsize)
		{
			nodestatus[(ncur-1)+2] = -1;
		}
		
		treemap[(kbuild-1)*2+(1-1)] = ncur + 1;
		treemap[(kbuild-1)*2+(2-1)] = ncur + 2;
		nodestatus[kbuild-1] = 1;
		
		ncur = ncur + 2;
		
		if (ncur >= nrnodes) 
		{
			//printf("----------build tree ------ reach the max number of nodes (%d)\n", nrnodes);
			break;
		}
    } 
	
	*ndbigtree = nrnodes;
	for (int k = nrnodes; k >= 1; k--)
	{
		if (nodestatus[k-1] == 0) 
		{
			*ndbigtree = *ndbigtree - 1;
		}
		if (nodestatus[k-1] == 2) 
		{
			nodestatus[k-1] = -1;
		}
	}
	
	/* === Form prediction in terminal nodes === */
	for( int kn = 1; kn <= (*ndbigtree); kn++)
	{
		if(nodestatus[kn-1] == -1) 
		{
            double pp = 0.0;
            ntie = 1;
            for(int j = 1; j <= nclass; j++)
			{
				if (classpop[(kn-1)*nclass+(j-1)] > pp) 
				{
					nodeclass[kn-1] = j;
					pp = classpop[(kn-1)*nclass+(j-1)];
				}
				// Break ties at random:
				if (classpop[(kn-1)*nclass+(j-1)] == pp) 
				{
					ntie = ntie + 1;
					rrand(&xrand);   //library function
					if (xrand <  (1.0 / ntie)) 
					{
						nodeclass[kn-1] = j;
						pp = classpop[(kn-1)*nclass+(j-1)];
					}
				}
			}
	    }
    }
	
	delete [] root_crit_nv_src_discrete;
	delete [] root_crit_nv_src_continuous;
    
    delete [] purity_increase_nv_src;
    delete [] rank_nv_src;
    delete [] surr_best_splits;
	
	delete [] alpha_nv_src_discrete_final;
	delete [] alpha_nv_src_continuous_final;
	
	//double split_time = double(clock() - startTime);
	//printf(">>>>>> time elapsed = %f\n", split_time);
}

/*
	Funcition FINDBESTSPLIT
	For the best split, msplit is the variable split on. decsplit is the
	dec. in impurity.  If msplit is numerical, nsplit is the case number
	of value of msplit split on, and nsplitnext is the case number of the
	next larger value of msplit.  If msplit is categorical, then nsplit is
	the coding into an integer of the categories going left.
*/
void findbestsplit(int* a, int* b, int* cl, int mdim, int nsample, int nclass, int* cat,           \
	int* maxcat, int ndstart, int ndend, double* tclasspop, double* tclasscat, int* msplit,                       \
	double* decsplit, int* nbest, int* ncase, int* jstat, int mtry, double* win, double* wr, double* wl,              \
	int mred, int* mind, 
	int num_ori_samples,
	double root_crit_visual, double* root_crit_nv_src_discrete, double* root_crit_nv_src_continuous,
	int has_nv_src_mode,
	double alpha_visual,
	double* win_nv_src,
	int num_nv_src_discrete,
	int* Y_nv_src_discrete,
	double* alpha_nv_src_discrete,
	int max_num_class_nv_src_discrete,
	double* tclasspop_nv_src_discrete, double* classpop_nv_src_discrete_all_nodes,
	double* parent_crit_nv_src_discrete,
	double* pno_nv_src_discrete, double* pdo_nv_src_discrete, 
	double* rln_nv_src_discrete, double* rld_nv_src_discrete,
	double* rrn_nv_src_discrete, double* rrd_nv_src_discrete,
	double* wl_nv_src_discrete, double* wr_nv_src_discrete,
	int num_nv_src_continuous,
	double* Y_nv_src_continuous,
	double* alpha_nv_src_continuous,
	double* parent_crit_nv_src_continuous,
	double* sum_nv_src_continuous_all_nodes, 
	double* sum_nv_src_continuous, 
	double* suml_nv_src_continuous, double* sumr_nv_src_continuous,
	double* surr_var_assoc, int* num_measured_assoc, int* indicator_LR_child,
    double* purity_increase_nv_src, int* rank_nv_src, int* surr_best_splits,
	double* visual_info_gain, double* num_flat_region, double* num_split) 
{
	int ncmax = 10, nn = 0, n_old_data = 0, ntie=0, nhit=0;
	int ncsplit = 512;
	
	double pno = 0.0, pdo = 0.0, xrand=0.0, score=0.0;
	double dn[32];
		
	*jstat = 0;
	*msplit = 0;
	double visual_gain_max = 0.0;
	
	/* === Compute the Gini impurity of current parent node before splitting: both visual and non-visual data === */
	// -- The total Gini critieria
	double parent_crit = 0;
	
	// -- For visual data (two-class in case of clusering)
	for( int j = 1; j <= nclass; j++)
	{
		pno = pno + tclasspop[j-1] * tclasspop[j-1];
		pdo = pdo + tclasspop[j-1];
	}
	double parent_crit_visual = pno / pdo;
	parent_crit += parent_crit_visual * alpha_visual;
	
	if(has_nv_src_mode == 1)
	{
		// -- For discrete non-visual data 
		zeroDouble(pno_nv_src_discrete, num_nv_src_discrete);
		zeroDouble(pdo_nv_src_discrete, num_nv_src_discrete);
		for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
		{	
			// not weighted
			double parent_crit_nv_src = 0;
			for( int idx_cl = 0; idx_cl < max_num_class_nv_src_discrete; idx_cl++)
			{
				double pop_cl = tclasspop_nv_src_discrete[idx_nv*max_num_class_nv_src_discrete + idx_cl];
				pno_nv_src_discrete[idx_nv] += pop_cl * pop_cl;
				pdo_nv_src_discrete[idx_nv] += pop_cl;
			}
			if(pdo_nv_src_discrete[idx_nv] > 0)
			{
				parent_crit_nv_src = pno_nv_src_discrete[idx_nv] / pdo_nv_src_discrete[idx_nv];
			}
			parent_crit_nv_src_discrete[idx_nv] = parent_crit_nv_src;
			// printf("-- parent_crit_nv_src_discrete[%d] = %f \n", idx_nv, parent_crit_nv_src_discrete[idx_nv]);
			// Added into $parent_crit$
			parent_crit += parent_crit_nv_src * alpha_nv_src_discrete[idx_nv];
		}
		
		
		// -- For continous non-visual data
		double nodecnt = (double) (ndend - ndstart + 1);
		for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
		{
			double parent_crit_nv_src = sum_nv_src_continuous[idx_nv] / nodecnt;
			parent_crit_nv_src *= sum_nv_src_continuous[idx_nv];
			// not weighted
			parent_crit_nv_src_continuous[idx_nv] = parent_crit_nv_src;
			// printf("-- parent_crit_nv_src_continuous[%d] = %f, sum_nv_src_continuous[%d] = %f \n", idx_nv, parent_crit_nv_src_continuous[idx_nv], idx_nv, sum_nv_src_continuous[idx_nv]);
			// Added into $parent_crit$
			parent_crit += parent_crit_nv_src * alpha_nv_src_continuous[idx_nv];
		}
	}	
	
	
	/* === initialise the data structure for keeping purity increase of non-visual data === */
    int num_nv_src = num_nv_src_discrete + num_nv_src_continuous;
	zeroDouble(purity_increase_nv_src, num_nv_src);
    zeroInt(surr_best_splits, mtry * 2);
    
	/* === Start the main loop through variables to find the best split === */
	double critmax = -1.0e25;  
	for( int k = 1; k <= mred; k++)
	{
		mind[k-1] = k;
	}
	nn = mred;
	
	// sampling mtry variables w/o replacement.
	for(int mt = 1; mt <= mtry; mt++)
	{
		// 1. Randomly select a variable/dimension
		rrand(&xrand);   // a library function
		int j = (int)(nn * xrand) + 1;  
		int mvar = mind[j-1]; // start from 1 
		mind[j-1] = mind[nn-1];
		mind[nn-1] = mvar;
		nn = nn - 1;
		int lcat = cat[mvar-1];
		
		surr_best_splits[(mt-1)*2] = mvar;
		double critmax_var = -1.0e25; // the max crit for current selected variable/predictor
		if (lcat == 1) // Split on a numerical predictor.
		{
			/* ====== 1. Do initialization: all nodes are sent into the right child ====== */
			// -- Deal with visual data (Classification)
			double rrn = pno;
			double rrd = pdo;
			double rln = 0.0;     
			double rld = 0.0;
			
			zeroDouble(wl, nclass);
			for(int j = 1; j <= nclass; j++)
			{
				wr[j-1] = tclasspop[j-1];
			}
			
			if(has_nv_src_mode == 1)
			{
				// -- Deal with discrete non-visual data
				
				// **** Deal with numerators and denominator 
				zeroDouble(rln_nv_src_discrete, num_nv_src_discrete);
				zeroDouble(rld_nv_src_discrete, num_nv_src_discrete);
				zeroDouble(rrn_nv_src_discrete, num_nv_src_discrete);
				zeroDouble(rrd_nv_src_discrete, num_nv_src_discrete);
				for (int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
				{
					rrn_nv_src_discrete[idx_nv] = pno_nv_src_discrete[idx_nv];
					rrd_nv_src_discrete[idx_nv] = pdo_nv_src_discrete[idx_nv];
				}
				// **** Deal with weights
				zeroDouble(wl_nv_src_discrete, num_nv_src_discrete * max_num_class_nv_src_discrete);
				for (int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
				{
					for(int idx_cl = 0; idx_cl < max_num_class_nv_src_discrete; idx_cl++)
					{
						wr_nv_src_discrete[idx_nv*max_num_class_nv_src_discrete + idx_cl] = tclasspop_nv_src_discrete[idx_nv*max_num_class_nv_src_discrete + idx_cl];
					}
				}
			
				// -- Deal with continuous non-visual data
				for (int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
				{					
					suml_nv_src_continuous[idx_nv] = 0;
					sumr_nv_src_continuous[idx_nv] = sum_nv_src_continuous[idx_nv];
				}
			}
			
			/* ====== 2. Loop: search the best splitting on the selected variable over the coming data points ====== */
			ntie = 1;
			for(int nsp = ndstart; nsp <= (ndend-1); nsp++)
			{
				// ---- Update varibles about visual features for computing Gini impurity
				int nc = a[(nsp-1)*mdim+(mvar-1)]; // the index of the chosen sample
				double u = win[nc-1];
				int k = cl[nc-1];
				rln = rln + u * (2 * wl[k-1] + u);
				rrn = rrn + u * ((-2) * wr[k-1] + u);
				rld = rld + u;
				rrd = rrd - u;
				wl[k-1] = wl[k-1] + u;
				wr[k-1] = wr[k-1] - u;
				
				if(has_nv_src_mode == 1 && nc <= num_ori_samples) 
				{	
					double u_nv = win_nv_src[(nc-1)]; // the weight of current sample
					
					// ---- Update varibles about discrete non-visual data for computing Gini impurity
					for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
					{
						int k_nv = Y_nv_src_discrete[idx_nv * nsample + (nc-1)]; // class label of nv-src
						if(k_nv > 0) 
						{
							// weighted
							rln_nv_src_discrete[idx_nv] += u_nv * (2 * wl_nv_src_discrete[idx_nv*max_num_class_nv_src_discrete + k_nv-1] + u_nv);
							rrn_nv_src_discrete[idx_nv] += u_nv * ((-2) * wr_nv_src_discrete[idx_nv*max_num_class_nv_src_discrete + k_nv-1] + u_nv);
							rld_nv_src_discrete[idx_nv] += u_nv;
							rrd_nv_src_discrete[idx_nv] -= u_nv;
							wl_nv_src_discrete[idx_nv*max_num_class_nv_src_discrete + k_nv-1] += u_nv;
							wr_nv_src_discrete[idx_nv*max_num_class_nv_src_discrete + k_nv-1] -= u_nv;
						}
						else // the label is missing
						{
							//printf("The %d-th sample missing %d-th discrete nv-src data\n", nc, idx_nv+1);
						}
					}
					
					// ---- Update varibles about continuous non-visual data for computing Gini impurity
					for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
					{
						// weighted
						double y_nv = Y_nv_src_continuous[idx_nv * nsample + (nc-1)] * u_nv;
						if(y_nv != y_nv) // the label is missing
						{
							//printf("The %d-th sample missing %d(%d)-th (continuous) nv-src data\n", nc, num_nv_src_discrete + idx_nv+1, idx_nv+1);
						}
						else
						{
							suml_nv_src_continuous[idx_nv] += y_nv;
							sumr_nv_src_continuous[idx_nv] -= y_nv;
						}
					}
				}
				
				// if the next sample has the same value with that of current split sample on the selected variable
				int abuf = a[((nsp-1) + 1)*mdim + (mvar-1)];		   
				if (b[(nc-1)*mdim+(mvar-1)] < b[(abuf-1)*mdim + (mvar-1)]) // not equal, <
				{
					*num_split += 1;
					// If neither nodes is empty, compute the information gain of current split.
					if (min(rrd, rld) > 1.0e-5)   
					{
						/* ====== The total Gini (purity), larger is better ====== */
						double crit = 0;
						
						// Compute the Gini for both the left and right child using different sources
						
						// ------ With visual data
						double crit_visual = (rln / rld) + (rrn / rrd);
						crit_visual -= parent_crit_visual;
						crit_visual /= root_crit_visual; // Normalisation
						double visual_gain_split = crit_visual;
						crit_visual *= alpha_visual; // weighting
						
						crit += crit_visual;
						
						if(has_nv_src_mode == 1)
						{
							// ------- With discrete non-visual data
							for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
							{
								if(alpha_nv_src_discrete[idx_nv] <= 0) continue;
								
								double crit_nv_src_left = 0;
								double crit_nv_src_right = 0;
								if(rld_nv_src_discrete[idx_nv] > 0)
								{
									crit_nv_src_left = rln_nv_src_discrete[idx_nv] / rld_nv_src_discrete[idx_nv];
								}
								if(rrd_nv_src_discrete[idx_nv] > 0) 
								{
									crit_nv_src_right = rrn_nv_src_discrete[idx_nv] / rrd_nv_src_discrete[idx_nv];
								}
								double crit_nv_src = crit_nv_src_left + crit_nv_src_right;
								crit_nv_src -= parent_crit_nv_src_discrete[idx_nv];
								crit_nv_src /= root_crit_nv_src_discrete[idx_nv];   // Normalisation
								crit_nv_src *= alpha_nv_src_discrete[idx_nv];       // weighting
								crit += crit_nv_src;
							}
							
							// ------- With continuous non-visual data
							for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
							{
								if(alpha_nv_src_continuous[idx_nv] <= 0) continue;
								
								double crit_nv_src_left = suml_nv_src_continuous[idx_nv] * suml_nv_src_continuous[idx_nv] / (nsp - ndstart + 1);
								double crit_nv_src_right = sumr_nv_src_continuous[idx_nv] * sumr_nv_src_continuous[idx_nv] / (ndend - nsp);
								double crit_nv_src = crit_nv_src_left + crit_nv_src_right;
								crit_nv_src -= parent_crit_nv_src_continuous[idx_nv];
								crit_nv_src /= root_crit_nv_src_continuous[idx_nv];     // Normalisation
								crit_nv_src *= alpha_nv_src_continuous[idx_nv];         // weighting
								crit += crit_nv_src;
							}
						}
						// Update the critmax
						if (crit > critmax) 
						{
							*nbest = nsp;
							critmax = crit;
							*msplit = mvar;
							visual_gain_max = visual_gain_split;
							
							// Update the critmax_var, for computing surr_var_assoc
							critmax_var = crit;
							surr_best_splits[(mt-1)*2+1] = nsp;
							
							if(has_nv_src_mode == 1)
							{
								// ------- With discrete non-visual data
								for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
								{								
									double crit_nv_src_left = 0;
									double crit_nv_src_right = 0;
									if(rld_nv_src_discrete[idx_nv] > 0)
									{
										crit_nv_src_left = rln_nv_src_discrete[idx_nv] / rld_nv_src_discrete[idx_nv];
									}
									if(rrd_nv_src_discrete[idx_nv] > 0) 
									{
										crit_nv_src_right = rrn_nv_src_discrete[idx_nv] / rrd_nv_src_discrete[idx_nv];
									}
									double crit_nv_src = crit_nv_src_left + crit_nv_src_right;
									crit_nv_src -= parent_crit_nv_src_discrete[idx_nv];
									crit_nv_src /= root_crit_nv_src_discrete[idx_nv];   // Normalisation
									crit_nv_src *= alpha_nv_src_discrete[idx_nv];       // weighting
									purity_increase_nv_src[idx_nv] = crit_nv_src;
								}
								
								// ------- With continuous non-visual data
								for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
								{
									double crit_nv_src_left = suml_nv_src_continuous[idx_nv] * suml_nv_src_continuous[idx_nv] / (nsp - ndstart + 1);
									double crit_nv_src_right = sumr_nv_src_continuous[idx_nv] * sumr_nv_src_continuous[idx_nv] / (ndend - nsp);
									double crit_nv_src = crit_nv_src_left + crit_nv_src_right;
									crit_nv_src -= parent_crit_nv_src_continuous[idx_nv];
									crit_nv_src /= root_crit_nv_src_continuous[idx_nv];     // Normalisation
									crit_nv_src *= alpha_nv_src_continuous[idx_nv];         // weighting
									purity_increase_nv_src[num_nv_src_discrete+idx_nv] = crit_nv_src;
								}
							}	
						}
						else
						{	
							// Update the critmax_var, for computing surr_var_assoc
							if (crit > critmax_var)
							{
								critmax_var = crit;
								surr_best_splits[(mt-1)*2+1] = nsp;
							}
						}
						
						// Break ties at random:
						if (crit == critmax) 
						{
							ntie = ntie + 1;
							rrand(&xrand);   
							if (xrand < (1.0 / ntie)) 
							{
								*nbest = nsp;
								critmax = crit;
								*msplit = mvar;
								
								// Update the critmax_var, for computing surr_var_assoc
								critmax_var = crit;
								surr_best_splits[(mt-1)*2+1] = nsp;
								if(has_nv_src_mode == 1)
								{
									// ------- With discrete non-visual data
									for(int idx_nv = 0; idx_nv < num_nv_src_discrete; idx_nv++)
									{								
										double crit_nv_src_left = 0;
										double crit_nv_src_right = 0;
										if(rld_nv_src_discrete[idx_nv] > 0)
										{
											crit_nv_src_left = rln_nv_src_discrete[idx_nv] / rld_nv_src_discrete[idx_nv];
										}
										if(rrd_nv_src_discrete[idx_nv] > 0) 
										{
											crit_nv_src_right = rrn_nv_src_discrete[idx_nv] / rrd_nv_src_discrete[idx_nv];
										}
										double crit_nv_src = crit_nv_src_left + crit_nv_src_right;
										crit_nv_src -= parent_crit_nv_src_discrete[idx_nv];
										crit_nv_src /= root_crit_nv_src_discrete[idx_nv];   // Normalisation										
										crit_nv_src *= alpha_nv_src_discrete[idx_nv];       // weighting
										purity_increase_nv_src[idx_nv] = crit_nv_src;
									}
									
									// ------- With continuous non-visual data
									for(int idx_nv = 0; idx_nv < num_nv_src_continuous; idx_nv++)
									{
										double crit_nv_src_left = suml_nv_src_continuous[idx_nv] * suml_nv_src_continuous[idx_nv] / (nsp - ndstart + 1);
										double crit_nv_src_right = sumr_nv_src_continuous[idx_nv] * sumr_nv_src_continuous[idx_nv] / (ndend - nsp);
										double crit_nv_src = crit_nv_src_left + crit_nv_src_right;
										crit_nv_src -= parent_crit_nv_src_continuous[idx_nv];
										crit_nv_src /= root_crit_nv_src_continuous[idx_nv];     // Normalisation
										crit_nv_src *= alpha_nv_src_continuous[idx_nv];         // weighting
										purity_increase_nv_src[num_nv_src_discrete+idx_nv] = crit_nv_src;
									}
								}
							}
						}
					}
				}
				else
				{
					*num_flat_region += 1;
				}
			}
		}
		else // Split on a categorical predictor.
		{
			// printf(">>>>>>>>>>>>>>>>>>>>> Split on a categorical predictor <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
			zeroDouble2(tclasscat, nclass, 32);
			for(int nsp = ndstart; nsp <= ndend; nsp++)
			{
			   int nc = ncase[nsp-1];
			   int l = a[(ncase[nsp-1]-1)*mdim + (mvar-1)];
			   tclasscat[(l-1)*nclass+(cl[nc-1]-1)] = tclasscat[(l-1)*nclass+(cl[nc-1]-1)] + win[nc-1];
			}
			int nnz = 0;   
			for(int i = 1; i<=lcat; i++)
			{
				double su = 0.0; 
				for(int j = 1; j<=nclass; j++)
				{
					su = su + tclasscat[(i-1)*nclass+(j-1)];
				}
				dn[i-1] = su;
				if(su > 0) nnz = nnz + 1;
			}
			nhit = 0;
			if (nnz > 1) 
			{
			   if ((nclass == 2) && (lcat > ncmax))   
			   {
				/*   catmaxb(pdo, tclasscat, tclasspop, nclass,
	 &                 lcat, nbest, critmax, nhit, dn)*/
				  //catmaxb(pdo, tclasscat, tclasspop, nclass, lcat, *nbest, critmax, nhit, dn);	//TODO 
			   }
			   else
			   {
				  /*call catmax(pdo, tclasscat, tclasspop, nclass, lcat,
	 &                 nbest, critmax, nhit, maxcat, ncmax, ncsplit)*/
				  //catmax(pdo, tclasscat, tclasspop, nclass, lcat, *nbest, critmax, nhit, maxcat, ncmax, ncsplit);	 //TODO
			   }  
			   if (nhit == 1) *msplit = mvar;
	//            else
	//               critmax = -1.0e25
			}
		}
	}
	
	if (*msplit == 0) // One stop criterion: no optimal split var on this node
	{ 
		*jstat = -1;
	}
	
	if (critmax < -1.0e10) // One stop criterion: almost no information gain (not use this criterion)
	{ 
		// *jstat = -1;
	}
	
	*decsplit = critmax; // used for computing Gini importance of variable
	
	// To compute var association
    if (*msplit > 0)
	{ 
        /* ====== The association between vars ========= */
		// To compute association between visual vars, only positive association
		// See: http://www.mathworks.co.uk/help/stats/compactregressiontree.meansurrvarassoc.html 		
		// -- The node probabilities for the optimal split of node i into Left and Right nodes respectively
		double num_sample = double(ndend - ndstart + 1);
		double P_l = (double(*nbest - ndstart + 1)) / num_sample;
		double P_r = (double(ndend - *nbest)) / num_sample;
		double P_lr_min = min(P_l, P_r);
		
		// Update the side indicator of each sample
		for (int nsp = ndstart; nsp <= *nbest; nsp++)
		{
			int idx_samp = a[(nsp-1)*mdim+(*msplit-1)];
			indicator_LR_child[idx_samp-1] = 1; // left child
		}
		for (int nsp = *nbest+1; nsp <= ndend; nsp++)
		{
			int idx_samp = a[(nsp-1)*mdim+(*msplit-1)];
			indicator_LR_child[idx_samp-1] = 2; // left child
		}

		// Compute association between vars
		for (int idx_try = 0; idx_try < mtry; idx_try++)
		{
			int surr_var_idx = surr_best_splits[idx_try*2]; // start from 1
			int best_split_surr_var = surr_best_splits[idx_try*2 + 1];
			if(best_split_surr_var == 0)
			{
				continue;
			}
			/* // A bug found, for different vars, the sample order should be different
			int com_split_left = min(best_split_surr_var, *nbest);
			int com_split_right = max(best_split_surr_var, *nbest);
			// -- The probability that both (optimal) node i and (surrogate) node j send an observation to the Left/Right.
			double P_L_ij = (double(com_split_left - ndstart + 1)) / (double(ndend - ndstart + 1));
			double P_R_ij = (double(ndend - com_split_right)) / (double(ndend - ndstart + 1)); */
			
			double P_L_ij = 0;
			for (int nsp = ndstart; nsp <= best_split_surr_var; nsp++)
			{
				int idx_samp = a[(nsp-1)*mdim+(surr_var_idx-1)];
				if(indicator_LR_child[idx_samp-1] == 1)
				{
					P_L_ij++;
				}
			}
			P_L_ij = P_L_ij / num_sample;
			
			double P_R_ij = 0;
			for (int nsp = best_split_surr_var+1; nsp <= ndend; nsp++)
			{
				int idx_samp = a[(nsp-1)*mdim+(surr_var_idx-1)];
				if(indicator_LR_child[idx_samp-1] == 2)
				{
					P_R_ij++;
				}
			}
			P_R_ij = P_R_ij / num_sample;
			
			// The association
			double lambda_assoc = (P_lr_min - (1 - P_L_ij - P_R_ij)) / P_lr_min;
			//printf("+++++ num_sample=%f, best_split_surr_var=%d, ----> (P_L_ij=%f, P_R_ij=%f)\n", num_sample, best_split_surr_var, P_L_ij, P_R_ij);
			// Only use positive association, no negative, 
			lambda_assoc = max(lambda_assoc, 0);
			
			/* This average is computed by summing positive values of the predictive measure of association 
			   over optimal splits on predictor i and surrogate splits on predictor j and 
			   dividing by the total number of optimal splits on predictor i, including splits 
			   for which the predictive measure of association between predictors i and j is negative.*/
			surr_var_assoc[(surr_var_idx-1)*mdim + (*msplit - 1)] += lambda_assoc;
            num_measured_assoc[(surr_var_idx-1)*mdim + (*msplit - 1)] += 1;
		}
		
		// To compute association between visual and non-visual vars
		if(has_nv_src_mode == 1)
		{
			/* The algorithm based on rank of information
			// -- To rank non-visual source based on their purity increase, start from 1
			for (int idx = 0; idx < num_nv_src; idx++)
			{
				rank_nv_src[idx] = idx + 1;
			}
			// printf("-------- To sort purity increase of nv sources \n");
			R_qsort_I(purity_increase_nv_src, rank_nv_src, 1, num_nv_src);

			// -- The rank number, the higher the more correlated since they are sorted ascendly
			for (int idx_rank = 0; idx_rank < num_nv_src; idx_rank++)
			{
				int idx_nv = rank_nv_src[idx_rank] - 1; // 0, 1, ..., num_nv_src - 1
				surr_var_assoc[(mdim+idx_nv)*mdim + (*msplit-1)] += idx_rank + 1; // add rank order: 1, 2, ..., num_nv_src 
				num_measured_assoc[(mdim+idx_nv)*mdim + (*msplit-1)] += 1;
			} */
			
			/* The algorithm based on information gain */
			for (int idx_nv = 0; idx_nv < num_nv_src; idx_nv++)
			{	
				// Only consider positive cases
				if(purity_increase_nv_src[idx_nv] > 0)
				{
					surr_var_assoc[(mdim+idx_nv)*mdim + (*msplit-1)] += purity_increase_nv_src[idx_nv];
					//printf("-----$$$$$ purity_increase_nv_src[%d] = %f\n", idx_nv, purity_increase_nv_src[idx_nv]);
				}
				
				num_measured_assoc[(mdim+idx_nv)*mdim + (*msplit-1)] += 1; // if moved to if branch???
				//else
				//{
					//printf("-----$$$$$ purity_increase_nv_src[%d] = %f\n", idx_nv, purity_increase_nv_src[idx_nv]);
				//}
			}
		}
    
		/* ====== update information gain (visual part) ====== */
		if(visual_gain_max > 0)
		{
			*visual_info_gain += visual_gain_max;
		}
	}

	return;
}


// ==============================================================
/*
	SUBROUTINE MOVEDATA
	This subroutine is the heart of the buildtree construction. Based on the 
	best split the data in the part of the a matrix corresponding to the 
	current node is moved to the left if it belongs to the left child and 
	right if it belongs to the right child.
*/
void movedata(	int* a, int* ta, int mdim, int nsample, int ndstart, int ndend, int* idmove,       
				int* ncase, int msplit, int* cat, int nbest, int* ndendl)
{
	/* implicit double precision(a-h, o-z)
	   integer a(mdim, nsample), ta(nsample), idmove(nsample),
	   ncase(nsample), cat(mdim), icat(32) */
    
	// Compute idmove = indicator of case nos. going left
	int icat[32];	
    if (cat[msplit-1] == 1)
	{
		for(int nsp = ndstart; nsp <= nbest; nsp++) // to the left child
		{
			int nc = a[(nsp-1)*mdim+(msplit-1)];
			idmove[nc-1] = 1;
		}
		for(int nsp = (nbest+1); nsp <= ndend; nsp++) // to the right child
		{
			int nc = a[(nsp-1)*mdim+(msplit-1)];
			idmove[nc-1] = 0;
		}
		*ndendl = nbest;
    }
    else
	{
		*ndendl = ndstart-1;
		int l = cat[msplit-1];
		myunpack(l, nbest, icat);
		for(int nsp = ndstart; nsp <= ndend; nsp++)
		{
			int nc = ncase[nsp-1];
			int abuf = (nc-1)*mdim+(msplit-1);
			if (icat[abuf-1] == 1) 
			{
			   idmove[nc-1] = 1;
			   *ndendl = *ndendl+1;
			}
			else
			{
			   idmove[nc-1] = 0;
			}
		}
	}
      
	// Shift case. nos. right and left for numerical variables.        
    for(int msh = 1; msh <= mdim; msh++)
	{
        if (cat[msh-1] == 1) 
		{
            int k = ndstart-1;
            for(int n = ndstart; n <= ndend; n++)
			{
               int ih = a[(n-1)*mdim+(msh-1)];
               if (idmove[ih-1] == 1) 
			   {
                  k = k+1;
                  ta[k-1] = a[(n-1)*mdim+(msh-1)];
               }
			}
            for(int n = ndstart; n <= ndend; n++)
			{
               int ih = a[(n-1)*mdim+(msh-1)];
               if (idmove[ih-1]==0) 
			   {			   
                  k = k+1;
                  ta[k-1] = a[(n-1)*mdim+(msh-1)];
               }
			}          
            for(int k = ndstart; k <= ndend; k++)
			{
               a[(k-1)*mdim+(msh-1)] = ta[k-1];
			}
        }        
	}
      
	// Compute case nos. for right and left nodes.
	if (cat[msplit-1] == 1) 
	{
		for(int n = ndstart; n <= ndend; n++)
		{
			ncase[n-1] = a[(n-1)*mdim+(msplit-1)];
		}
	}
	else
	{
		int k = ndstart-1;
		for(int n=ndstart; n <= ndend; n++)
		{
			if (idmove[ncase[n-1]-1]==1) 
			{
				k = k+1;
				ta[k-1] = ncase[n-1];
			}
		}
		for(int n = ndstart; n <= ndend; n++)
		{
			if (idmove[ncase[n-1]-1]==0) 
			{
				k = k+1;
				ta[k-1] = ncase[n-1];
			}
		}
		for(int k = ndstart;k <= ndend; k++)
		{
			ncase[k-1] = ta[k-1];
		}
    }
}
	
//---------------------------------------------------------------------------
void myunpack(int l,int npack,int* icat)
{    
	// npack is a long integer.  The sub. returns icat, an integer of zeroes and
	// ones corresponding to the coefficients in the binary expansion of npack.
  
	// integer icat(32), npack
	for(int j=1; j<=32; j++)
	{
		icat[j-1]=0;
	}
	int n = npack;
	icat[1-1] = n%2;    
	for(int k = 2; k <= l; k++)
	{
		n=(n-icat[(k-1)-1])/2;
		icat[k-1]=n%2;  
	}  
}
  
//---------------------------------------------------------------------------
/*void zeroInt(int* ix, int m1)
{
	//integer ix(m1)
	for(int n=1; n<=m1; n++)
	{
		ix[n-1]=0;
	}
}

//---------------------------------------------------------------------------
void zeroDouble(double* rx, int m1)
{
	// double precision rx(m1)
	for(int n=1; n <= m1; n++)
	{
		rx[n-1]=0.0;   
	}
}*/

//---------------------------------------------------------------------------
void zeroInt2(int* mx, int m1, int m2)
{
	//integer mx(m1,m2)
	for(int i = 1; i <= m1; i++)
	{
		for(int j = 1; j <= m2; j++)
		{
			mx[(j-1)*m1+(i-1)] = 0;
		}
	}
}
  
//---------------------------------------------------------------------------
void zeroDouble2(double* rx, int m1, int m2)
{
	//double precision rx(m1,m2)
	for(int i=1; i<=m1; i++)
	{
		for(int j=1; j<=m2; j++)
		{
			rx[(j-1)*m1+(i-1)]=0.0;   
		}
	}
}

//---------------------------------------------------------------------------
void zeroDouble3(double* rx, int m1, int m2, int m3)
{
	//double precision rx(m1,m2)
	for(int i = 1; i <= m1; i++)
	{
		for(int j = 1; j <= m2; j++)
		{
			for(int k = 1; k <= m3; k++)
				rx[(k-1)*m1*m2+(j-1)*m1+(i-1)] = 0.0;   
		}
	}
}


//---------------------------------------------------------------------------
void rrand(double* xrand)
{
	//*xrand = (((double)randomMT())/((double)MAX_UINT_COKUS));
	double tmp = double(rand())/double(RAND_MAX);
	if(tmp >= 1)
	{
		tmp = 1 - 0.000000001;
	}
	*xrand = tmp;
}

//---------------------------------------------------------------------------
double min(double num1, double num2)
{
	if(num1 > num2)
		return num2;
	else
		return num1;
}

//---------------------------------------------------------------------------
double max(double num1, double num2)
{
	if(num1 < num2)
		return num2;
	else
		return num1;
}
