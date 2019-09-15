#ifndef CLASSRF_H
#define CLASSRF_H

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
		double * visual_info_gain, double *num_flat_region, double *num_split,
		int *node_size_table, double *pseudo_X);
		
void classForest(int *mdim, int *ntest, int *nclass, int *maxcat,
        int *nrnodes, int *ntree, double *x, double *xbestsplit,
        double *pid, double *cutoff, double *countts, int *treemap,
        int *nodestatus, int *cat, int *nodeclass, int *jts,
        int *jet, int *bestvar, int *node, int *treeSize,
        int *keepPred, int *prox, double *proxMat, int *nodes);
		
void oob(int nsample, int nclass, int *jin, int *cl, int *jtr, int *jerr,
        int *counttr, int *out, double *errtr, int *jest,
        double *cutoff);
		
void TestSetError(double *countts, int *jts, int *clts, int *jet, int ntest,
        int nclass, int nvote, double *errts,
        int labelts, int *nclts, double *cutoff);



#endif