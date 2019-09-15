#ifndef BUILD_TREE_H
#define BUILD_TREE_H


//-----------------
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

void buildtree(int* a, int* b, int* cl, int* cat, int* maxcat, int mdim, int nsample,                          \
        int nclass, int* treemap, int* bestvar, int* bestsplit, int* bestsplitnext, double* tgini,                 \
        int *nodestatus,int *nodepop, int *nodestart, double *classpop, double *tclasspop,                          \
        double* tclasscat,int* ta,int nrnodes, int * idmove, int  ndsize, int * ncase, int  mtry, int *iv,       \
        int *nodeclass, int* ndbigtree, double *win, double *wr, double *wl, int mred, int nuse, int *mind,       \
        int num_ori_samples,
		double *win_nv_src,
		int num_nv_src_discrete,
		int *Y_nv_src_discrete, 
		double* alpha_nv_src_discrete,
		int max_num_class_nv_src_discrete,			
		double* tclasspop_nv_src_discrete, double* classpop_nv_src_discrete_all_nodes,
		double* parent_crit_nv_src_discrete,
		double* pno_nv_src_discrete, double* pdo_nv_src_discrete, 
		double* rln_nv_src_discrete, double* rld_nv_src_discrete,
		double* rrn_nv_src_discrete, double* rrd_nv_src_discrete,
		double* wl_nv_src_discrete, double* wr_nv_src_discrete,
		int num_nv_src_continuous,
		double *Y_nv_src_continuous, 
		double* alpha_nv_src_continuous,
		double* parent_crit_nv_src_continuous,
		double* sum_nv_src_continuous_all_nodes, 
		double* sum_nv_src_continuous, 
		double* suml_nv_src_continuous, double* sumr_nv_src_continuous,
        double* surr_var_assoc, int* num_measured_assoc, int* indicator_LR_child,
		double* ratio_missing_nv_src,
		double* visual_info_gain, double* num_flat_region, double* num_split,
		int* node_size_table);


void findbestsplit( int* a, int* b, int* cl, int mdim, int nsample, int nclass, int* cat,
        int* maxcat, int ndstart, int ndend, double* tclasspop, double* tclasscat, int* msplit,
        double* decsplit, int* nbest, int* ncase, int* jstat, int mtry, double* win, double* wr, double* wl,
        int mred, int* mind,\
        int num_ori_samples,
        double crit_root_visual, double* crit_root_nv_src_discrete, double* crit_root_nv_src_continuous,
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
        double *Y_nv_src_continuous,
        double* alpha_nv_src_continuous,
        double* parent_crit_nv_src_continuous,
        double* sum_nv_src_continuous_all_nodes, 
        double* sum_nv_src_continuous, 
        double* suml_nv_src_continuous, double* sumr_nv_src_continuous,
        double* surr_var_assoc, int* num_measured_assoc, int* indicator_LR_child,
        double* purity_increase_nv_src, int* rank_nv_src, int* surr_best_splits,
		double* visual_info_gain, double* num_flat_region, double* num_split);

void movedata(int* a,int* ta, int mdim,int nsample,int ndstart,int ndend,int* idmove,    \
			  int* ncase,int msplit, int* cat,int nbest,int* ndendl);

void myunpack(int l,int npack,int* icat);
// void zeroInt(int* ix,int m1);
// void zeroDouble(double* rx, int m1);
void zeroInt2(int* mx,int m1, int m2);
void zeroDouble2(double* rx, int m1, int m2);
void zeroDouble3(double* rx, int m1, int m2, int m3);
void rrand(double* xrand);
double min(double num1, double num2);
double max(double num1, double num2);
#endif