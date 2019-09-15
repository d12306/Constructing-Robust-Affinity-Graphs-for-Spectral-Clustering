#ifndef STARTC_H
#define STARTC_H

#include <iostream>
#include <fstream>
using namespace std;

void start(); 
void readConstParam(int *dimx, int *ncl, int *maxcat,int *sampsize, int *strata, int *Options, 
                    int *ntree, int *nvar,int *ipi, int *nodesize, double *impsd, int *nrnodes,
                    int *testdat, double *xts, int *clts, int *nts, int *outclts, int* labelts, 
					double *proxts, int* print_verbose_tree_progression, int*nsample, 
					int* num_nv_src_discrete, int* max_num_class_nv_src_discrete,
					int* num_nv_src_continuous);
void readParam(	double *x, int *cl, int *cat, double *classwt, double *cut,
				double *countts, 
				int num_nv_src_discrete, int *Y_nv_src_discrete, double *alpha_nv_src_discrete,
				int num_nv_src_continuous, double *Y_nv_src_continuous, double *alpha_nv_src_continuous, 
				int dimx0, int nsample, int ncl, int nts);
void read_db_matrix(double* matrix, int row, int col, ifstream& intFile);
void read_db_array(double* array, int x, ifstream& intFile);
void read_db_num(double* num, ifstream& intFile);
void read_int_matrix(int* matrix, int row, int col, ifstream& intFile);
void read_int_array(int* array, int x, ifstream& intFile);
void read_int_num(int* num, ifstream& intFile);
void saveParam(double *x, int *dimx, int *cl, int *ncl, int *cat, int *maxcat,
        int *sampsize, int *strata, int *Options, int *ntree, int *nvar,
        int *ipi, double *classwt, double *cut, int *nodesize,
        int *outcl, int *counttr, double *prox,
        double *imprt, double *impsd, double *impmat, int *nrnodes,
        int *ndbigtree, int *nodestatus, int *bestvar, int *treemap,
        int *nodeclass, double *xbestsplit, double *errtr,
        int *testdat, double *xts, int *clts, int *nts, double *countts,
        int *outclts, int* labelts, double *proxts, double *errts, 
        int *inbag, int* print_verbose_tree_progression,
		int num_nv_src_discrete, int *Y_nv_src_discrete, 
		double *alpha_nv_src_discrete, int max_num_class_nv_src_discrete,
		int num_nv_src_continuous, double *Y_nv_src_continuous, 
		double *alpha_nv_src_continuous, int nsample);
void save_db_matrix(double* matrix, int row, int col, ofstream& outFile);
void save_db_array(double* array, int x, ofstream& outFile);
void save_db_num(double* num, ofstream& outFile);
void save_int_matrix(int* matrix, int row, int col, ofstream& outFile);
void save_int_array(int* array, int x, ofstream& outFile);
void save_int_num(int* num, ofstream& outFile);
void zeroInt_(int* a, int num);
void zeroDouble_(double* a, int num);

#endif