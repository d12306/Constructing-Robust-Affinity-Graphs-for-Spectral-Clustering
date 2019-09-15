//add xk 12-9-19

#include "startC.h"   
#include "classRF.h"   

/* The details of argues: type and size
 *-------------------Input----------------------
 double x(dimx[0]*nsample), classwt(nclass), cut(nclass), countts[nclass*nts],conswt[nclass]
 *      
 *double impsd, xts, proxts,     
 *
 int dimx[2], cl[nsample], cat[dimx[0]], Options[10], consY[nsample*ncons]
 *
 int ncl, nsample, maxcat, sampsize, strata, ntree, nvar,ipi, nodesize,
 *   nrnodes, print_verbose_tree_progression, ncons
 *   testdat, clts, nts, outclts, labelts
 *   
 *
 *-------------------Output----------------------
 * double prox[dimx[1]*dimx[1]], imprt(dimx[0]), impmat[nsample*dimx[1]],
 *        xbestsplit[nrnodes*ntree], errts[(nclass+1)*ntree], errtr[(ntree+1)*ntree]
 *        
 *         
 * int counttr[nclass*nsample], outcl[nsample], ndbigtree[nrnodes*ntree],
 *     nodestatus[nrnodes*ntree], bestvar[nrnodes*ntree], treemap[nrnodes*ntree*2],
 *     nodeclass[nrnodes*ntree], inbag[dimx[1]],
 *     
 *     
 *     
 *
 *
 *** note: strata may vary with some combination of options       /
 *   格式：定义=调用  cl=y, ncl=nclass, nvar = mtry, cut=cutoff, imprt=impout, impsd=impSD, 
 *               nodeclass=nodepred, clts=yts
 *                   
 */
//-----------------------------------------------------------
/*C++调试入口函数*/
void start()
{
	/*读取空间不变数组/常数*/
	//1.申明数据空间
	int nclass=0, maxcat=0, sampsize=0, strata=0, ntree=0, mtry=0, ipi=0, \
	nodesize=0, nrnodes=0, print_verbose_tree_progression=0,\
	testdat=0, yts=0, nts=0, outclts=0, labelts=0, nsample=0,
	num_nv_src_discrete=0, max_num_class_nv_src_discrete, num_nv_src_continuous=0;     
	double impSD=0.0, xts=0.0, proxts=0.0;     
	int* dimx = new int[2];
	int* Options = new int[10];
	//2.初始化
	zeroInt_(dimx, 2);
	zeroInt_(Options, 10);
	//3.文件中读出参数
	readConstParam(dimx, &nclass, &maxcat, &sampsize, &strata, Options, &ntree, &mtry,&ipi, 
				&nodesize, &impSD, &nrnodes, &testdat, &xts, &yts, &nts,&outclts, &labelts, 
				&proxts,&print_verbose_tree_progression, &nsample,
				&num_nv_src_discrete, &max_num_class_nv_src_discrete, &num_nv_src_continuous);
	

	
	/*读取空间可变数组*/
	//1.申明数据空间
	double* x = new double[dimx[0]*nsample];
	int* y = new int[nsample];
	int* cat = new int[dimx[0]];
	double* classwt = new double[nclass];
	double* cutoff = new double[nclass];
	double* countts = new double[nclass*nts];

	int* Y_nv_src_discrete = new int[nsample * num_nv_src_discrete];
	double* alpha_nv_src_discrete = new double[num_nv_src_discrete];
	double* Y_nv_src_continuous = new double[nsample * num_nv_src_continuous];
	double* alpha_nv_src_continuous = new double[num_nv_src_continuous];

	//2.初始化
	zeroDouble_(x, dimx[0]*nsample);
	zeroInt_(y, nsample);
	zeroInt_(cat, dimx[0]);
	zeroDouble_(classwt, nclass);
	zeroDouble_(cutoff, nclass);
	zeroDouble_(countts, nclass*nts);
   
	//3.文件中读出参数
	readParam(x, y, cat, classwt, cutoff, countts, 
			num_nv_src_discrete, Y_nv_src_discrete, alpha_nv_src_discrete,
			num_nv_src_continuous, Y_nv_src_continuous, alpha_nv_src_continuous, 
            dimx[0], nsample, nclass, nts);
 
	/*输出参数*/
	//1.申明数据空间
	int* outcl = new int[nsample];
	int* counttr = new int[nclass*nsample];
	double* prox = new double[dimx[1]*dimx[1]];
	double* impout = new double[dimx[0]];
	double* impmat = new double[nsample*dimx[1]];
	int* ndbigtree = new int[nrnodes*ntree];
	int* nodestatus = new int[nrnodes*ntree];
	int* bestvar = new int[nrnodes*ntree];
	int* treemap = new int[nrnodes*ntree*2];
	int* nodepred = new int[nrnodes*ntree];
	double* xbestsplit = new double[nrnodes*ntree];
	double* errtr = new double[(ntree+1)*ntree];
	int* inbag = new int[dimx[1]];
	double* errts = new double[(nclass+1)*ntree];
	int num_nv_src = num_nv_src_discrete + num_nv_src_continuous;
	double *meanSurrVarAssoc = new double[dimx[0] * (dimx[0] + num_nv_src)];
	double visual_info_gain = 0;
	double num_flat_region = 0;
	double num_split = 0;
	int* node_size_table = new int[nrnodes*ntree];
	double pseudo_X = 0;

   //2.初始化
   zeroInt_(outcl, nsample);
   zeroInt_(counttr, nclass*nsample);
   zeroDouble_(prox, dimx[1]*dimx[1]);
   zeroDouble_(impout, dimx[0]);
   zeroDouble_(impmat, nsample*dimx[1]);
   zeroInt_(ndbigtree, nrnodes*ntree);
   zeroInt_(nodestatus, nrnodes*ntree);
   zeroInt_(bestvar, nrnodes*ntree);
   zeroInt_(treemap, nrnodes*ntree*2);
   zeroInt_(nodepred, nrnodes*ntree);    
   zeroDouble_(xbestsplit, nrnodes*ntree);
   zeroDouble_(errtr, (ntree+1)*ntree);
   zeroInt_(inbag, dimx[1]);
   zeroDouble_(errts, (nclass+1)*ntree);   
			 
   //调用classRF
   classRF(x, dimx, y, &nclass, cat, &maxcat,
			&sampsize, &strata, Options, &ntree, &mtry,&ipi, 
			classwt, cutoff, &nodesize,outcl, counttr, prox,
			impout, &impSD, impmat, &nrnodes,ndbigtree, nodestatus, 
			bestvar, treemap,nodepred, xbestsplit, errtr,&testdat, 
			&xts, &yts, &nts, countts,&outclts, labelts, 
			&proxts, errts,inbag,print_verbose_tree_progression,
			num_nv_src_discrete, Y_nv_src_discrete, 
			alpha_nv_src_discrete, max_num_class_nv_src_discrete,
			num_nv_src_continuous, Y_nv_src_continuous, 
			alpha_nv_src_continuous, meanSurrVarAssoc,
			&visual_info_gain, &num_flat_region, &num_split,
			node_size_table, &pseudo_X);

	delete [] dimx;
	delete [] Options;
	delete [] x;
	delete [] y;
	delete [] cat;
	delete [] classwt;
	delete [] cutoff;
	delete [] countts;
	delete [] Y_nv_src_discrete;
	delete [] alpha_nv_src_discrete;
	delete [] Y_nv_src_continuous;
	delete [] alpha_nv_src_continuous;
	delete [] outcl;
	delete [] counttr;
	delete [] prox;
	delete [] impout;
	delete [] impmat;
	delete [] ndbigtree;
	delete [] nodestatus;
	delete [] bestvar;
	delete [] treemap;
	delete [] nodepred;
	delete [] xbestsplit;
	delete [] errtr;
	delete [] inbag;
	delete [] errts;
    delete [] meanSurrVarAssoc;
	delete [] node_size_table;
}



/**********************************
*读文件
***********************************/
//文件中读出参数
void readConstParam(int *dimx, int *ncl, int *maxcat,int *sampsize, int *strata, int *Options, 
                    int *ntree, int *nvar,int *ipi, int *nodesize, double *impsd, int *nrnodes,
                    int *testdat, double *xts, int *clts, int *nts, int *outclts, int* labelts, 
					double *proxts, int* print_verbose_tree_progression, int*nsample, 
					int* num_nv_src_discrete, int* max_num_class_nv_src_discrete,
					int* num_nv_src_continuous)
{
   //申明读出数据流
   ifstream intFile("c:\\para_data_CCF.dat", ios::binary );   
   if(!intFile)
   {
	  cout<<"para_data_CCF File open error!\n";
	  return;
   }

	//读数据
	read_int_array(dimx, 2, intFile);           
	read_int_num(ncl, intFile);
	read_int_num(maxcat, intFile);
	read_int_num(sampsize, intFile);
	read_int_num(strata, intFile);
	read_int_array(Options, 10, intFile);
	read_int_num(ntree, intFile);
	read_int_num(nvar, intFile);
	read_int_num(ipi, intFile);
	read_int_num(nodesize, intFile);
	read_db_num(impsd, intFile);
	read_int_num(nrnodes, intFile);
	read_int_num(testdat, intFile);
	read_db_num(xts, intFile);  
	read_int_num(clts, intFile);
	read_int_num(nts, intFile);
	read_int_num(outclts, intFile);
	read_int_num(labelts, intFile);
	read_db_num(proxts, intFile);
	read_int_num(print_verbose_tree_progression, intFile);
	
	read_int_num(num_nv_src_discrete, intFile);
	read_int_num(max_num_class_nv_src_discrete, intFile);
	read_int_num(num_nv_src_continuous, intFile); 
	read_int_num(nsample,intFile);
   //关闭文件
   intFile.close();
}
//----------------------------
//文件中读出参数
void readParam(	double *x, int *cl, int *cat, double *classwt, double *cut,
				double *countts, 
				int num_nv_src_discrete, int *Y_nv_src_discrete, double *alpha_nv_src_discrete,
				int num_nv_src_continuous, double *Y_nv_src_continuous, double *alpha_nv_src_continuous, 
				int dimx0, int nsample, int ncl, int nts)		
{
	//申明读出数据流
	ifstream intFile("c:\\training_data_CCF.dat", ios::binary );  
	if(!intFile)
	{
	  cout<<"training_data_CCF File open error!\n";
	  return;
	}
	//读数据
	read_db_array(x, dimx0 * nsample, intFile);   
	read_int_array(cl, nsample, intFile);          
	read_int_array(cat, dimx0, intFile);
	read_db_array(classwt,ncl, intFile); 
	read_db_array(cut,ncl, intFile); 
	read_db_array(countts,ncl*nts, intFile);
	read_int_array(Y_nv_src_discrete, nsample * num_nv_src_discrete, intFile);
	read_db_array(alpha_nv_src_discrete, num_nv_src_discrete, intFile);
	read_db_array(Y_nv_src_continuous, nsample * num_nv_src_continuous, intFile);
	read_db_array(alpha_nv_src_continuous, num_nv_src_continuous, intFile);   

	//关闭文件
	intFile.close();

	//show data
	//for(i=0; i<200; i++)
	//   cout<<temp[j];
}
//-------------------------------------
//读取double型矩阵
void read_db_matrix(double* matrix, int row, int col, ifstream& intFile)
{  
   for(int i=0; i<col; i++)
	 for(int j=0; j<row; j++)  
        intFile.read((char*)&matrix[i*row+j], sizeof(double));
}
//-------------------------------------
//读取double型数组
void read_db_array(double* _array, int x, ifstream& intFile)
{
   for(int i=0; i<x; i++)
	  intFile.read((char*)&_array[i], sizeof(double));
}
//-------------------------------------
//读取double型数字
void read_db_num(double* num, ifstream& intFile)
{
	intFile.read((char*)num, sizeof(double));
}
//--------------------------------------
//读取int型矩阵
void read_int_matrix(int* matrix, int row, int col, ifstream& intFile)
{
   for(int i=0; i<col; i++)
	 for(int j=0; j<row; j++)
	   intFile.read((char*)&matrix[i*row+j], sizeof(int));
}
//-------------------------------------
//读取int型数组
void read_int_array(int* _array, int x, ifstream& intFile)
{
   for(int i=0; i<x; i++)
	  intFile.read((char*)&_array[i], sizeof(int));
}
//-------------------------------------
//读取int型数字
void read_int_num(int* num, ifstream& intFile)
{
	intFile.read((char*)num, sizeof(int));
}
/**********************************
*写文件
***********************************/
//存储参数至文件
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
		double *alpha_nv_src_continuous, int nsample)
{
	/*存储空间不变数组/常数*/
	//申明写入数据流
	ofstream outFile("c:\\para_data_CCF.dat", ios::binary );
	if(!outFile)
	{
		cout<<"para_data_CCF File open error!\n";
		return;
	}
	//写数据（空间固定）   
	save_int_array(dimx, 2, outFile);  
	save_int_num(ncl, outFile);
	save_int_num(maxcat, outFile);
	save_int_num(sampsize, outFile);
	save_int_num(strata, outFile);    
	save_int_array(Options, 10, outFile);
	save_int_num(ntree, outFile);
	save_int_num(nvar, outFile);
	save_int_num(ipi, outFile);
	save_int_num(nodesize, outFile);
	save_db_num(impsd, outFile);
	save_int_num(nrnodes, outFile);
	save_int_num(testdat, outFile);
	save_db_num(xts, outFile);    
	save_int_num(clts, outFile);
	save_int_num(nts, outFile);
	save_int_num(outclts, outFile); 
	save_int_num(labelts, outFile);
	save_db_num(proxts, outFile);
	save_int_num(print_verbose_tree_progression, outFile);
	save_int_num(&num_nv_src_discrete, outFile);
	save_int_num(&max_num_class_nv_src_discrete, outFile);
	save_int_num(&num_nv_src_continuous, outFile);
	save_int_num(&nsample, outFile);   
	//关闭文件
	outFile.close();
   
	/*存储空间可变数组*/
	//申明写入数据流
	ofstream outFile1("c:\\training_data_CCF.dat", ios::binary );
	if(!outFile1)
	{
		cout<<"training_data_CCF File open error!\n";
		return;
	}
	//写数据（空间不固定） 
	save_db_array(x, dimx[0]* nsample, outFile1); 
	save_int_array(cl, nsample, outFile1);  
	save_int_array(cat, dimx[0], outFile1);
	save_db_array(classwt, *ncl, outFile1);
	save_db_array(cut, *ncl, outFile1); 
	save_db_array(countts, (*ncl)*(*nts), outFile1); 
	save_int_array(Y_nv_src_discrete, nsample * num_nv_src_discrete, outFile1);
	save_db_array(alpha_nv_src_discrete, num_nv_src_discrete, outFile1);
	save_db_array(Y_nv_src_continuous, nsample * num_nv_src_continuous, outFile1);
	save_db_array(alpha_nv_src_continuous, num_nv_src_continuous, outFile1);

	//关闭文件
	outFile1.close();
}


//--------------------------------------
//存储double型矩阵
void save_db_matrix(double* matrix, int row, int col, ofstream& outFile)
{  
   for(int i=0; i<col; i++)
	 for(int j=0; j<row; j++)
        outFile.write((char*)&matrix[i*row+j], sizeof(double));
}
//-------------------------------------
//存储double型数组
void save_db_array(double* _array, int x, ofstream& outFile)
{
   for(int i=0; i<x; i++)
       outFile.write((char*)&_array[i], sizeof(double));
}
//-------------------------------------
//存储double型数字
void save_db_num(double* num, ofstream& outFile)
{
	outFile.write((char*)num, sizeof(double));
}
//--------------------------------------
//存储int型矩阵
void save_int_matrix(int* matrix, int row, int col, ofstream& outFile)
{
   for(int i=0; i<col; i++)
	 for(int j=0; j<row; j++)
        outFile.write((char*)&matrix[i*row+j], sizeof(int));
}
//-------------------------------------
//存储int型数组
void save_int_array(int* _array, int x, ofstream& outFile)
{
   for(int i=0; i<x; i++)
       outFile.write((char*)&_array[i], sizeof(int));
}
//-------------------------------------
//存储int型数字
void save_int_num(int* num, ofstream& outFile)
{
	//outFile << *num;
	outFile.write((char*)num, sizeof(int));
}
//----------------------------------
void zeroInt_(int* a, int num)
{
	for(int i=0; i<num; i++)
         a[i]=0;
}
void zeroDouble_(double* a, int num)
{
	for(int i=0; i<num; i++)
		a[i]=0.0;
}
