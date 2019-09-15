#ifndef NEW_MEX_H
#define NEW_MEX_H

#include <stdlib.h>
 #include <stdio.h> 


#define mexPrintf printf
#define Rprintf printf


void *mxCalloc( unsigned int n, unsigned int size );

void mxFree(void *ptr);

#endif