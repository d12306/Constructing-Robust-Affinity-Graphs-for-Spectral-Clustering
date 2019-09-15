#include "mex_stub.h"


void *mxCalloc( unsigned int n, unsigned int size )
{
	return calloc( (size_t) n, (size_t) size );
}

void mxFree(void *ptr)
{
	delete [] ptr;
}
