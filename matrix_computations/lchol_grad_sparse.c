
#include "mex.h"
#include "string.h"
#include "math.h"

/* TODO: THIS IMPLEMENTATION IS QUICK AND DIRTY, NOT OPTIMIZED NOR MADE EASY TO USE. */

int compare
(
 const void * a, 
 const void * b
)
{
  return (int)*(double*)a - (int)*(double*)b;
}

/*
 * Cholesky decomposition and its derivative for a sparse matrix K.
 * Decomposes the symmetric positive definite sparse matrix into a
 * lower triangular matrix L and its derivatives dL.
 *
 * ind: The indices of the symbolic structure of L and dL. First index
 * is zero.
 *
 * k: The values of matrix K for the indices.
 *
 * dk: The values of the derivative of K for the indices. This may
 * contain several derivatives.
 *
 * length: The number of indices.
 *
 * d: The number of derivative vectors in dk.
 *
 * n: The dimensionality of the matrix K. 
 */
int cholesky
(
 double *array,
 double *cols,
 int length,
 int D,
 int N
)
{
  int i, j, ind, jnd, ind_pivot, d, n;
  double *item, key;
  double *p, *k, *dk;
  p = array;
  k = array + length;
  dk = array + 2*length;

  ind = 0;
  for (n=1; n<=N; n++) {

    /* Check positive definiteness */
    if (k[ind] <= 0.0)
      return -2;

    /* Define pivot */
    k[ind] = sqrt(k[ind]);
    for (d=0; d<D; d++) {
      dk[ind+d*length] = 0.5 * dk[ind+d*length] / k[ind];
    }
    ind_pivot = ind;

    /* Adjust lead column */
    for (ind = ind_pivot+1; ind < length && p[ind] <= n*N; ind++) {
      k[ind] = k[ind] / k[ind_pivot];
      for (d=0; d<D; d++) {
	dk[ind+d*length] = (dk[ind+d*length] - k[ind]*dk[ind_pivot+d*length]) / k[ind_pivot];
      }
    }

    /* Row operations */
    item = p+ind;
    for (ind = ind_pivot+1; ind < length && p[ind] <= n*N; ind++) {
      /* current row to column index */
      i = ((int)(p[ind]-1) % N) + 1; 
      for (jnd = ind; jnd < length && p[jnd] <= n*N; jnd++) {
	/* current row to row index */
	j = ((int)(p[jnd]-1) % N) + 1; 
	/* find (j,i) from p using binary search */
	key = (i-1)*N+j;
	/*item = (double*)bsearch(&key, item+1, length-(item-p)-1, sizeof(double), compare);*/
	item = (double*)bsearch(&key, p+(int)cols[i-1], cols[i]-cols[i-1], sizeof(double), compare);	
	/* item = (double*)bsearch(&key, p+jnd+1, length-jnd-1, sizeof(double), compare);*/
	if (item==NULL)
	  return -1;
	/*item = p+length-1;*/
	
	*(item+length) = *(item+length) - k[ind]*k[jnd];
	for (d=0; d<D; d++) {
	  *(item+(d+2)*length) = *(item+(d+2)*length) - dk[ind+d*length]*k[jnd] - dk[jnd+d*length]*k[ind];
	}
      }
    }
    
  }


  return 0;
  

  /*
  for (i=0; i<length; i++) {
    *(k+i) = *(k+i) + 1;
    for (j=0; j<d; j++) {
      *(dk+i+j*length) = *(dk+i+j*length) * (j+2);
    }
  }
  */

}

void mexFunction
(
 int nout, 
 mxArray *pout[],
 int nin, 
 const mxArray *pin[]
 )
{
  mwSize mrows, ncols;
  double *x, *y, *cols;
  int *p;
  int N;

  if (nin != 3)
    mexErrMsgTxt("Three inputs required.");
  if (nout != 1)
    mexErrMsgTxt("One output required.");

  mrows = mxGetM(pin[0]);
  ncols = mxGetN(pin[0]);

  /*pout[0] = pin[0];*/
  pout[0] = mxCreateDoubleMatrix(mrows, ncols, mxREAL);

  y = mxGetPr(pout[0]);
  x = mxGetPr(pin[0]);
  cols = mxGetPr(pin[1]);

  memcpy(y, x, mrows*ncols*sizeof(double));

  /*
  p = malloc(mrows*sizeof(int));
  for (int m=0; m<mrows; m++) {
    p[m] = (int)x[m] - 1;
  }
  */

  N = (int)mxGetScalar(pin[2]);
  int retval = cholesky(y, cols, mrows, ncols-2, N);
  /*
  int retval = cholesky(y, y+mrows, y+(2*mrows), mrows, ncols-2, N);
  */
  
  /*free(p);*/

  if (retval == -1) 
    mexErrMsgTxt("Symbolic structure not consistent.");

  if (retval == -2)
    mexErrMsgTxt("Matrix not positive definite.");

  return;
}
