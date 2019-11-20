#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *out, *in, *arr1, *arr2;
    mwSize numRows, numColumns, arrLength;
    int i, j;
     
    in = mxGetPr(prhs[0]);
    arr1 = mxGetPr(prhs[1]);
    arr2 = mxGetPr(prhs[2]);
    
    numRows = mxGetM(prhs[0]);
    numColumns = mxGetN(prhs[0]);
    arrLength = mxGetM(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(numRows, numColumns, mxREAL);
    out = mxGetPr(plhs[0]);
    
    for (i=0; i<numRows*numColumns; i++) {
        out[i] = in[i];
        for (j=0; j<arrLength; j++) {
            if (in[i]==arr1[j]) {
                out[i] = arr2[j];
            }
        }
    }

    return;
}