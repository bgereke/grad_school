#include "mex.h"
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize *triangles, *sorted;
    double *ord;
    int numTriangles, numTrans;
    int dims[3];
    int i, d;
    long indx;
    
    triangles = mxGetPr(prhs[0]);
    numTriangles = mxGetScalar(prhs[1]);
    numTrans = mxGetScalar(prhs[2]);
    ord = mxGetPr(prhs[3]);
    
    dims[0] = 3;
    dims[1] = numTriangles;
    dims[2] = numTrans;
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    sorted = mxGetPr(plhs[0]);
    
    for (i = 0; i < numTriangles; i++) {
        for (d = 0; d < numTrans; d++) {
            indx = 3*ord[i]-3;
            sorted[3*numTriangles*d+3*i] = triangles[3*numTriangles*d+indx];
            sorted[3*numTriangles*d+3*i+1] = triangles[3*numTriangles*d+indx+1];
            sorted[3*numTriangles*d+3*i+2] = triangles[3*numTriangles*d+indx+2];
        }
    }
    
    return;
}