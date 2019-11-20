#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *vert, *znew, *vertNew;
    mwSize nDim, nVert;
    int dims[3];
    int i, d;
    
    vert = mxGetPr(prhs[0]);
    znew = mxGetPr(prhs[1]);
    nVert = mxGetScalar(prhs[2]);
    nDim = mxGetScalar(prhs[3]);
    
    dims[0] = 3;
    dims[1] = nVert;
    dims[2] = nDim;
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    vertNew = mxGetPr(plhs[0]);
    
    for (i = 0; i < nVert; i++) {
        for (d = 0; d < nDim; d++) {
            vertNew[3*nVert*d+3*i] = vert[3*nVert*d+3*i];
            vertNew[3*nVert*d+3*i+1] = vert[3*nVert*d+3*i+1];
            vertNew[3*nVert*d+3*i+2] = znew[i];
        }
    }
    
    return;
}