#include "mex.h"
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *vertexArray, *nDim, *nVert;
    bool *isVisible;
    mwSize *tria, nTria, index, d, nCoord;
    mwSize *newTria;
    int dims[3];
    
    nTria = mxGetN(prhs[0]);
    
    tria = mxGetPr(prhs[0]);
    vertexArray = mxGetPr(prhs[1]);
    nDim = mxGetPr(prhs[2]);
    nVert = mxGetPr(prhs[3]);
    isVisible = mxGetPr(prhs[4]);
    
    dims[0] = 3;
    dims[1] = nTria;
    dims[2] = nDim[0];
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    newTria = mxGetPr(plhs[0]);
    
    nCoord = 3*nVert[0];
    
    for ( d = 0; d < nDim[0]; d++ ) {
        for ( index = 0; index < nTria; index++ ) {
            newTria[3*nTria*d+3*index] = 0;
            newTria[3*nTria*d+3*index+1] = 0;
            newTria[3*nTria*d+3*index+2] = 0;
            if (isVisible[index]==1 && (vertexArray[nCoord*d+3*tria[3*index]+2]==1 || vertexArray[nCoord*d+3*tria[3*index+1]+2]==1 || vertexArray[nCoord*d+3*tria[3*index+2]+2]==1)) {
                newTria[3*nTria*d+3*index] = tria[3*index];
                newTria[3*nTria*d+3*index+1] = tria[3*index+1];
                newTria[3*nTria*d+3*index+2] = tria[3*index+2];
            }
        }
    }
    
    return;
}