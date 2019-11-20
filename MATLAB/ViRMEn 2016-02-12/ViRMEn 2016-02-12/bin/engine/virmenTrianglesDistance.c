#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize *triangles;
    double *dist, *minDist, *output;
    mwSize numTriangles;
    double *ord;
    int i;
    
    dist = mxGetPr(prhs[0]);
    triangles = (long *)mxGetData(prhs[1]);
    
    numTriangles = mxGetN(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(1, numTriangles, mxREAL);
    minDist = mxGetPr(plhs[0]);
    
    for (i = 0; i < numTriangles; i++) {
        if (dist[triangles[3*i]] < dist[triangles[3*i+1]]) {
            minDist[i] = dist[triangles[3*i]];
        }
        else {
            minDist[i] = dist[triangles[3*i+1]];
        }
        if (dist[triangles[3*i+2]] < minDist[i]) {
            minDist[i] = dist[triangles[3*i+2]];
        }
    }
    
    return;
}