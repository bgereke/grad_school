#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *coords, *colors, *cnt, *lettGrid, *lett, *lettSize, *lettPos, *lettNum, *color;
    mwSize indx, sz, row, lettNdx, coordNdx;
       
    coords = mxGetPr(prhs[0]);
    colors = mxGetPr(prhs[1]);
    cnt = mxGetPr(prhs[2]);
    lettGrid = mxGetPr(prhs[3]);
    lett = mxGetPr(prhs[4]);
    lettSize = mxGetPr(prhs[5]);
    lettPos = mxGetPr(prhs[6]);
    lettNum = mxGetPr(prhs[7]);
    color = mxGetPr(prhs[8]);
    
    sz = mxGetN(prhs[4]);
    
    for ( indx = 0; indx < sz; indx++ ) {
        lettNdx = lett[indx];
        coordNdx = indx+cnt[0];
        for ( row = 0; row < 4; row++ ) {
            coords[4*coordNdx+row] = lettGrid[4*(lettNdx-1)+row];
        }
        coords[4*coordNdx] = (coords[4*coordNdx]+lettNum[0]-1)*lettSize[0]+lettPos[0];
        coords[4*coordNdx+1] = coords[4*coordNdx+1]*lettSize[0]+lettPos[1];
        coords[4*coordNdx+2] = (coords[4*coordNdx+2]+lettNum[0]-1)*lettSize[0]+lettPos[0];
        coords[4*coordNdx+3] = coords[4*coordNdx+3]*lettSize[0]+lettPos[1];
        colors[6*coordNdx] = color[0];
        colors[6*coordNdx+1] = color[1];
        colors[6*coordNdx+2] = color[2];
        colors[6*coordNdx+3] = color[0];
        colors[6*coordNdx+4] = color[1];
        colors[6*coordNdx+5] = color[2];
    }
    
    cnt[0] = cnt[0] + sz;
    return;
}