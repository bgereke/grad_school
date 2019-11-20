#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *coord3, *coord3new, *pos, *z;
    mwSize ncols, index, row;
    double c, s;
    double temp;
    
    ncols = mxGetN(prhs[0]);
    
    coord3 = mxGetPr(prhs[0]);
    pos = mxGetPr(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(3,ncols,mxREAL);
    coord3new = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1,ncols,mxREAL);
    z = mxGetPr(plhs[1]);
    
    for ( index = 0; index < ncols; index++ ) {
        for ( row = 0; row < 3; row++ ) {
            coord3new[3*index+row] = coord3[3*index+row]-pos[row];
        }
        z[index] = sqrt(coord3new[3*index]*coord3new[3*index] + coord3new[3*index+1]*coord3new[3*index+1] + coord3new[3*index+2]*coord3new[3*index+2]);
    }
    
    if (pos[3] != 0) {
        c = cos(-pos[3]);
        s = sin(-pos[3]);
        for ( index = 0; index < ncols; index++ ) {
            temp = c*coord3new[3*index] - s*coord3new[3*index+1];
            coord3new[3*index+1] = s*coord3new[3*index] + c*coord3new[3*index+1];
            coord3new[3*index] = temp;
        }
    }
    
    return;
}