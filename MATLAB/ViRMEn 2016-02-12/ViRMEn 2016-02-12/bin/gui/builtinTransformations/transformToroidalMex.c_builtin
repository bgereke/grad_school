#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize ncols, index;
    double *coord3new, *coord3;
    double r, rinv, rnew;
    
    ncols = mxGetN(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(3,ncols,mxREAL);

    coord3new = mxGetPr(plhs[0]);
    coord3 = mxGetPr(prhs[0]);
    
    for ( index = 0; index < ncols; index++ ) {
        coord3new[3*index+2] = 1;
        r = sqrt(coord3[3*index]*coord3[3*index]+coord3[3*index+1]*coord3[3*index+1]);
        rnew = 0.4625*atan(coord3[3*index+2]/r)+0.4929;
        rinv = 1;
        if ( rnew < 0 ) {
            rnew = 0;
            coord3new[3*index+2] = 0;
        }
        if ( rnew > 1 ) {
            rnew = 1;
            coord3new[3*index+2] = 0;
        }
        
        coord3new[3*index] = rnew*coord3[3*index]/r;
        coord3new[3*index+1] = rnew*coord3[3*index+1]/r;
    }
    return;
}
