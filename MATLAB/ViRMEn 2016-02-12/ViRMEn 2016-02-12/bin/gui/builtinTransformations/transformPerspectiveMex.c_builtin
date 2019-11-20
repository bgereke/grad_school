#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize ncols, index;
    double *coord3new, *coord3;
    float aspectRatio;
    int s, p, ndims, xSign, zSign;
    int dims[2];
    
    // Get aspect ratio of window
    aspectRatio = 1.8;
    
    // viewing parameters
    s = 1;
    p = 1;
    
    // inputs
    ncols = mxGetN(prhs[0]);
    coord3 = mxGetPr(prhs[0]);
    
    // outputs
    ndims = 2;
    dims[0] = 3;
    dims[1] = ncols;
    plhs[0] = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
    coord3new = mxGetPr(plhs[0]);
    
    for (index = 0; index < ncols; index++) {
        
        // perspective transformation
        // -----------------------------
        coord3new[3*index] = s * coord3[3*index] / coord3[3*index+1];      // monitor x value
        coord3new[3*index+1] = s * coord3[3*index+2] / coord3[3*index+1];  // monitor y value 
        // check if point is visible and if not, clip the point
        if (coord3[3*index+1] <= 0) {
            coord3new[3*index+2] = 0; // invisible (clipping needed)
            if (coord3[3*index] > 0) {
                xSign = 1;
            }
            else {
                xSign = -1;
            }
            if (coord3[3*index+2] > 0) {
                zSign = 1;
            }
            else {
                zSign = -1;
            }
            // here is the clipping
            coord3new[3*index] = p * xSign * aspectRatio;
            coord3new[3*index+1] = p * zSign;
        }
        else {
            coord3new[3*index+2] = 1; // visible (no clipping needed)
        }    
    }
    
    return;
}