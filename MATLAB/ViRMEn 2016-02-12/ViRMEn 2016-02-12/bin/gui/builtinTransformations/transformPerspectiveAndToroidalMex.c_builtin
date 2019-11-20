#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize ncols, index;
    double *coord3new, *coord3;
    float aspectRatio;
    int s, p, ndims, xSign, zSign, offset;
    double r, rinv, rnew;
    int dims[3];
    
    // Get aspect ratio of window
    aspectRatio = 1.8;
    
    // viewing parameters
    s = 1;
    p = 1;
    
    // inputs
    ncols = mxGetN(prhs[0]);
    coord3 = mxGetPr(prhs[0]);
    
    // outputs
    ndims = 3;
    dims[0] = 3;
    dims[1] = ncols;
    dims[2] = 2;
    plhs[0] = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
    coord3new = mxGetPr(plhs[0]);
    offset = 3 * ncols;
    
    // coordinate transformations (toroidal + perspective)
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
    
        // toroidal transformation
        // -------------------------
        coord3new[offset+3*index+2] = 1;
        r = sqrt(coord3[3*index]*coord3[3*index]+coord3[3*index+1]*coord3[3*index+1]);
        rnew = 0.4625*atan(coord3[3*index+2]/r)+0.4929;
        rinv = 1;
        if ( rnew < 0 ) {
            rnew = 0;
            coord3new[offset+3*index+2] = 0;
        }
        if ( rnew > 1 ) {
            rnew = 1;
            coord3new[offset+3*index+2] = 0;
        }
        
        coord3new[offset+3*index] = rnew*coord3[3*index]/r;
        coord3new[offset+3*index+1] = rnew*coord3[3*index+1]/r;
    }
    
    return;
}