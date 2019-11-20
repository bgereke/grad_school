#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *out, *in;
    mwSize numRows, numColumns;
    double currentColor, hasChanged, isVacancy;
    int row, column, indx;
    
    in = mxGetPr(prhs[0]);
    
    numRows = mxGetM(prhs[0]);
    numColumns = mxGetN(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(numRows, numColumns, mxREAL);
    out = mxGetPr(plhs[0]);
    
    currentColor = 0;
    isVacancy = 1;
    while (isVacancy == 1) {
        currentColor++;
        isVacancy = 0;
        for (indx = 0; indx < numRows*numColumns; indx++) {
            if (out[indx]==0 && isVacancy==0) {
                out[indx] = currentColor;
                isVacancy = 1;
            }
        }
        hasChanged = 1;
        while (hasChanged == 1) {
            hasChanged = 0;
            for (row = 0; row < numRows; row++) {
                for (column = 0; column < numColumns; column++) {
                    indx = row+numRows*column;
                    if (out[indx] == currentColor) {
                        if (row > 0 && in[indx] == in[(row-1)+numRows*column] && out[(row-1)+numRows*column]==0) {
                            out[(row-1)+numRows*column] = currentColor;
                            hasChanged = 1;
                        }
                        if (row < numRows-1 && in[indx] == in[(row+1)+numRows*column] && out[(row+1)+numRows*column]==0) {
                            out[(row+1)+numRows*column] = currentColor;
                            hasChanged = 1;
                        }
                        if (column > 0 && in[indx] == in[row+numRows*(column-1)] && out[row+numRows*(column-1)]==0) {
                            out[row+numRows*(column-1)] = currentColor;
                            hasChanged = 1;
                        }
                        if (column < numColumns-1 && in[indx] == in[row+numRows*(column+1)] && out[row+numRows*(column+1)]==0) {
                            out[row+numRows*(column+1)] = currentColor;
                            hasChanged = 1;
                        }
                    }
                }
            }
        }
    }
    
    return;
}