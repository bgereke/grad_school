#include <mex.h>
#include <math.h>
#include <matrix.h>

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((Y) < (X) ? (X) : (Y))

/*mex getBestSpotStubLearner.c to compile */

/* A function that computes the optimum covariate and corresponding threshold
 * for a stub learner. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* start by reading the arguments */
    /* replaces this code:
            cumsums = cumsum(r(this.cache.ranks),1);
            %Convert to correlation
            scaledccs = abs(bsxfun(@times,cumsums,this.cache.ips).*this.cache.mask);
            [~,bestspot] = max(scaledccs(:));
       with bestspot = getBestSpotStubLearner(r,this.cache.rank,this.cache.ips,this.cache.mask);
    */
    
    int ncols, nrows, i, j, bestidx, offset, ipslarge;
    double bestres, csum, potentialres;
    double *r;
    int *rank;
    double *ips;
    mxLogical *mask;
    
    r = mxGetPr(prhs[0]);
    rank = (int *) mxGetPr(prhs[1]);
    ncols = mxGetN(prhs[1]);
    nrows = mxGetM(prhs[1]);
    ips = mxGetPr(prhs[2]);
    ipslarge = mxGetN(prhs[2]) > 1;
    mask = (mxLogical *) mxGetPr(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix( ncols, 1, mxREAL);
    double *bestidxs = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(ncols, 1, mxREAL);
    double *bestress = mxGetPr(plhs[1]);
    
    for( i = 0; i < ncols; i++)
    {
        csum = 0;
        bestidx = 0;
        bestres = 0;
        offset = i*nrows;
        for(j = 0; j < nrows; j++)
        {
            csum += r[rank[j+offset]-1];
            potentialres = mask[j+offset] ? (ipslarge ? csum * ips[j + offset] : csum * ips[j]) : 0;
            potentialres = potentialres > 0 ? potentialres : -potentialres;
            if(potentialres > bestres)
            {
                bestidx = j;
                bestres = potentialres;
            }
            
        }
        bestidxs[i] = (double) bestidx;
        bestress[i] = bestres;
    }
}
