/*
 * _coder_permMC_api.h
 *
 * Code generation for function 'permMC'
 *
 */

#ifndef ___CODER_PERMMC_API_H__
#define ___CODER_PERMMC_API_H__
/* Include files */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef struct_emxArray_creal_T
#define struct_emxArray_creal_T
struct emxArray_creal_T
{
    creal_T *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_creal_T*/
#ifndef typedef_emxArray_creal_T
#define typedef_emxArray_creal_T
typedef struct emxArray_creal_T emxArray_creal_T;
#endif /*typedef_emxArray_creal_T*/
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T
{
    double *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_real_T*/
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /*typedef_emxArray_real_T*/
#ifndef struct_emxArray_real32_T
#define struct_emxArray_real32_T
struct emxArray_real32_T
{
    float *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_real32_T*/
#ifndef typedef_emxArray_real32_T
#define typedef_emxArray_real32_T
typedef struct emxArray_real32_T emxArray_real32_T;
#endif /*typedef_emxArray_real32_T*/

/* Function Declarations */
extern void permMC_initialize(emlrtContext *aContext);
extern void permMC_terminate(void);
extern void permMC_atexit(void);
extern void permMC_api(const mxArray *prhs[5], const mxArray *plhs[1]);
extern void permMC(emxArray_creal_T *fCS, emxArray_real_T *delta, double sigma, emxArray_real32_T *start, float numperms, emxArray_creal_T *pMC);
extern void permMC_xil_terminate(void);

#endif
/* End of code generation (_coder_permMC_api.h) */
