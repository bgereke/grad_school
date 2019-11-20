/*
 * permMC.h
 *
 * Code generation for function 'permMC'
 *
 */

#ifndef __PERMMC_H__
#define __PERMMC_H__

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "permMC_types.h"

/* Function Declarations */
extern void permMC(const emxArray_creal_T *fCS, const emxArray_real_T *delta,
                   double sigma, const emxArray_real32_T *start, float numperms,
                   emxArray_creal_T *pMC);

#endif

/* End of code generation (permMC.h) */
