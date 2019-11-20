/*
 * permMC.c
 *
 * Code generation for function 'permMC'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "permMC.h"
#include "permMC_emxutil.h"

/* Function Definitions */
void permMC(const emxArray_creal_T *fCS, const emxArray_real_T *delta, double
            sigma, const emxArray_real32_T *start, float numperms,
            emxArray_creal_T *pMC)
{
  int fCS_idx_0;
  int i0;
  int loop_ub;
  int p;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  emxArray_real_T *c_y;
  emxArray_creal_T *b_fCS;
  emxArray_creal_T *d_y;
  int i1;
  unsigned int uv0[2];
  int i2;
  double c;
  int cr;
  int i3;
  int br;
  int ib;
  unsigned int unnamed_idx_1;
  int ic;
  double re;
  double im;
  double fCS_im;
  fCS_idx_0 = fCS->size[0];
  i0 = pMC->size[0] * pMC->size[1] * pMC->size[2];
  pMC->size[0] = fCS_idx_0;
  emxEnsureCapacity((emxArray__common *)pMC, i0, (int)sizeof(creal_T));
  fCS_idx_0 = delta->size[1];
  i0 = pMC->size[0] * pMC->size[1] * pMC->size[2];
  pMC->size[1] = fCS_idx_0;
  emxEnsureCapacity((emxArray__common *)pMC, i0, (int)sizeof(creal_T));
  i0 = pMC->size[0] * pMC->size[1] * pMC->size[2];
  pMC->size[2] = (int)numperms;
  emxEnsureCapacity((emxArray__common *)pMC, i0, (int)sizeof(creal_T));
  loop_ub = fCS->size[0] * delta->size[1] * (int)numperms;
  for (i0 = 0; i0 < loop_ub; i0++) {
    pMC->data[i0].re = 0.0;
    pMC->data[i0].im = 0.0;
  }

  p = 0;
  emxInit_real_T(&y, 2);
  emxInit_real_T(&b_y, 2);
  emxInit_real_T(&c_y, 2);
  emxInit_creal_T(&b_fCS, 2);
  emxInit_creal_T(&d_y, 2);
  while (p <= (int)numperms - 1) {
    if ((double)start->data[(int)(1.0F + (float)p) - 1] > fCS->size[1]) {
      i0 = 0;
      i1 = 0;
    } else {
      i0 = (int)start->data[(int)(1.0F + (float)p) - 1] - 1;
      i1 = fCS->size[1];
    }

    if (1.0F > start->data[(int)(1.0F + (float)p) - 1] - 1.0F) {
      loop_ub = -1;
    } else {
      loop_ub = (int)(start->data[(int)(1.0F + (float)p) - 1] - 1.0F) - 1;
    }

    for (i2 = 0; i2 < 2; i2++) {
      uv0[i2] = (unsigned int)fCS->size[i2];
    }

    c = sigma * sigma;
    i2 = y->size[0] * y->size[1];
    y->size[0] = delta->size[0];
    y->size[1] = delta->size[1];
    emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
    fCS_idx_0 = delta->size[0] * delta->size[1];
    for (i2 = 0; i2 < fCS_idx_0; i2++) {
      y->data[i2] = -0.5 * delta->data[i2] * delta->data[i2] / c;
    }

    i2 = y->size[0] * y->size[1];
    for (fCS_idx_0 = 0; fCS_idx_0 < i2; fCS_idx_0++) {
      y->data[fCS_idx_0] = exp(y->data[fCS_idx_0]);
    }

    c = sigma * sigma;
    i2 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = delta->size[0];
    b_y->size[1] = delta->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i2, (int)sizeof(double));
    fCS_idx_0 = delta->size[0] * delta->size[1];
    for (i2 = 0; i2 < fCS_idx_0; i2++) {
      b_y->data[i2] = -0.5 * delta->data[i2] * delta->data[i2] / c;
    }

    i2 = b_y->size[0] * b_y->size[1];
    for (fCS_idx_0 = 0; fCS_idx_0 < i2; fCS_idx_0++) {
      b_y->data[fCS_idx_0] = exp(b_y->data[fCS_idx_0]);
    }

    if (((int)uv0[1] == 1) || (b_y->size[0] == 1)) {
      i2 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = (int)uv0[0];
      c_y->size[1] = b_y->size[1];
      emxEnsureCapacity((emxArray__common *)c_y, i2, (int)sizeof(double));
      fCS_idx_0 = (int)uv0[0];
      for (i2 = 0; i2 < fCS_idx_0; i2++) {
        cr = b_y->size[1];
        for (i3 = 0; i3 < cr; i3++) {
          c_y->data[i2 + c_y->size[0] * i3] = 0.0;
          br = (int)uv0[1];
          for (ib = 0; ib < br; ib++) {
            c_y->data[i2 + c_y->size[0] * i3] += b_y->data[ib + b_y->size[0] *
              i3];
          }
        }
      }
    } else {
      unnamed_idx_1 = (unsigned int)b_y->size[1];
      i2 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = (int)uv0[0];
      emxEnsureCapacity((emxArray__common *)c_y, i2, (int)sizeof(double));
      i2 = c_y->size[0] * c_y->size[1];
      c_y->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)c_y, i2, (int)sizeof(double));
      fCS_idx_0 = (int)uv0[0] * (int)unnamed_idx_1;
      for (i2 = 0; i2 < fCS_idx_0; i2++) {
        c_y->data[i2] = 0.0;
      }

      if (((int)uv0[0] == 0) || (b_y->size[1] == 0)) {
      } else {
        fCS_idx_0 = (int)uv0[0] * (b_y->size[1] - 1);
        cr = 0;
        while (((int)uv0[0] > 0) && (cr <= fCS_idx_0)) {
          i2 = cr + (int)uv0[0];
          for (ic = cr; ic + 1 <= i2; ic++) {
            c_y->data[ic] = 0.0;
          }

          cr += (int)uv0[0];
        }

        br = 0;
        cr = 0;
        while (((int)uv0[0] > 0) && (cr <= fCS_idx_0)) {
          i2 = br + (int)uv0[1];
          for (ib = br; ib + 1 <= i2; ib++) {
            if (b_y->data[ib] != 0.0) {
              i3 = cr + (int)uv0[0];
              for (ic = cr; ic + 1 <= i3; ic++) {
                c_y->data[ic] += b_y->data[ib];
              }
            }
          }

          br += (int)uv0[1];
          cr += (int)uv0[0];
        }
      }
    }

    fCS_idx_0 = fCS->size[0];
    cr = fCS->size[0] - 1;
    i2 = b_fCS->size[0] * b_fCS->size[1];
    b_fCS->size[0] = fCS_idx_0;
    b_fCS->size[1] = ((i1 - i0) + loop_ub) + 1;
    emxEnsureCapacity((emxArray__common *)b_fCS, i2, (int)sizeof(creal_T));
    br = i1 - i0;
    for (i2 = 0; i2 < br; i2++) {
      for (i3 = 0; i3 < fCS_idx_0; i3++) {
        b_fCS->data[i3 + b_fCS->size[0] * i2] = fCS->data[i3 + fCS->size[0] *
          (i0 + i2)];
      }
    }

    for (i2 = 0; i2 <= loop_ub; i2++) {
      for (i3 = 0; i3 <= cr; i3++) {
        b_fCS->data[i3 + b_fCS->size[0] * ((i2 + i1) - i0)] = fCS->data[i3 +
          fCS->size[0] * i2];
      }
    }

    i0 = d_y->size[0] * d_y->size[1];
    d_y->size[0] = y->size[0];
    d_y->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)d_y, i0, (int)sizeof(creal_T));
    loop_ub = y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      fCS_idx_0 = y->size[0];
      for (i1 = 0; i1 < fCS_idx_0; i1++) {
        d_y->data[i1 + d_y->size[0] * i0].re = y->data[i1 + y->size[0] * i0];
        d_y->data[i1 + d_y->size[0] * i0].im = 0.0;
      }
    }

    loop_ub = b_fCS->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      fCS_idx_0 = d_y->size[1];
      for (i1 = 0; i1 < fCS_idx_0; i1++) {
        re = 0.0;
        im = 0.0;
        cr = b_fCS->size[1];
        for (i2 = 0; i2 < cr; i2++) {
          c = b_fCS->data[i0 + b_fCS->size[0] * i2].re * d_y->data[i2 +
            d_y->size[0] * i1].re - b_fCS->data[i0 + b_fCS->size[0] * i2].im *
            d_y->data[i2 + d_y->size[0] * i1].im;
          fCS_im = b_fCS->data[i0 + b_fCS->size[0] * i2].re * d_y->data[i2 +
            d_y->size[0] * i1].im + b_fCS->data[i0 + b_fCS->size[0] * i2].im *
            d_y->data[i2 + d_y->size[0] * i1].re;
          re += c;
          im += fCS_im;
        }

        c = c_y->data[i0 + c_y->size[0] * i1];
        if (im == 0.0) {
          pMC->data[(i0 + pMC->size[0] * i1) + pMC->size[0] * pMC->size[1] *
            ((int)(1.0F + (float)p) - 1)].re = re / c;
          pMC->data[(i0 + pMC->size[0] * i1) + pMC->size[0] * pMC->size[1] *
            ((int)(1.0F + (float)p) - 1)].im = 0.0;
        } else if (re == 0.0) {
          pMC->data[(i0 + pMC->size[0] * i1) + pMC->size[0] * pMC->size[1] *
            ((int)(1.0F + (float)p) - 1)].re = 0.0;
          pMC->data[(i0 + pMC->size[0] * i1) + pMC->size[0] * pMC->size[1] *
            ((int)(1.0F + (float)p) - 1)].im = im / c;
        } else {
          pMC->data[(i0 + pMC->size[0] * i1) + pMC->size[0] * pMC->size[1] *
            ((int)(1.0F + (float)p) - 1)].re = re / c;
          pMC->data[(i0 + pMC->size[0] * i1) + pMC->size[0] * pMC->size[1] *
            ((int)(1.0F + (float)p) - 1)].im = im / c;
        }
      }
    }

    p++;
  }

  emxFree_creal_T(&d_y);
  emxFree_creal_T(&b_fCS);
  emxFree_real_T(&c_y);
  emxFree_real_T(&b_y);
  emxFree_real_T(&y);
}

/* End of code generation (permMC.c) */
