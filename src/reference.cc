// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <iostream>

#include "reference.h"
#include "floatexp.h"
#include "complex.h"

template <typename real>
count_t reference(complex<real> *Z, const count_t MaximumReferenceIterations, const mpfr_t Cx, const mpfr_t Cy, progress_t *progress, bool *running)
{
  count_t M = MaximumReferenceIterations - 1;
  mpfr_prec_t prec = std::max(mpfr_get_prec(Cx), mpfr_get_prec(Cy));
  mpfr_t Zx, Zy, Zx2, Zy2, z_0, Z2, ER2;
  // allocate variables
  mpfr_init2(Zx, prec);
  mpfr_init2(Zy, prec);
  mpfr_init2(Zx2, prec);
  mpfr_init2(Zy2, prec);
  mpfr_init2(z_0, prec);
  mpfr_init2(Z2, prec);
  mpfr_init2(ER2, 53);
  // initialize
  mpfr_set_d(ER2, 4, MPFR_RNDN);
  mpfr_set_d(Zx, 0, MPFR_RNDN);
  mpfr_set_d(Zy, 0, MPFR_RNDN);
  mpfr_sqr(Zx2, Zx, MPFR_RNDN);
  mpfr_sqr(Zy2, Zy, MPFR_RNDN);
  mpfr_add(Z2, Zx2, Zy2, MPFR_RNDN);
  // calculate reference in high precision
  for (count_t i = 0; i < MaximumReferenceIterations; ++i)
  {
    // store low precision orbit
    long ex = 0, ey = 0;
    double zx = mpfr_get_d_2exp(&ex, Zx, MPFR_RNDN);
    double zy = mpfr_get_d_2exp(&ey, Zy, MPFR_RNDN);
    Z[i] = complex<real>(real(floatexp(zx, ex)), real(floatexp(zy, ey)));
    // z = z^2 + c
    mpfr_add(z_0, Zx, Zy, MPFR_RNDN);
    mpfr_sqr(z_0, z_0, MPFR_RNDN);
    mpfr_sub(Zx, Zx2, Zy2, MPFR_RNDN);
    mpfr_sub(Zy, z_0, Z2, MPFR_RNDN);
    mpfr_add(Zx, Zx, Cx, MPFR_RNDN);
    mpfr_add(Zy, Zy, Cy, MPFR_RNDN);
    mpfr_sqr(Zx2, Zx, MPFR_RNDN);
    mpfr_sqr(Zy2, Zy, MPFR_RNDN);
    mpfr_add(Z2, Zx2, Zy2, MPFR_RNDN);
    // if |z|^2 > er2
    if (mpfr_greater_p(Z2, ER2) || ! *running)
    {
      M = i;
      break;
    }
    *progress = (i + 1) / progress_t(MaximumReferenceIterations);
  }
  // deallocate variables
  mpfr_clear(Zx);
  mpfr_clear(Zy);
  mpfr_clear(Zx2);
  mpfr_clear(Zy2);
  mpfr_clear(z_0);
  mpfr_clear(Z2);
  mpfr_clear(ER2);
  return M;
}

template count_t reference(complex<float> *Z, const count_t MaximumReferenceIterations, const mpfr_t Cx, const mpfr_t Cy, progress_t *progress, bool *running);
template count_t reference(complex<double> *Z, const count_t MaximumReferenceIterations, const mpfr_t Cx, const mpfr_t Cy, progress_t *progress, bool *running);
template count_t reference(complex<long double> *Z, const count_t MaximumReferenceIterations, const mpfr_t Cx, const mpfr_t Cy, progress_t *progress, bool *running);
template count_t reference(complex<floatexp> *Z, const count_t MaximumReferenceIterations, const mpfr_t Cx, const mpfr_t Cy, progress_t *progress, bool *running);
