// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <iostream>

#include "complex.h"
#include "floatexp.h"
#include "formula.h"
#include "reference.h"

template <typename real>
count_t run_reference(complex<real> *Z, const count_t MaximumReferenceIterations, reference *ref, progress_t *progress, bool *running)
{
  count_t M = MaximumReferenceIterations - 1;
  // calculate reference in high precision
  for (count_t i = 0; i < MaximumReferenceIterations; ++i)
  {
    // store low precision orbit
    complex<floatexp> z = ref->get();
    Z[i] = complex<real>(real(z.x), real(z.y));
    // step
    ref->step();
    // escape check
    if (ref->escaped() || ! *running)
    {
      M = i;
      break;
    }
    *progress = (i + 1) / progress_t(MaximumReferenceIterations);
  }
  return M;
}

template count_t run_reference(complex<float> *Z, const count_t MaximumReferenceIterations, reference *ref, progress_t *progress, bool *running);
template count_t run_reference(complex<double> *Z, const count_t MaximumReferenceIterations, reference *ref, progress_t *progress, bool *running);
template count_t run_reference(complex<long double> *Z, const count_t MaximumReferenceIterations, reference *ref, progress_t *progress, bool *running);
template count_t run_reference(complex<floatexp> *Z, const count_t MaximumReferenceIterations, reference *ref, progress_t *progress, bool *running);
