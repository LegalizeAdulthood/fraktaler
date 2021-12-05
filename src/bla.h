// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "complex.h"
#include "types.h"

template <typename real>
struct bla
{
  complex<real> A, B;
  real r2;
  count_t l;
};

template<typename real>
struct blas
{
  count_t M;
  count_t L;
  struct bla<real> **b;

  blas(const count_t M, const complex<real> *Z, const formulaC *formula, const real h, const real k, const real L, progress_t *progress, bool *running);

  inline ~blas()
  {
    if (b)
    {
      delete[] b[0];
    }
    delete[] b;
  }

  const struct bla<real> *lookup(const count_t m, const real z2) const;
};
