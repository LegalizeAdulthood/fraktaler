// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "complex.h"
#include "matrix.h"
#include "types.h"
#include "float128.h"

template <typename real>
struct blaC
{
  complex<real> A, B;
  real r2;
  count_t l;
};

template<typename real>
struct blasC
{
  count_t M;
  count_t L;
  struct blaC<real> **b;

  blasC(const count_t M, const complex<real> *Z, const formulaCbase *formula, const real h, const real k, const real L, progress_t *progress, bool *running);

  inline ~blasC()
  {
    if (b)
    {
      delete[] b[0];
    }
    delete[] b;
  }

  const struct blaC<real> *lookup(const count_t m, const real z2) const;
};

template <typename real>
struct blaR2
{
  mat2<real> A, B;
  real r2;
  count_t l;
};

template<typename real>
struct blasR2
{
  count_t M;
  count_t L;
  struct blaR2<real> **b;

  blasR2(const count_t M, const complex<real> *Z, const formulaR2base *formula, const real h, const real k, const real L, progress_t *progress, bool *running);

  inline ~blasR2()
  {
    if (b)
    {
      delete[] b[0];
    }
    delete[] b;
  }

  const struct blaR2<real> *lookup(const count_t m, const real z2) const;
};
