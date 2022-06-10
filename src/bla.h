// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <vector>

#include "complex.h"
#include "matrix.h"
#include "types.h"

struct phybrid;

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
  std::vector<std::vector<blaR2<real>>> b;

  blasR2(const std::vector<complex<real>> &Z, const phybrid &H, const count_t phase, const real h, const real k, const real stepcount, volatile progress_t *progress, volatile bool *running);

  const struct blaR2<real> *lookup(const count_t m, const real z2) const noexcept;
};
