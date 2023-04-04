// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
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

template <typename real>
blaR2<real> merge(const blaR2<real> &y, const blaR2<real> &x, const real &c)
{
  using std::min;
  using std::max;
  using std::sqrt, ::sqrt;
  const count_t l = x.l + y.l;
  const mat2<real> A = y.A * x.A;
  const mat2<real> B = y.A * x.B + y.B;
  const real xA = sup(x.A);
  const real xB = sup(x.B);
  const real r = min(sqrt(x.r2), max(real(0), (sqrt(y.r2) - xB * c) / xA));
  const real r2 = r * r;
  blaR2<real> b = { A, B, r2, l };
  return b;
}

template<typename real>
struct blasR2
{
  count_t M;
  count_t L;
  std::vector<std::vector<blaR2<real>>> b;

  blasR2(const std::vector<complex<real>> &Z, const std::vector<std::vector<opcode>> &opss, const std::vector<int> &degrees, const count_t phase, const real h, const real k, const real stepcount, int skip_levels, volatile progress_t *progress, volatile bool *running);

  const struct blaR2<real> *lookup(const count_t m, const real z2) const noexcept;
};
