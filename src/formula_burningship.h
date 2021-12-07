// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "bla.h"
#include "complex.h"
#include "floatexp.h"

const char burningship_name[] = "Burning Ship";

template <typename T>
complex<T> burningship_plain(const complex<T> &C, const complex<T> &Z)
{
  return sqr(complex<T>(abs(Z.x), abs(Z.y))) + C;
}

template <typename T, typename t>
complex<t> burningship_perturb(const complex<T> &C, const complex<T> &Z, const complex<t> &c, const complex<t> &z) noexcept
{
  (void) C;
  t x = (2 * Z.x + z.x) * z.x - (2 * Z.y + z.y) * z.y + c.x;
  t y = 2 * diffabs(Z.x * Z.y, Z.x * z.y + z.x * (Z.y + z.y)) + c.y;
  return complex<t>(x, y);
}

template <typename real>
blaR2<real> burningship_bla(const real &h, const real &k, const real &L, const complex<real> &Z) noexcept
{
  using std::abs;
  using std::min;
  using std::max;
  const int s = 2 * sgn(Z.x) * sgn(Z.y);
  const mat2<real> A(2 * Z.x, -2 * Z.y, s * Z.y, s * Z.x);
  const mat2<real> B(1);
  const real mZ = min(abs(Z.x), abs(Z.y)) / 2; // FIXME arbitrary factor
  const real mA = abs(A);
  const real mB = abs(B);
  const real r = max(real(0), (mZ - mB * h * k) / (mA + 1)) / L;
  const real r2 = r * r;
  const count_t l = 1;
  blaR2<real> b = { A, B, r2, l };
  return b;
}
