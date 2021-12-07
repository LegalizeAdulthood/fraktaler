// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "bla.h"
#include "complex.h"

const char mandelbrot_name[] = "Mandelbrot";

template <typename T>
T mandelbrot_plain(const T &C, const T &Z)
{
  return sqr(Z) + C;
}

template <typename T, typename t>
t mandelbrot_perturb(const T &C, const T &Z, const t &c, const t &z) noexcept
{
  (void) C;
  return (2 * Z + z) * z + c;
}

template <typename real>
blaC<real> mandelbrot_bla(const real &h, const real &k, const real &L, const complex<real> &Z) noexcept
{
  using std::max;
  const complex<real> A(2 * Z);
  const complex<real> B(1);
  const real mZ = abs(Z);
  const real mA = abs(A);
  const real mB = abs(B);
  const real r = max(real(0), (mZ - mB * h * k) / (mA + 1)) / L;
  const real r2 = r * r;
  const count_t l = 1;
  blaC<real> b = { A, B, r2, l };
  return b;
}
