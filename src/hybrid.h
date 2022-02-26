// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <vector>

#include "dual.h"
#include "param.h"

// http://www.burtleburtle.net/bob/hash/integer.html
inline CONSTEXPR uint32_t burtle_hash(uint32_t a) noexcept
{
  a = (a+0x7ed55d16) + (a<<12);
  a = (a^0xc761c23c) ^ (a>>19);
  a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);
  a = (a+0xfd7046c5) + (a<<3);
  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}

inline double radical_inverse(coord_t a, const coord_t base) noexcept
{
  constexpr double one_minus_epsilon = 0.99999999999999989;
  const double base1 = 1.0 / base;
  coord_t reversed = 0;
  double base1n = 1;
  while (a)
  {
    const coord_t next  = a / base;
    const coord_t digit = a - base * next;
    reversed = reversed * base + digit;
    base1n *= base1;
    a = next;
  }
  return std::min(reversed * base1n, one_minus_epsilon);
}

inline double wrap(const double v) noexcept
{
  return v - std::floor(v);
}

inline double triangle(const double a) noexcept
{
  const double b = a * 2 - 1;
  const double c = std::sqrt(std::abs(b));
  const double e = b > 0 ? c - 1 : 1 - c;
  return e;
}

inline void jitter(const coord_t width, const coord_t height, const coord_t i, const coord_t j, const coord_t k, double &x, double &y) noexcept
{
  coord_t ix = (k * height + j) * width + i;
  double h = burtle_hash(ix) / double(0x100000000LL);
  x = triangle(wrap(radical_inverse(k, 2) + h));
  y = triangle(wrap(radical_inverse(k, 3) + h));
}

template <typename T>
inline complex<T> hybrid_plain(const struct phybrid1 &H, const complex<T> &C, const complex<T> &Z)
{
  using std::abs;
  complex<T> W = Z;
  if (H.abs_x) W.x = abs(W.x);
  if (H.abs_y) W.y = abs(W.y);
  if (H.neg_x) W.x = -W.x;
  if (H.neg_y) W.y = -W.y;
  return pow(W, H.power) + C;
}

template <typename T, typename t>
inline constexpr complex<t> hybrid_perturb(const struct phybrid1 &H, const complex<T> &C, const complex<T> &Z, const complex<t> &c, const complex<t> &z) noexcept
{
  (void) C;
  using std::abs;
  T X = Z.x;
  T Y = Z.y;
  t x = z.x;
  t y = z.y;
  complex<t> W = Z + z;
  complex<T> B = Z;
  if (H.abs_x)
  {
    x = diffabs(X, x);
    W.x = abs(W.x);
    B.x = abs(B.x);
  }
  if (H.abs_y)
  {
    y = diffabs(Y, y);
    W.y = abs(W.y);
    B.y = abs(B.y);
  }
  if (H.neg_x)
  {
    x = -x;
    W.x = -W.x;
    B.x = -B.x;
  }
  if (H.neg_y)
  {
    y = -y;
    W.y = -W.y;
    B.y = -B.y;
  }
  complex<t> P(x, y);
  complex<t> S(0);
  for (int i = 0; i <= H.power - 1; ++i)
  {
    int j = H.power - 1 - i;
    S += pow(W, i) * pow(B, j);
  }
  return P * S + c;
}

template <typename real>
inline constexpr blaR2<real> hybrid_bla(const struct phybrid1 &H, const real &h, const real &k, const real &L, const complex<real> &Z) noexcept
{
  using std::abs, ::abs;
  using std::min;
  using std::max;
  dual<2, real> x(Z.x); x.dx[0] = 1;
  dual<2, real> y(Z.y); y.dx[1] = 1;
  complex<dual<2, real>> W(x, y);
  complex<dual<2, real>> C(0, 0);
  W = hybrid_plain(H, C, W);
  const mat2<real> A(W.x.dx[0], W.x.dx[1], W.y.dx[0], W.y.dx[1]);
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

template <typename t> void hybrid_blas(std::vector<blasR2<t>> &B, const std::vector<std::vector<complex<t>>> &Z, const phybrid &H, t h, t k, t L, volatile progress_t *progress, volatile bool *running);
template <typename t> count_t hybrid_reference(complex<t> *Zp, const struct phybrid &H, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template <typename t> void hybrid_references(std::vector<std::vector<complex<t>>> &Zp, const struct phybrid &H, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template <typename t> void hybrid_render(map &out, stats &sta, const phybrid &H, const std::vector<blasR2<t>> &bla, const count_t subframe, const param &par, const t Zoom, const complex<t> offset, const std::vector<std::vector<complex<t>>> &Zp, volatile progress_t *progress, volatile bool *running);
template <typename t> count_t hybrid_reference(complex<t> *Zp, const struct phybrid &H, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template <typename t> count_t hybrid_period(const phybrid &H, const std::vector<std::vector<complex<t>>> &Zp, const complex<floatexp> &c0, const count_t &N, const count_t &ReferencePeriod, const floatexp &s, const mat2<double> &K, volatile progress_t *progress, volatile bool *running);
bool hybrid_center(const phybrid &h, complex<mpreal> &C0, const count_t period, volatile progress_t *progress, volatile bool *running);
bool hybrid_size(floatexp &s, mat2<double> &K, const phybrid &h, const complex<mpreal> &C, count_t period, volatile progress_t *progress, volatile bool *running);
