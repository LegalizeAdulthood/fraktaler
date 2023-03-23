// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <vector>

#include "dual.h"
#include "param.h"

struct tile;

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

inline void jitter(const coord_t width, const coord_t height, const coord_t frame, const coord_t i, const coord_t j, const coord_t k, double &x, double &y) noexcept
{
  coord_t ix = (frame * height + j) * width + i;
  double h = burtle_hash(ix) / double(0x100000000LL);
  x = triangle(wrap(radical_inverse(k, 2) + h));
  y = triangle(wrap(radical_inverse(k, 3) + h));
}

template <typename T>
inline void hybrid_plain(const opcode &op, const complex<T> &C, complex<T> &Z, complex<T> &Z_stored)
{
  using std::abs;
  switch (op)
  {
    case op_add: Z += C; break;
    case op_store: Z_stored = Z; break;
    case op_sqr: Z = sqr(Z); break;
    case op_mul: Z *= Z_stored; break;
    case op_absx: Z.x = abs(Z.x); break;
    case op_absy: Z.y = abs(Z.y); break;
    case op_negx: Z.x = -(Z.x); break;
    case op_negy: Z.y = -(Z.y); break;
  }
}

template <typename T>
inline complex<T> hybrid_plain(const std::vector<opcode> &ops, const complex<T> &C, const complex<T> &Z0)
{
  complex<T> Z = Z0;
  complex<T> Z_stored = Z;
  for (const auto & op : ops)
  {
    hybrid_plain(op, C, Z, Z_stored);
  }
  return Z;
}

template <typename T, typename t>
inline constexpr void hybrid_perturb(const opcode &op, const complex<T> &C, complex<T> &Z, complex<T> &Z_stored, const complex<t> &c, complex<t> &z, complex<t> &z_stored) noexcept
{
  using std::abs;
  switch (op)
  {
    case op_add: z += c; Z += C; break;
    case op_store: z_stored = z; Z_stored = Z; break;
    case op_sqr: z = (2 * Z + z) * z; Z = sqr(Z); break;
    case op_mul: z = Z * z_stored + Z_stored * z + z_stored * z; Z *= Z_stored; break;
    case op_absx: z.x = diffabs(Z.x, z.x); Z.x = abs(Z.x); break;
    case op_absy: z.y = diffabs(Z.y, z.y); Z.y = abs(Z.y); break;
    case op_negx: z.x = -(z.x); Z.x = -(Z.x); break;
    case op_negy: z.y = -(z.y); Z.y = -(Z.y); break;
  }
}

template <typename T, typename t>
inline constexpr complex<t> hybrid_perturb(const std::vector<opcode> &ops, const complex<T> &C, const complex<T> &Z0, const complex<t> &c, const complex<t> &z0, bool &rebased) noexcept
{
  complex<T> Z = Z0;
  complex<t> z = z0;
  complex<T> Z_stored = Z;
  complex<t> z_stored = z;
  for (const auto & op : ops)
  {
    complex<t> Zz = Z + z;
    if (norm(Zz) < norm(z))
    {
      z = Zz;
      Z = 0;
      rebased = true;
    }
    hybrid_perturb(op, C, Z, Z_stored, c, z, z_stored);
  }
  return z;
}

template <typename real>
inline constexpr blaR2<real> hybrid_bla(const opcode &op, const real &e, complex<real> &Z, complex<real> &Z_stored) noexcept
{
  using std::abs;
  dual<2, real> x(Z.x); x.dx[0] = 1;
  dual<2, real> y(Z.y); y.dx[1] = 1;
  complex<dual<2, real>> W(x, y);
  dual<2, real> x_stored(Z_stored.x); x_stored.dx[0] = 1;
  dual<2, real> y_stored(Z_stored.y); y_stored.dx[1] = 1;
  complex<dual<2, real>> W_stored(x_stored, y_stored);
  complex<dual<2, real>> C(0, 0); // FIXME
  hybrid_plain(op, C, W, W_stored);
  Z.x = W.x.x;
  Z.y = W.y.x;
  Z_stored.x = W_stored.x.x;
  Z_stored.y = W_stored.y.x;
  const mat2<real> A(W.x.dx[0], W.x.dx[1], W.y.dx[0], W.y.dx[1]);
  const mat2<real> O(0);
  const mat2<real> I(1);
  real infinity = real(1) / real(0);
  switch (op)
  {
    case op_add: return blaR2<real>{ A, I, infinity, 1 };
    case op_store:return blaR2<real>{ A, O, infinity, 1 };
    case op_sqr: return blaR2<real>{ A, O, sqr(e * inf(A)), 1 };
    case op_mul: return blaR2<real>{ A, O, sqr(e * inf(A)), 1 }; // FIXME verify
    case op_absx: return blaR2<real>{ A, O, sqr(abs(Z.x) / 2), 1 }; // FIXME arbitrary factor
    case op_absy: return blaR2<real>{ A, O, sqr(abs(Z.y) / 2), 1 }; // FIXME arbitrary factor
    case op_negx: return blaR2<real>{ A, O, infinity, 1 };
    case op_negy: return blaR2<real>{ A, O, infinity, 1 };
  }
  assert(! "reachable");
}

template <typename real>
inline constexpr blaR2<real> hybrid_bla(const std::vector<opcode> &ops, const real &h, const real &k, const real &L, const complex<real> &Z) noexcept
{
  using std::abs;
  complex<real> W (Z);
  complex<real> W_stored (Z);
  const mat2<real> O(0);
  const mat2<real> I(1);
  real infinity = real(1) / real(0);  
  blaR2<real> b = { I, O, infinity, 1 };
  real hk = h * k;
  real e = 1 / L;
  for (const auto & op : ops)
  {
    b = merge(hybrid_bla(op, e, W, W_stored), b, hk);
  }
  b.l = 1;
  return b;
}

template <typename t> bool hybrid_blas(std::vector<blasR2<t>> &B, const std::vector<std::vector<complex<t>>> &Z, const std::vector<std::vector<opcode>> &opss, t h, t k, t L, volatile progress_t *progress, volatile bool *running);
template <typename t> count_t hybrid_reference(complex<t> *Zp, const std::vector<std::vector<opcode>> &opss, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template <typename t> void hybrid_references(std::vector<std::vector<complex<t>>> &Zp, const std::vector<std::vector<opcode>> &opss, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template <typename t> count_t hybrid_reference(complex<t> *Zp, const std::vector<std::vector<opcode>> &opss, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template <typename t> count_t hybrid_period(const std::vector<std::vector<opcode>> &opss, const std::vector<std::vector<complex<t>>> &Zp, const complex<floatexp> &c0, const count_t &N, const floatexp &s, const mat2<double> &K, volatile progress_t *progress, volatile bool *running);
bool hybrid_center(const std::vector<std::vector<opcode>> &opss, complex<mpreal> &C0, const count_t period, volatile progress_t *progress, volatile bool *running);
bool hybrid_size(floatexp &s, mat2<double> &K, const std::vector<std::vector<opcode>> &opss, const std::vector<int> &degrees, const complex<mpreal> &C, count_t period, volatile progress_t *progress, volatile bool *running);

std::string hybrid_perturb_opencl(const std::vector<std::vector<opcode>> &opss, const std::vector<int> &degrees);

template <typename T>
bool hybrid_render(coord_t frame, coord_t x0, coord_t y0, coord_t x1, coord_t y1, coord_t subframe, tile *data, const param &par, const std::vector<std::vector<complex<T>>> &ref, const std::vector<blasR2<T>> &bla, volatile bool *running);
