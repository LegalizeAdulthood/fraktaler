// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

template <typename real>
struct mat2
{
  real x[2][2];

  // POD
  inline ~mat2() = default;
  inline mat2() = default;
  inline constexpr mat2(const mat2 &m) = default;
  inline constexpr mat2(mat2 &&m) = default;
  inline constexpr mat2 &operator=(const mat2 &m) = default;
  
  inline constexpr mat2(const real &s) noexcept
  {
    x[0][0] = s;
    x[0][1] = 0;
    x[1][0] = 0;
    x[1][1] = s;
  }
  inline constexpr mat2(const complex<real> &z) noexcept
  {
    x[0][0] = z.x;
    x[0][1] = -z.y;
    x[1][0] = z.y;
    x[1][1] = z.x;
  }
  inline constexpr mat2(const real &a, const real &b, const real &c, const real &d) noexcept
  {
    x[0][0] = a;
    x[0][1] = b;
    x[1][0] = c;
    x[1][1] = d;
  }
};

template <typename real>
inline constexpr mat2<real> operator+(const mat2<real> &a, const mat2<real> &b)
{
  return mat2<real>
    ( a.x[0][0] + b.x[0][0]
    , a.x[0][1] + b.x[0][1]
    , a.x[1][0] + b.x[1][0]
    , a.x[1][1] + b.x[1][1]
    );
}

template <typename real>
inline constexpr mat2<real> operator*(const mat2<real> &a, const mat2<real> &b)
{
  return mat2<real>
    ( a.x[0][0] * b.x[0][0] + a.x[0][1] * b.x[1][0]
    , a.x[0][0] * b.x[0][1] + a.x[0][1] * b.x[1][1]
    , a.x[1][0] * b.x[0][0] + a.x[1][1] * b.x[1][0]
    , a.x[1][0] * b.x[0][1] + a.x[1][1] * b.x[1][1]
    );
}

template <typename S, typename T>
inline constexpr complex<T> operator*(const mat2<S> &m, const complex<T> &z) noexcept
{
  return complex<T>(m.x[0][0] * z.x + m.x[0][1] * z.y, m.x[1][0] * z.x + m.x[1][1] * z.y);
}

template <typename real>
inline constexpr mat2<real> transpose(const mat2<real> &a)
{
  return mat2<real>
    ( a.x[0][0]
    , a.x[1][0]
    , a.x[0][1]
    , a.x[1][1]
    );
}

template <typename real>
inline constexpr real trace(const mat2<real> &a)
{
  return a.x[0][0] + a.x[1][1];
}

template <typename real>
inline constexpr real determinant(const mat2<real> &a)
{
  return a.x[0][0] * a.x[1][1] - a.x[0][1] * a.x[1][0];
}

template <typename real>
inline constexpr real norm(const mat2<real> &a)
{
  using std::max;
  using std::sqrt;
  const mat2<real> aTa = transpose(a) * a;
  const real T = trace(aTa);
  const real D = determinant(aTa);
  return (T + sqrt(max(real(0), sqr(T) - 4 * D))) / 2;
}

template <typename real>
inline constexpr real abs(const mat2<real> &a)
{
  using std::sqrt;
  return sqrt(norm(a));
}
