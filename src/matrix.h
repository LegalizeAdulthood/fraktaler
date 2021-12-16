// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <algorithm>

#include <glm/glm.hpp>

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
  inline constexpr mat2(const glm::mat3 &m) // FIXME check transpose?
  {
    x[0][0] = m[0][0];
    x[0][1] = m[1][0];
    x[1][0] = m[0][1];
    x[1][1] = m[1][1];
  }

  inline constexpr mat2 &operator+=(const mat2 &b)
  {
    x[0][0] += b.x[0][0];
    x[0][1] += b.x[0][1];
    x[1][0] += b.x[1][0];
    x[1][1] += b.x[1][1];
    return *this;
  }
  inline constexpr mat2 &operator/=(const real &b)
  {
    x[0][0] /= b;
    x[0][1] /= b;
    x[1][0] /= b;
    x[1][1] /= b;
    return *this;
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

template <typename S, typename T>
inline constexpr complex<T> operator*(const complex<T> &z, const mat2<S> &m) noexcept
{
  return complex<T>(m.x[0][0] * z.x + m.x[1][0] * z.y, m.x[0][1] * z.x + m.x[1][1] * z.y);
}

template <typename S, typename T>
inline constexpr mat2<T> operator/(const mat2<T> &m, const S &b) noexcept
{
  return mat2<T>(m.x[0][0] / b, m.x[0][1] / b, m.x[1][0] / b, m.x[1][1] / b);
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

template <typename real>
inline constexpr mat2<real> inverse(const mat2<real> &a)
{
  return mat2<real>(a.x[1][1], -a.x[0][1], -a.x[1][0], a.x[0][0]) / determinant(a);
}

template <typename real>
inline constexpr mat2<real> rotation(const real &radians)
{
  using std::cos;
  using std::sin;
  return mat2<real>(complex<real>(cos(radians), sin(radians)));
}
