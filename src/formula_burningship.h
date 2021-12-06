// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <string>

#include "bla.h"
#include "complex.h"
#include "dual.h"
#include "floatexp.h"
#include "formula.h"

struct reference_burningship : public reference
{
  mpfr_t Zx2, Zy2, Z2, z_0;
  reference_burningship(const mpfr_t Cx0, const mpfr_t Cy0);
  virtual ~reference_burningship();
  virtual void step();
};

reference_burningship::reference_burningship(const mpfr_t Cx0, const mpfr_t Cy0)
: reference(Cx0, Cy0)
{
  mpfr_prec_t prec = std::max(mpfr_get_prec(Cx0), mpfr_get_prec(Cy0));
  mpfr_init2(Zx2, prec); mpfr_set_ui(Zx2, 0, MPFR_RNDN);
  mpfr_init2(Zy2, prec); mpfr_set_ui(Zy2, 0, MPFR_RNDN);
  mpfr_init2(Z2,  prec); mpfr_set_ui(Z2,  0, MPFR_RNDN);
  mpfr_init2(z_0, prec); mpfr_set_ui(z_0, 0, MPFR_RNDN);
}

reference_burningship::~reference_burningship()
{
  mpfr_clear(Zx2);
  mpfr_clear(Zy2);
  mpfr_clear(Z2);
  mpfr_clear(z_0);
}

void reference_burningship::step()
{
  // z = |x| + i |y|
  mpfr_abs(Zx, Zx, MPFR_RNDN);
  mpfr_abs(Zy, Zy, MPFR_RNDN);
  // z = z^2 + c
  mpfr_add(z_0, Zx, Zy, MPFR_RNDN);
  mpfr_sqr(z_0, z_0, MPFR_RNDN);
  mpfr_sub(Zx, Zx2, Zy2, MPFR_RNDN);
  mpfr_sub(Zy, z_0, Z2, MPFR_RNDN);
  mpfr_add(Zx, Zx, Cx, MPFR_RNDN);
  mpfr_add(Zy, Zy, Cy, MPFR_RNDN);
  mpfr_sqr(Zx2, Zx, MPFR_RNDN);
  mpfr_sqr(Zy2, Zy, MPFR_RNDN);
  mpfr_add(Z2, Zx2, Zy2, MPFR_RNDN);
}

template <typename real>
inline constexpr blaR2<real> bla_burningship(const real &h, const real &k, const real &L, const complex<real> &Z) noexcept
{
  using std::abs;
  using std::min;
  using std::max;
  const int s = 2 * sgn(Z.x) * sgn(Z.y);
  const mat2<real> A(2 * Z.x, -2 * Z.y, s * Z.y, s * Z.x); // FIXME verify
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

template <typename T, typename t>
inline constexpr complex<t> ptb_burningship(const complex<T> &C, const complex<T> &Z, const complex<t> &c, const complex<t> &z) noexcept
{
  (void) C;
  t xn = (2 * Z.x + z.x) * z.x - (2 * Z.y + z.y) * z.y + c.x;
  t yn = 2 * diffabs(Z.x * Z.y, Z.x * z.y + z.x * (Z.y + z.y)) + c.y;
  return complex<t>(xn, yn);
}

struct formulaR2_burningship : public formulaR2
{
  formulaR2_burningship()
  : formulaR2()
  {
  }
  virtual ~formulaR2_burningship()
  {
  }
  virtual std::string name() const noexcept
  {
    return "Burning Ship";
  }
  virtual reference *new_reference(const mpfr_t Cx, const mpfr_t Cy) const
  {
    return new reference_burningship(Cx, Cy);
  }

  virtual blaR2<float> bla1(const float h, const float k, const float L, const complex<float> Z) const noexcept
  {
    return bla_burningship(h, k, L, Z);
  }
  virtual blaR2<double> bla1(const double h, const double k, const double L, const complex<double> Z) const noexcept
  {
    return bla_burningship(h, k, L, Z);
  }
  virtual blaR2<long double> bla1(const long double h, const long double k, const long double L, const complex<long double> Z) const noexcept
  {
    return bla_burningship(h, k, L, Z);
  }
  virtual blaR2<floatexp> bla1(const floatexp h, const floatexp k, const floatexp L, const complex<floatexp> Z) const noexcept
  {
    return bla_burningship(h, k, L, Z);
  }

  virtual complex<float> perturb(const complex<float> &C, const complex<float> &Z, const complex<float> &c, const complex<float> &z) const noexcept
  {
    return ptb_burningship(C, Z, c, z);
  }
  virtual complex<double> perturb(const complex<double> &C, const complex<double> &Z, const complex<double> &c, const complex<double> &z) const noexcept
  {
    return ptb_burningship(C, Z, c, z);
  }
  virtual complex<long double> perturb(const complex<long double> &C, const complex<long double> &Z, const complex<long double> &c, const complex<long double> &z) const noexcept
  {
    return ptb_burningship(C, Z, c, z);
  }
  virtual complex<floatexp> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const complex<floatexp> &c, const complex<floatexp> &z) const noexcept
  {
    return ptb_burningship(C, Z, c, z);
  }

  virtual complex<dual<2, float>> perturb(const complex<float> &C, const complex<float> &Z, const complex<dual<2, float>> &c, const complex<dual<2, float>> &z) const noexcept
  {
    return ptb_burningship(C, Z, c, z);
  }
  virtual complex<dual<2, double>> perturb(const complex<double> &C, const complex<double> &Z, const complex<dual<2, double>> &c, const complex<dual<2, double>> &z) const noexcept
  {
    return ptb_burningship(C, Z, c, z);
  }
  virtual complex<dual<2, long double>> perturb(const complex<long double> &C, const complex<long double> &Z, const complex<dual<2, long double>> &c, const complex<dual<2, long double>> &z) const noexcept
  {
    return ptb_burningship(C, Z, c, z);
  }
  virtual complex<dual<2, floatexp>> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const complex<dual<2, floatexp>> &c, const complex<dual<2, floatexp>> &z) const noexcept
  {
    return ptb_burningship(C, Z, c, z);
  }
};
