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

struct reference_mandelbrot : public reference
{
  mpfr_t Zx2, Zy2, Z2, z_0;
  reference_mandelbrot(const mpfr_t Cx0, const mpfr_t Cy0);
  virtual ~reference_mandelbrot();
  virtual void step();
};

reference_mandelbrot::reference_mandelbrot(const mpfr_t Cx0, const mpfr_t Cy0)
: reference(Cx0, Cy0)
{
  mpfr_prec_t prec = std::max(mpfr_get_prec(Cx0), mpfr_get_prec(Cy0));
  mpfr_init2(Zx2, prec); mpfr_set_ui(Zx2, 0, MPFR_RNDN);
  mpfr_init2(Zy2, prec); mpfr_set_ui(Zy2, 0, MPFR_RNDN);
  mpfr_init2(Z2,  prec); mpfr_set_ui(Z2,  0, MPFR_RNDN);
  mpfr_init2(z_0, prec); mpfr_set_ui(z_0, 0, MPFR_RNDN);
}

reference_mandelbrot::~reference_mandelbrot()
{
  mpfr_clear(Zx2);
  mpfr_clear(Zy2);
  mpfr_clear(Z2);
  mpfr_clear(z_0);
}

void reference_mandelbrot::step()
{
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

template <typename real> bla<real> bla_mandelbrot(const real &h, const real &k, const real &L, const complex<real> &Z)
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
  bla<real> b = { A, B, r2, l };
  return b;
}

template <typename T, typename t>
t ptb_mandelbrot(const T &C, const T &Z, const t &c, const t &z)
{
  return (2 * Z + z) * z + c;
}

struct formulaC_mandelbrot : public formulaC
{
  formulaC_mandelbrot()
  : formulaC()
  {
  }
  virtual ~formulaC_mandelbrot()
  {
  }
  virtual std::string name() const
  {
    return "Mandelbrot";
  }
  virtual reference *new_reference(const mpfr_t Cx, const mpfr_t Cy) const
  {
    return new reference_mandelbrot(Cx, Cy);
  }

  virtual bla<float> bla1(const float h, const float k, const float L, const complex<float> Z) const
  {
    return bla_mandelbrot(h, k, L, Z);
  }
  virtual bla<double> bla1(const double h, const double k, const double L, const complex<double> Z) const
  {
    return bla_mandelbrot(h, k, L, Z);
  }
  virtual bla<long double> bla1(const long double h, const long double k, const long double L, const complex<long double> Z) const
  {
    return bla_mandelbrot(h, k, L, Z);
  }
  virtual bla<floatexp> bla1(const floatexp h, const floatexp k, const floatexp L, const complex<floatexp> Z) const
  {
    return bla_mandelbrot(h, k, L, Z);
  }

  virtual complex<float> perturb(const complex<float> &C, const complex<float> &Z, const complex<float> &c, const complex<float> &z) const
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual complex<double> perturb(const complex<double> &C, const complex<double> &Z, const complex<double> &c, const complex<double> &z) const
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual complex<long double> perturb(const complex<long double> &C, const complex<long double> &Z, const complex<long double> &c, const complex<long double> &z) const
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual complex<floatexp> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const complex<floatexp> &c, const complex<floatexp> &z) const
  {
    return ptb_mandelbrot(C, Z, c, z);
  }

  virtual dual<1, complex<float>> perturb(const complex<float> &C, const complex<float> &Z, const dual<1, complex<float>> &c, const dual<1, complex<float>> &z) const
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual dual<1, complex<double>> perturb(const complex<double> &C, const complex<double> &Z, const dual<1, complex<double>> &c, const dual<1, complex<double>> &z) const
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual dual<1, complex<long double>> perturb(const complex<long double> &C, const complex<long double> &Z, const dual<1, complex<long double>> &c, const dual<1, complex<long double>> &z) const
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual dual<1, complex<floatexp>> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const dual<1, complex<floatexp>> &c, const dual<1, complex<floatexp>> &z) const
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
};
