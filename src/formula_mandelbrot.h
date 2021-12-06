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

template <typename real>
count_t period_mandelbrot(const complex<real> *Zp, const count_t M, const complex<real> c, const count_t N, const real &s, const mat2<real> &K, progress_t *progress, bool *running)
{
  // perturbed version of knighty's Taylor ball with Zhuoran's rebasing
  complex<real> z(0), dz(0), Z(0), Zz(0);
  const real r(s);
  real rdz = 0, rz = 0, rzd = 0, rr = 0, rrdz = 0, Ei = 0, ER = 65536;
  count_t n = 0, m = 0;
  while (n < N && rz < ER)
  {
    progress[0] = n / progress_t(N);
    if (! *running)
    {
      break;
    }
    // step
    Ei = rdz * rdz + (2 * rz + r * (2 * rdz + r * Ei)) * Ei;
    dz = 2 * Zz * dz + 1;
    z = (2 * Z + z) * z + c;
    n++;
    m++;
    Z = Zp[m];
    Zz = Z + z;
    // calculate ball
    rdz = abs(dz);
    rz = abs(Zz);
    rzd = abs(z);
    rr = r * (rdz + r * Ei);
    rrdz = r * (rdz - r  * Ei);
    // rebase
    if (rz < rzd || m == M - 1)
    {
      z = Zz;
      m = 0;
      Z = Zp[m];
    }
    if (rz - rr > 2)
    {
      // escaped
      break;
    }
    if (rz <= rr)
    {
      count_t period = n;
      if (! (rz <= rrdz))
      {
        period = -period;
      }
      return period;
    }
  }
  return 0;
}

template <typename real>
blaC<real> bla_mandelbrot(const real &h, const real &k, const real &L, const complex<real> &Z) noexcept
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

template <typename T, typename t>
t ptb_mandelbrot(const T &C, const T &Z, const t &c, const t &z) noexcept
{
  (void) C;
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
  virtual std::string name() const noexcept
  {
    return "Mandelbrot";
  }
  virtual reference *new_reference(const mpfr_t Cx, const mpfr_t Cy) const
  {
    return new reference_mandelbrot(Cx, Cy);
  }

  virtual count_t period(const complex<float> *Zp, const count_t M, const complex<float> c, const count_t N, const float &s, const mat2<float> &K, progress_t *progress, bool *running) const noexcept
  {
    return period_mandelbrot(Zp, M, c, N, s, K, progress, running);
  }
  virtual count_t period(const complex<double> *Zp, const count_t M, const complex<double> c, const count_t N, const double &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept
  {
    return period_mandelbrot(Zp, M, c, N, s, K, progress, running);
  }
  virtual count_t period(const complex<long double> *Zp, const count_t M, const complex<long double> c, const count_t N, const long double &s, const mat2<long double> &K, progress_t *progress, bool *running) const noexcept
  {
    return period_mandelbrot(Zp, M, c, N, s, K, progress, running);
  }
  virtual count_t period(const complex<floatexp> *Zp, const count_t M, const complex<floatexp> c, const count_t N, const floatexp &s, const mat2<floatexp> &K, progress_t *progress, bool *running) const noexcept
  {
    return period_mandelbrot(Zp, M, c, N, s, K, progress, running);
  }

  virtual blaC<float> bla1(const float h, const float k, const float L, const complex<float> Z) const noexcept
  {
    return bla_mandelbrot(h, k, L, Z);
  }
  virtual blaC<double> bla1(const double h, const double k, const double L, const complex<double> Z) const noexcept
  {
    return bla_mandelbrot(h, k, L, Z);
  }
  virtual blaC<long double> bla1(const long double h, const long double k, const long double L, const complex<long double> Z) const noexcept
  {
    return bla_mandelbrot(h, k, L, Z);
  }
  virtual blaC<floatexp> bla1(const floatexp h, const floatexp k, const floatexp L, const complex<floatexp> Z) const noexcept
  {
    return bla_mandelbrot(h, k, L, Z);
  }

  virtual complex<float> perturb(const complex<float> &C, const complex<float> &Z, const complex<float> &c, const complex<float> &z) const noexcept
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual complex<double> perturb(const complex<double> &C, const complex<double> &Z, const complex<double> &c, const complex<double> &z) const noexcept
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual complex<long double> perturb(const complex<long double> &C, const complex<long double> &Z, const complex<long double> &c, const complex<long double> &z) const noexcept
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual complex<floatexp> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const complex<floatexp> &c, const complex<floatexp> &z) const noexcept
  {
    return ptb_mandelbrot(C, Z, c, z);
  }

  virtual dual<1, complex<float>> perturb(const complex<float> &C, const complex<float> &Z, const dual<1, complex<float>> &c, const dual<1, complex<float>> &z) const noexcept
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual dual<1, complex<double>> perturb(const complex<double> &C, const complex<double> &Z, const dual<1, complex<double>> &c, const dual<1, complex<double>> &z) const noexcept
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual dual<1, complex<long double>> perturb(const complex<long double> &C, const complex<long double> &Z, const dual<1, complex<long double>> &c, const dual<1, complex<long double>> &z) const noexcept
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
  virtual dual<1, complex<floatexp>> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const dual<1, complex<floatexp>> &c, const dual<1, complex<floatexp>> &z) const noexcept
  {
    return ptb_mandelbrot(C, Z, c, z);
  }
};
