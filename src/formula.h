// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <vector>

#include <mpfr.h>

#include "types.h"

struct reference
{
  mpfr_t Cx, Cy, Zx, Zy;
  reference(const mpfr_t Cx, const mpfr_t Cy);
  virtual ~reference();
  complex<floatexp> get() const;
  bool escaped() const;
  virtual void step() = 0;
};

struct formula
{
  formula() { }
  virtual ~formula() { }
  virtual bool complex_analytic() const noexcept = 0;
  virtual std::string name() const noexcept = 0;
  virtual reference *new_reference(const mpfr_t Cx, const mpfr_t Cy) const = 0;
};

struct formulaC : public formula
{
  formulaC() { }
  virtual ~formulaC() { }
  virtual bool complex_analytic() const noexcept { return true; }
  virtual std::string name() const noexcept = 0;
  virtual reference *new_reference(const mpfr_t Cx, const mpfr_t Cy) const = 0;

  virtual blaC<float> bla1(const float h, const float k, const float L, const complex<float> Z) const noexcept = 0;
  virtual blaC<double> bla1(const double h, const double k, const double L, const complex<double> Z) const noexcept = 0;
  virtual blaC<long double> bla1(const long double h, const long double k, const long double L, const complex<long double> Z) const noexcept = 0;
  virtual blaC<floatexp> bla1(const floatexp h, const floatexp k, const floatexp L, const complex<floatexp> Z) const noexcept = 0;

  virtual complex<float> perturb(const complex<float> &C, const complex<float> &Z, const complex<float> &c, const complex<float> &z) const noexcept = 0;
  virtual complex<double> perturb(const complex<double> &C, const complex<double> &Z, const complex<double> &c, const complex<double> &z) const noexcept = 0;
  virtual complex<long double> perturb(const complex<long double> &C, const complex<long double> &Z, const complex<long double> &c, const complex<long double> &z) const noexcept = 0;
  virtual complex<floatexp> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const complex<floatexp> &c, const complex<floatexp> &z) const noexcept = 0;

  virtual dual<1, complex<float>> perturb(const complex<float> &C, const complex<float> &Z, const dual<1, complex<float>> &c, const dual<1, complex<float>> &z) const noexcept = 0;
  virtual dual<1, complex<double>> perturb(const complex<double> &C, const complex<double> &Z, const dual<1, complex<double>> &c, const dual<1, complex<double>> &z) const noexcept = 0;
  virtual dual<1, complex<long double>> perturb(const complex<long double> &C, const complex<long double> &Z, const dual<1, complex<long double>> &c, const dual<1, complex<long double>> &z) const noexcept = 0;
  virtual dual<1, complex<floatexp>> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const dual<1, complex<floatexp>> &c, const dual<1, complex<floatexp>> &z) const noexcept = 0;
};

struct formulaR2 : public formula
{
  formulaR2() { }
  virtual ~formulaR2() { }
  virtual bool complex_analytic() const noexcept { return false; }
  virtual std::string name() const noexcept = 0;
  virtual reference *new_reference(const mpfr_t Cx, const mpfr_t Cy) const = 0;

  virtual blaR2<float> bla1(const float h, const float k, const float L, const complex<float> Z) const noexcept = 0;
  virtual blaR2<double> bla1(const double h, const double k, const double L, const complex<double> Z) const noexcept = 0;
  virtual blaR2<long double> bla1(const long double h, const long double k, const long double L, const complex<long double> Z) const noexcept = 0;
  virtual blaR2<floatexp> bla1(const floatexp h, const floatexp k, const floatexp L, const complex<floatexp> Z) const noexcept = 0;

  virtual complex<float> perturb(const complex<float> &C, const complex<float> &Z, const complex<float> &c, const complex<float> &z) const noexcept = 0;
  virtual complex<double> perturb(const complex<double> &C, const complex<double> &Z, const complex<double> &c, const complex<double> &z) const noexcept = 0;
  virtual complex<long double> perturb(const complex<long double> &C, const complex<long double> &Z, const complex<long double> &c, const complex<long double> &z) const noexcept = 0;
  virtual complex<floatexp> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const complex<floatexp> &c, const complex<floatexp> &z) const noexcept = 0;

  virtual complex<dual<2, float>> perturb(const complex<float> &C, const complex<float> &Z, const complex<dual<2, float>> &c, const complex<dual<2, float>> &z) const noexcept = 0;
  virtual complex<dual<2, double>> perturb(const complex<double> &C, const complex<double> &Z, const complex<dual<2, double>> &c, const complex<dual<2, double>> &z) const noexcept = 0;
  virtual complex<dual<2, long double>> perturb(const complex<long double> &C, const complex<long double> &Z, const complex<dual<2, long double>> &c, const complex<dual<2, long double>> &z) const noexcept = 0;
  virtual complex<dual<2, floatexp>> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const complex<dual<2, floatexp>> &c, const complex<dual<2, floatexp>> &z) const noexcept = 0;
};

extern std::vector<formula *> formulas;
void formulas_init();
