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

struct formulaC
{
  formulaC();
  virtual ~formulaC();
  virtual std::string name() const = 0;
  virtual reference *new_reference(const mpfr_t Cx, const mpfr_t Cy) const = 0;

  virtual bla<float> bla1(const float h, const float k, const float L, const complex<float> Z) const = 0;
  virtual bla<double> bla1(const double h, const double k, const double L, const complex<double> Z) const = 0;
  virtual bla<long double> bla1(const long double h, const long double k, const long double L, const complex<long double> Z) const = 0;
  virtual bla<floatexp> bla1(const floatexp h, const floatexp k, const floatexp L, const complex<floatexp> Z) const = 0;

  virtual complex<float> perturb(const complex<float> &C, const complex<float> &Z, const complex<float> &c, const complex<float> &z) const = 0;
  virtual complex<double> perturb(const complex<double> &C, const complex<double> &Z, const complex<double> &c, const complex<double> &z) const = 0;
  virtual complex<long double> perturb(const complex<long double> &C, const complex<long double> &Z, const complex<long double> &c, const complex<long double> &z) const = 0;
  virtual complex<floatexp> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const complex<floatexp> &c, const complex<floatexp> &z) const = 0;

  virtual dual<1, complex<float>> perturb(const complex<float> &C, const complex<float> &Z, const dual<1, complex<float>> &c, const dual<1, complex<float>> &z) const = 0;
  virtual dual<1, complex<double>> perturb(const complex<double> &C, const complex<double> &Z, const dual<1, complex<double>> &c, const dual<1, complex<double>> &z) const = 0;
  virtual dual<1, complex<long double>> perturb(const complex<long double> &C, const complex<long double> &Z, const dual<1, complex<long double>> &c, const dual<1, complex<long double>> &z) const = 0;
  virtual dual<1, complex<floatexp>> perturb(const complex<floatexp> &C, const complex<floatexp> &Z, const dual<1, complex<floatexp>> &c, const dual<1, complex<floatexp>> &z) const = 0;
};

extern std::vector<formulaC *> formulas;
void formulas_init();
