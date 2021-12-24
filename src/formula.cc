// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "formula.h"
#include "formula_mandelbrot.h"
#include "formula_burningship.h"

std::vector<formula *> formulas;

void formulas_init()
{
  #define C(stem) formulas.push_back(new formulaC< \
    stem ## _name, \
    stem ## _plain <complex<mpreal>>, \
    stem ## _plain <dual<1, complex<mpreal>>>, \
    stem ## _perturb <complex<float>, dual<1, complex<float>>>, \
    stem ## _perturb <complex<double>, dual<1, complex<double>>>, \
    stem ## _perturb <complex<long double>, dual<1, complex<long double>>>, \
    stem ## _perturb <complex<floatexp>, dual<1, complex<floatexp>>>, \
    stem ## _perturb <complex<softfloat>, dual<1, complex<softfloat>>>, \
    stem ## _perturb <complex<float128>, dual<1, complex<float128>>>, \
    stem ## _bla <float>, \
    stem ## _bla <double>, \
    stem ## _bla <long double>, \
    stem ## _bla <floatexp>, \
    stem ## _bla <softfloat>, \
    stem ## _bla <float128> \
  >());
  #define R2(stem) formulas.push_back(new formulaR2< \
    stem ## _name, \
    stem ## _plain <mpreal>, \
    stem ## _plain <dual<2, mpreal>>, \
    stem ## _perturb <float, dual<2, float>>, \
    stem ## _perturb <double, dual<2, double>>, \
    stem ## _perturb <long double, dual<2, long double>>, \
    stem ## _perturb <floatexp, dual<2, floatexp>>, \
    stem ## _perturb <softfloat, dual<2, softfloat>>, \
    stem ## _perturb <float128, dual<2, float128>>, \
    stem ## _bla <float>, \
    stem ## _bla <double>, \
    stem ## _bla <long double>, \
    stem ## _bla <floatexp>, \
    stem ## _bla <softfloat>, \
    stem ## _bla <float128> \
  >());
  C(mandelbrot)
  R2(burningship)
}
