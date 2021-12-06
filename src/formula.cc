// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "complex.h"
#include "floatexp.h"
#include "formula.h"

reference::reference(const mpfr_t Cx0, const mpfr_t Cy0)
{
  mpfr_prec_t prec = std::max(mpfr_get_prec(Cx0), mpfr_get_prec(Cy0));
  mpfr_init2(Cx, prec); mpfr_set(Cx, Cx0, MPFR_RNDN);
  mpfr_init2(Cy, prec); mpfr_set(Cy, Cy0, MPFR_RNDN);
  mpfr_init2(Zx, prec); mpfr_set_ui(Zx, 0, MPFR_RNDN);
  mpfr_init2(Zy, prec); mpfr_set_ui(Zy, 0, MPFR_RNDN);
}

reference::~reference()
{
  mpfr_clear(Cx);
  mpfr_clear(Cy);
  mpfr_clear(Zx);
  mpfr_clear(Zy);
}

complex<floatexp> reference::get() const
{
  long ex = 0, ey = 0;
  double x = mpfr_get_d_2exp(&ex, Zx, MPFR_RNDN);
  double y = mpfr_get_d_2exp(&ey, Zy, MPFR_RNDN);
  return complex<floatexp>(floatexp(x, ex), floatexp(y, ey));
}

bool reference::escaped() const
{
  return 4 < norm(get()); // FIXME hardcoded
}

#include "formula_mandelbrot.h"
#include "formula_burningship.h"

std::vector<formula *> formulas;

void formulas_init()
{
  formulas.push_back(new formulaC_mandelbrot());
  formulas.push_back(new formulaR2_burningship());
}
