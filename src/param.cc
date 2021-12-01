// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <iostream>

#include "param.h"

void zoom(param &par, double x, double y, double g)
{
  floatexp d = (1 - 1 / g) * 2 / par.Zoom;
  floatexp u = d * x * par.Width / par.Height;
  floatexp v = d * y;
  par.Zoom *= g;
  mpfr_prec_t prec = std::max(24, 24 + (par.Zoom * par.Height).exp);
  mpfr_t dx, dy;
  mpfr_init2(dx, 53);
  mpfr_init2(dy, 53);
  mpfr_set_d(dx, u.val, MPFR_RNDN);
  mpfr_set_d(dy, v.val, MPFR_RNDN);
  mpfr_mul_2si(dx, dx, u.exp, MPFR_RNDN);
  mpfr_mul_2si(dy, dy, v.exp, MPFR_RNDN);
  mpfr_prec_round(par.Cx, prec, MPFR_RNDN);
  mpfr_prec_round(par.Cy, prec, MPFR_RNDN);
  mpfr_add(par.Cx, par.Cx, dx, MPFR_RNDN);
  mpfr_add(par.Cy, par.Cy, dy, MPFR_RNDN);
  mpfr_clear(dx);
  mpfr_clear(dy);
  mpfr_fprintf(stderr, "Re: %Re\nIm: %Re\n", par.Cx, par.Cy);
  std::cerr << "Zoom: " << par.Zoom << std::endl;
}
