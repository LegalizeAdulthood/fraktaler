// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <iostream>

#include "param.h"

void restring(param &par)
{
  par.sRe = par.C.x.toString();
  par.sIm = par.C.y.toString();
  { std::ostringstream s; s << par.Zoom; par.sZoom = s.str(); }
  { std::ostringstream s; s << par.Iterations; par.sIterations = s.str(); }
  { std::ostringstream s; s << par.MaxRefIters; par.sMaxRefIters = s.str(); }
  { std::ostringstream s; s << par.MaxPtbIters; par.sMaxPtbIters = s.str(); }
  { std::ostringstream s; s << par.EscapeRadius; par.sEscapeRadius = s.str(); }
}

void home(param &par)
{
  par.C.x.set_prec(24);
  par.C.y.set_prec(24);
  par.C = 0;
  par.Zoom = 1;
  par.Iterations = 1024;
  par.MaxRefIters = par.Iterations;
  par.MaxPtbIters = 1024;
  par.ReferencePeriod = 0;
  par.LockMaxRefItersToPeriod = false;
  par.ReuseReference = false;
  par.ReuseBLA = false;
  restring(par);
}

void zoom(param &par, double x, double y, double g, bool fixed_click)
{
  floatexp d = (fixed_click ? 1 - 1 / g : 1) * 2 / par.Zoom;
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
  par.C.x.set_prec(prec);
  par.C.y.set_prec(prec);
  mpfr_add(par.C.x.mpfr_ptr(), par.C.x.mpfr_srcptr(), dx, MPFR_RNDN);
  mpfr_add(par.C.y.mpfr_ptr(), par.C.y.mpfr_srcptr(), dy, MPFR_RNDN);
  mpfr_clear(dx);
  mpfr_clear(dy);
  restring(par);
}
