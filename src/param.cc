// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <iostream>

#include "param.h"

void restring(param &par)
{
  int re_bytes = mpfr_snprintf(nullptr, 0, "%Re", par.Cx);
  if (re_bytes >= 0)
  {
    char *sRe = new char[re_bytes + 1];
    mpfr_snprintf(sRe, re_bytes + 1, "%Re", par.Cx);
    par.sRe = sRe;
    delete[] sRe;
  }
  else
  {
    par.sRe = "0";
  }
  int im_bytes = mpfr_snprintf(nullptr, 0, "%Re", par.Cy);
  if (im_bytes >= 0)
  {
    char *sIm = new char[im_bytes + 1];
    mpfr_snprintf(sIm, im_bytes + 1, "%Re", par.Cy);
    par.sIm = sIm;
    delete[] sIm;
  }
  else
  {
    par.sIm = "0";
  }
  { std::ostringstream s; s << par.Zoom; par.sZoom = s.str(); }
  { std::ostringstream s; s << par.Iterations; par.sIterations = s.str(); }
  { std::ostringstream s; s << par.MaxRefIters; par.sMaxRefIters = s.str(); }
  { std::ostringstream s; s << par.MaxPtbIters; par.sMaxPtbIters = s.str(); }
  { std::ostringstream s; s << par.EscapeRadius; par.sEscapeRadius = s.str(); }
}

void home(param &par)
{
  mpfr_set_prec(par.Cx, 24);
  mpfr_set_prec(par.Cy, 24);
  mpfr_set_d(par.Cx, 0, MPFR_RNDN);
  mpfr_set_d(par.Cy, 0, MPFR_RNDN);
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
  mpfr_prec_round(par.Cx, prec, MPFR_RNDN);
  mpfr_prec_round(par.Cy, prec, MPFR_RNDN);
  mpfr_add(par.Cx, par.Cx, dx, MPFR_RNDN);
  mpfr_add(par.Cy, par.Cy, dy, MPFR_RNDN);
  mpfr_clear(dx);
  mpfr_clear(dy);
  restring(par);
}
