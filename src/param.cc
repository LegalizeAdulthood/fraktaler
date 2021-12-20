// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <iostream>

#include "map.h"
#include "param.h"

param::param()
: C(0)
, Zoom(1)
, Iterations(1024)
, MaxRefIters(1024)
, MaxPtbIters(1024)
, EscapeRadius(625)
, ReferencePeriod(0)
, LockMaxRefItersToPeriod(false)
, ReuseReference(false)
, ReuseBLA(false)
, MaxSubframes(1)
, ExponentialMap(false)
, ZoomOutSequence(false)
, Channels(Channels_default)
, Stem("fraktaler-3.exr")
, Width(1920)
, Height(1080)
, K(1)
, ForceNumberType(nt_none)
{
  home(*this);
  restring(*this);
}


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
  par.K = mat2<double>(1);
  restring(par);
}

void zoom(param &par, double x, double y, double g, bool fixed_click)
{
  complex<double> w (x * par.Width / par.Height, -y);
  w = par.K * w;
  floatexp d = (fixed_click ? 1 - 1 / g : 1) * 2 / par.Zoom;
  floatexp u = d * w.x;
  floatexp v = d * w.y;
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

void zoom(param &par, const mat3 &T, const mat3 &T0)
{
  // translate
  vec3 t = T * vec3(0.0f, 0.0f, 1.0f);
  vec2 w = vec2(t) / t.z;
  zoom(par, w.x, -w.y, 1, false);
  // zoom
  mat2<double> T2(T0);
  double g = std::sqrt(std::abs(determinant(T2)));
  par.Zoom *= g;
  mpfr_prec_t prec = std::max(24, 24 + (par.Zoom * par.Height).exp);
  par.C.x.set_prec(prec);
  par.C.y.set_prec(prec);
  // rotate, skew
  T2 /= g;
  T2 = inverse(T2);
  par.K = par.K * T2;
  // done
  restring(par);
}
