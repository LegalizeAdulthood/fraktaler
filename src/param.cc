// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <iostream>

#include "map.h"
#include "param.h"

param::param()
: p
  { { "0", "0", "1", 0 }
  , { 1024, 1024, 1024, 625.0 }
  , { true, false, false, std::vector<std::string>{ "float", "double", "long double", "softfloat"} }
  , { 1024, 576, 1, 1 }
  , { false, 0, 0, 0, false }
  , { "fraktaler-3.exr", false, 0, 0 }
  }
, center(0)
, zoom(1)
{
  home(*this);
  restring(*this);
}

void restring(param &par)
{
  par.p.location.real = par.center.x.toString();
  par.p.location.imag = par.center.y.toString();
  { std::ostringstream s; s << par.zoom; par.p.location.zoom = s.str(); }
  { std::ostringstream s; s << par.p.bailout.iterations; par.s_iterations = s.str(); }
  { std::ostringstream s; s << par.p.bailout.maximum_reference_iterations; par.s_maximum_reference_iterations = s.str(); }
  { std::ostringstream s; s << par.p.bailout.maximum_perturb_iterations; par.s_maximum_perturb_iterations = s.str(); }
  { std::ostringstream s; s << par.p.bailout.escape_radius; par.s_escape_radius = s.str(); }
  par.transform = mat2<double>(polar2<double>(par.p.transform.reflect ? -1 : 1, 1, par.p.transform.rotate, std::exp2(par.p.transform.stretch_amount / 100), par.p.transform.stretch_angle));
}

void unstring(param &par)
{
  par.zoom = floatexp(par.p.location.zoom);
  mpfr_prec_t prec = std::min(24, 24 + (par.zoom * par.p.image.height).exp);
  par.center.x.set_prec(prec);
  par.center.y.set_prec(prec);
  par.center.x = par.p.location.real;
  par.center.y = par.p.location.imag;
  par.p.bailout.iterations = std::atoll(par.s_iterations.c_str());
  par.p.bailout.maximum_reference_iterations = std::atoll(par.s_maximum_reference_iterations.c_str());
  par.p.bailout.maximum_perturb_iterations = std::atoll(par.s_maximum_perturb_iterations.c_str());
  par.p.bailout.escape_radius = std::atof(par.s_escape_radius.c_str());
  polar2<double> P(par.transform);
  par.p.transform.reflect = P.sign < 0 ? true : false;
  par.p.transform.rotate = P.rotate;
  par.p.transform.stretch_amount = 100 * std::log2(P.stretch_factor);
  par.p.transform.stretch_angle = P.stretch_angle;
}

void home(param &par)
{
  par.center.x.set_prec(24);
  par.center.y.set_prec(24);
  par.center = 0;
  par.zoom = 1;
  par.p.bailout.iterations = 1024;
  par.p.bailout.maximum_reference_iterations = par.p.bailout.iterations;
  par.p.bailout.maximum_perturb_iterations = 1024;
  par.p.location.period = 1;
  par.p.algorithm.lock_maximum_reference_iterations_to_period = true;
  par.p.algorithm.reuse_reference = true;
  par.p.algorithm.reuse_bilinear_approximation = true;
  par.transform = mat2<double>(1);
  restring(par);
}

void zoom(param &par, double x, double y, double g, bool fixed_click)
{
  complex<double> w (x * par.p.image.width / par.p.image.height, -y);
  w = par.transform * w;
  floatexp d = (fixed_click ? 1 - 1 / g : 1) * 2 / par.zoom;
  floatexp u = d * w.x;
  floatexp v = d * w.y;
  par.zoom *= g;
  mpfr_prec_t prec = std::max(24, 24 + (par.zoom * par.p.image.height).exp);
  mpfr_t dx, dy;
  mpfr_init2(dx, 53);
  mpfr_init2(dy, 53);
  mpfr_set_d(dx, u.val, MPFR_RNDN);
  mpfr_set_d(dy, v.val, MPFR_RNDN);
  mpfr_mul_2si(dx, dx, u.exp, MPFR_RNDN);
  mpfr_mul_2si(dy, dy, v.exp, MPFR_RNDN);
  par.center.x.set_prec(prec);
  par.center.y.set_prec(prec);
  mpfr_add(par.center.x.mpfr_ptr(), par.center.x.mpfr_srcptr(), dx, MPFR_RNDN);
  mpfr_add(par.center.y.mpfr_ptr(), par.center.y.mpfr_srcptr(), dy, MPFR_RNDN);
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
  par.zoom *= g;
  mpfr_prec_t prec = std::max(24, 24 + (par.zoom * par.p.image.height).exp);
  par.center.x.set_prec(prec);
  par.center.y.set_prec(prec);
  // rotate, skew
  T2 /= g;
  T2 = inverse(T2);
  par.transform = par.transform * T2;
  // done
  restring(par);
}

void param::parse(std::istream &in, const std::string &filename)
{
  // throws
  p = toml::get<pparam>(toml::parse(in, filename));
  restring(*this);
}
