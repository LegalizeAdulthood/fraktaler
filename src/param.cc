// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <iostream>
#include <toml.hpp>

#include "colour.h"
#include "formula.h"
#include "map.h"
#include "param.h"
#include "source.h"

param::param()
: center(0)
, zoom(1)
, reference(0)
{
  home(*this);
}

void restring_locs(param &par)
{
  par.p.location.real = par.center.x.toString();
  par.p.location.imag = par.center.y.toString();
  par.p.reference.real = par.reference.x.toString();
  par.p.reference.imag = par.reference.y.toString();
}

void restring_vals(param &par)
{
  { std::ostringstream s; s << par.p.reference.period; par.s_period = s.str(); }
  { std::ostringstream s; s << par.zoom; par.p.location.zoom = s.str(); }
  { std::ostringstream s; s << par.p.bailout.iterations; par.s_iterations = s.str(); }
  { std::ostringstream s; s << par.p.bailout.maximum_reference_iterations; par.s_maximum_reference_iterations = s.str(); }
  { std::ostringstream s; s << par.p.bailout.maximum_perturb_iterations; par.s_maximum_perturb_iterations = s.str(); }
  { std::ostringstream s; s << par.p.bailout.escape_radius; par.s_escape_radius = s.str(); }
  par.transform = mat2<double>(polar2<double>(par.p.transform.reflect ? -1 : 1, 1, par.p.transform.rotate * 2 * M_PI / 360, std::exp2(par.p.transform.stretch_amount / 100), par.p.transform.stretch_angle * 2 * M_PI / 360));
}

void unstring_locs(param &par)
{
  par.zoom = floatexp(par.p.location.zoom);
  mpfr_prec_t prec = std::max(24, 24 + (par.zoom * par.p.image.height).exp);
  par.center.x.set_prec(prec);
  par.center.y.set_prec(prec);
  par.center.x = par.p.location.real;
  par.center.y = par.p.location.imag;
  par.reference.x.set_prec(prec); // FIXME
  par.reference.y.set_prec(prec); // FIXME
  par.reference.x = par.p.reference.real;
  par.reference.y = par.p.reference.imag;
}

void unstring_vals(param &par)
{
  par.p.reference.period = std::atoll(par.s_period.c_str());
  par.p.bailout.iterations = std::atoll(par.s_iterations.c_str());
  par.p.bailout.maximum_reference_iterations = std::atoll(par.s_maximum_reference_iterations.c_str());
  par.p.bailout.maximum_perturb_iterations = std::atoll(par.s_maximum_perturb_iterations.c_str());
  par.p.bailout.escape_radius = std::atof(par.s_escape_radius.c_str());
  polar2<double> P(par.transform);
  par.p.transform.reflect = P.sign < 0 ? true : false;
  par.p.transform.rotate = P.rotate * 360 / (2 * M_PI);
  par.p.transform.stretch_amount = 100 * std::log2(P.stretch_factor);
  par.p.transform.stretch_angle = P.stretch_angle * 360 / (2 * M_PI);
}

void home(param &par)
{
  par.reference.x.set_prec(24);
  par.reference.y.set_prec(24);
  par.center.x.set_prec(24);
  par.center.y.set_prec(24);
  par.center = 0;
  par.zoom = 1;
  par.p.bailout.iterations = 1024;
  par.p.bailout.maximum_reference_iterations = par.p.bailout.iterations;
  par.p.bailout.maximum_perturb_iterations = 1024;
  par.p.reference.period = 0;
  par.p.algorithm.lock_maximum_reference_iterations_to_period = false;
  par.p.algorithm.reuse_reference = false;
  par.p.algorithm.reuse_bilinear_approximation = false;
  restring_locs(par);
  restring_vals(par);
  par.transform = mat2<double>(1);
  unstring_vals(par);
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
  restring_locs(par);
  restring_vals(par);
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
  // done
  restring_locs(par);
  restring_vals(par);
  par.transform = par.transform * T2;
  unstring_vals(par);
}

void param::load_toml(const std::string &filename)
{
  std::ifstream ifs(filename, std::ios_base::binary);
  ifs.exceptions(std::ifstream::badbit);
  ifs >> p;
  unstring_locs(*this);
  restring_vals(*this);
}

void param::from_string(const std::string &str)
{
  std::istringstream ifs(str);
  ifs >> p;
  unstring_locs(*this);
  restring_vals(*this);
}

std::istream &operator>>(std::istream &ifs, pparam &p)
{
  auto t = toml::parse(ifs);
  std::string formula_name = toml::find_or(t, "formula", "name", formulas[p.formula_id]->name());
  size_t formula_id = 0;
  for (const auto &f : formulas)
  {
    if (f->name() == formula_name)
    {
      break;
    }
    formula_id++;
  }
  if (formula_id < formulas.size())
  {
    p.formula_id = formula_id;
  }
  else
  {
    // FIXME formula not found
  }
  std::string colour_name = toml::find_or(t, "colour", "name", colours[p.colour_id]->name());
  size_t colour_id = 0;
  for (const auto &c : colours)
  {
    if (c->name() == colour_name)
    {
      break;
    }
    colour_id++;
  }
  if (colour_id < colours.size())
  {
    p.colour_id = colour_id;
  }
  else
  {
    // FIXME colour not found
  }
#define LOAD(a,b) p.a.b = toml::find_or(t, #a, #b, p.a.b);
  LOAD(location, real)
  LOAD(location, imag)
  LOAD(location, zoom)
  p.reference.real = p.location.real; // FIXME
  p.reference.imag = p.location.imag; // FIXME
  p.reference.period = 0; // FIXME
  LOAD(reference, real)
  LOAD(reference, imag)
  LOAD(reference, period)
  LOAD(bailout, iterations)
  LOAD(bailout, maximum_reference_iterations)
  LOAD(bailout, maximum_perturb_iterations)
  LOAD(bailout, escape_radius)
  LOAD(algorithm, lock_maximum_reference_iterations_to_period)
  LOAD(algorithm, reuse_reference)
  LOAD(algorithm, reuse_bilinear_approximation)
  LOAD(algorithm, number_types)
  LOAD(image, width)
  LOAD(image, height)
  LOAD(image, subsampling)
  LOAD(image, subframes)
  LOAD(transform, reflect)
  LOAD(transform, rotate)
  LOAD(transform, stretch_angle)
  LOAD(transform, stretch_amount)
  LOAD(transform, exponential_map)
  LOAD(render, filename)
  LOAD(render, zoom_out_sequence);
  LOAD(render, start_frame);
  LOAD(render, frame_count);
#undef LOAD
  return ifs;
}

std::ostream &operator<<(std::ostream &ofs, const pparam &p)
{
  pparam q;
  ofs << "program = " << toml::value("fraktaler-3") << "\n";
  ofs << "version = " << toml::value(fraktaler_3_version_string) << "\n";
  if (p.formula_id != q.formula_id)
  {
    ofs << "formula.name = " << toml::value(formulas[p.formula_id]->name()) << "\n";
  }
  if (p.colour_id != q.colour_id)
  {
    ofs << "colour.name = " << toml::value(colours[p.colour_id]->name()) << "\n";
  }
#define SAVE(a,b) if (p.a.b != q.a.b) { ofs << #a << "." << #b << " = " << std::setw(70) << toml::value(p.a.b) << "\n"; }
  SAVE(location, real)
  SAVE(location, imag)
  SAVE(location, zoom)
  if (p.location.real != p.reference.real) SAVE(reference, real)
  if (p.location.imag != p.reference.imag) SAVE(reference, real)
  SAVE(reference, period)
  SAVE(bailout, iterations)
  SAVE(bailout, maximum_reference_iterations)
  SAVE(bailout, maximum_perturb_iterations)
  SAVE(bailout, escape_radius)
  SAVE(algorithm, lock_maximum_reference_iterations_to_period)
  SAVE(algorithm, reuse_reference)
  SAVE(algorithm, reuse_bilinear_approximation)
  SAVE(algorithm, number_types)
  SAVE(image, width)
  SAVE(image, height)
  SAVE(image, subsampling)
  SAVE(image, subframes)
  SAVE(transform, reflect)
  SAVE(transform, rotate)
  SAVE(transform, stretch_angle)
  SAVE(transform, stretch_amount)
  SAVE(transform, exponential_map)
  SAVE(render, filename)
  SAVE(render, zoom_out_sequence);
  SAVE(render, start_frame);
  SAVE(render, frame_count);
#undef SAVE
  return ofs;
}

void param::save_toml(const std::string &filename) const
{
  std::ofstream ofs(filename, std::ios_base::binary);
  ofs.exceptions(std::ofstream::badbit);
  ofs << p;
}

std::string param::to_string() const
{
  std::ostringstream ofs;
  ofs << p;
  return ofs.str();
}
