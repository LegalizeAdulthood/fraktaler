// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <iostream>
#include <map>

#include <toml.hpp>

#include "param.h"
#include "source.h"

param::param()
: center(0)
, zoom(1)
, reference(0)
{
  post_edit_formula(*this);
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
  { std::ostringstream s; s << par.p.bailout.maximum_bla_steps; par.s_maximum_bla_steps = s.str(); }
  { std::ostringstream s; s << par.p.bailout.escape_radius; par.s_escape_radius = s.str(); }
  { std::ostringstream s; s << par.p.bailout.inscape_radius; par.s_inscape_radius = s.str(); }
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
  par.p.bailout.maximum_bla_steps = std::atoll(par.s_maximum_bla_steps.c_str());
  par.p.bailout.escape_radius = std::atof(par.s_escape_radius.c_str());
  par.p.bailout.inscape_radius = std::atof(par.s_inscape_radius.c_str());
  polar2<double> P(par.transform);
  par.p.transform.reflect = P.sign < 0 ? true : false;
  par.p.transform.rotate = P.rotate * 360 / (2 * M_PI);
  par.p.transform.stretch_amount = 100 * std::log2(P.stretch_factor);
  par.p.transform.stretch_angle = P.stretch_angle * 360 / (2 * M_PI);
}

void post_edit_formula(param &par)
{
  par.opss = compile_formula(par.p.formula);
  assert(validate_opcodess(par.opss));
  par.degrees = opcodes_degrees(par.opss);
  par.s_opss = print_opcodess(par.opss);
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
  par.p.bailout.maximum_bla_steps = 1024;
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

complex<floatexp> get_delta_c(const param &par, double x, double y)
{
  complex<double> w (x * par.p.image.width / par.p.image.height, -y);
  w = par.transform * w;
  return complex<floatexp>(floatexp(w.x / par.zoom), floatexp(w.y / par.zoom));
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
  post_edit_formula(*this);
}

void param::from_string(const std::string &str)
{
  std::istringstream ifs(str);
  ifs >> p;
  unstring_locs(*this);
  restring_vals(*this);
  post_edit_formula(*this);
}

std::istream &operator>>(std::istream &ifs, pparam &p)
{
  auto t = toml::parse(ifs);
  auto f = toml::find_or(t, "formula", std::vector<toml::table>());
  std::vector<phybrid1> per;
  for (auto f1 : f)
  {
    toml::value g(f1);
    phybrid1 h;
    h.abs_x = toml::find_or(g, "abs_x", h.abs_x);
    h.abs_y = toml::find_or(g, "abs_y", h.abs_y);
    h.neg_x = toml::find_or(g, "neg_x", h.neg_x);
    h.neg_y = toml::find_or(g, "neg_y", h.neg_y);
    h.power = toml::find_or(g, "power", h.power);
    h.opcodes = parse_opcodes(toml::find_or(g, "opcodes", print_opcodes(h.opcodes)));
    per.push_back(h);
  }
  if (! per.empty())
  {
    p.formula.per = per;
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
  LOAD(bailout, maximum_bla_steps)
  LOAD(bailout, escape_radius)
  LOAD(bailout, inscape_radius)
  LOAD(algorithm, lock_maximum_reference_iterations_to_period)
  LOAD(algorithm, reuse_reference)
  LOAD(algorithm, reuse_bilinear_approximation)
  LOAD(algorithm, bla_skip_levels)
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
  LOAD(render, zoom_out_factor);
  LOAD(render, start_frame);
  LOAD(render, frame_count);
  LOAD(newton, action);
  LOAD(newton, domain);
  LOAD(newton, absolute);
  LOAD(newton, power);
  LOAD(newton, factor);
  LOAD(opencl, platform);
  LOAD(opencl, device);
  LOAD(opencl, tile_width);
  LOAD(opencl, tile_height);
#undef LOAD
  return ifs;
}

std::ostream &operator<<(std::ostream &ofs, const pparam &p)
{
  pparam q;
  ofs << "program = " << toml::value("fraktaler-3") << "\n";
  ofs << "version = " << toml::value(fraktaler_3_version_string) << "\n";
#define SAVE(a,b) if (p.a.b != q.a.b) { ofs << #a << "." << #b << " = " << std::setw(70) << toml::value(p.a.b) << "\n"; }
  SAVE(location, real)
  SAVE(location, imag)
  SAVE(location, zoom)
  if (p.location.real != p.reference.real) SAVE(reference, real)
  if (p.location.imag != p.reference.imag) SAVE(reference, imag)
  SAVE(reference, period)
  SAVE(bailout, iterations)
  SAVE(bailout, maximum_reference_iterations)
  SAVE(bailout, maximum_perturb_iterations)
  SAVE(bailout, maximum_bla_steps)
  SAVE(bailout, escape_radius)
  SAVE(bailout, inscape_radius)
  SAVE(algorithm, lock_maximum_reference_iterations_to_period)
  SAVE(algorithm, reuse_reference)
  SAVE(algorithm, reuse_bilinear_approximation)
  SAVE(algorithm, bla_skip_levels)
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
  SAVE(render, zoom_out_factor);
  SAVE(render, start_frame);
  SAVE(render, frame_count);
  SAVE(newton, action);
  SAVE(newton, domain);
  SAVE(newton, absolute);
  SAVE(newton, power);
  SAVE(newton, factor);
  SAVE(opencl, platform);
  SAVE(opencl, device);
  SAVE(opencl, tile_width);
  SAVE(opencl, tile_height);
#undef SAVE
  if (! (p.formula == q.formula))
  {
    toml::array per;
    phybrid1 def;
    for (auto h : p.formula.per)
    {
      std::map<std::string, toml::value> f;
      if (h.opcodes.size())
      {
        f["opcodes"] = print_opcodes(h.opcodes);
      }
      else
      {
        if (h.abs_x != def.abs_x) f["abs_x"] = h.abs_x;
        if (h.abs_y != def.abs_y) f["abs_y"] = h.abs_y;
        if (h.neg_x != def.neg_x) f["neg_x"] = h.neg_x;
        if (h.neg_y != def.neg_y) f["neg_y"] = h.neg_y;
        if (h.power != def.power) f["power"] = h.power;
      }
      per.push_back(f);
    }
    std::map<std::string, toml::array> f;
    f["formula"] = per;
    ofs << toml::value(f) << "\n";
  }
  return ofs;
}

void param::save_toml(const std::string &filename) const
{
  std::ofstream ofs;
  ofs.exceptions(std::ofstream::badbit | std::ofstream::failbit);
  ofs.open(filename, std::ios_base::binary);
  ofs << p;
}

std::string param::to_string() const
{
  std::ostringstream ofs;
  ofs << p;
  return ofs.str();
}

std::vector<std::vector<opcode>> compile_formula(const phybrid &H)
{
  std::vector<std::vector<opcode>> result;
  for (const auto & h : H.per)
  {
    std::vector<opcode> current;
    if (h.opcodes.size())
    {
      bool need_store = false;
      for (const auto & op : h.opcodes)
      {
        if (op == op_mul)
        {
          need_store = true;
          break;
        }
        if (op == op_store)
        {
          break;
        }
      }
      if (need_store)
      {
        current.push_back(op_store);
      }
      current.insert(current.end(), h.opcodes.begin(), h.opcodes.end());
    }
    else
    {
      if (h.abs_x) current.push_back(op_absx);
      if (h.abs_y) current.push_back(op_absy);
      if (h.neg_x) current.push_back(op_negx);
      if (h.neg_y) current.push_back(op_negy);
      int p = h.power;
      bool is_power_of_two = (p & (p - 1)) == 0;
      if (! is_power_of_two)
      {
        current.push_back(op_store);
      }
      std::vector<opcode> power;
      while (p > 1)
      {
        if (p & 1)
        {
          assert(! is_power_of_two);
          power.push_back(op_mul);
        }
        power.push_back(op_sqr);
        p >>= 1;
      }
      std::reverse(power.begin(), power.end());
      // FIXME remove this validation code later
      int q = 1;
      for (const auto & op : power)
      {
        switch (op)
        {
          case op_mul: q += 1; break;
          case op_sqr: q <<= 1; break;
          case op_add:
          case op_store:
          case op_absx:
          case op_absy:
          case op_negx:
          case op_negy:
            assert(! "expected mul or sqr");
            break;
        }
      }
      assert(q == h.power);
      current.insert(current.end(), power.begin(), power.end());
      current.push_back(op_add);
    }
    result.push_back(current);
  }
  return result;
}

int opcodes_degree(const std::vector<opcode> &ops)
{
  int deg_stored = 0;
  int deg = 1;
  for (const auto & op : ops)
  {
    switch (op)
    {
      case op_store: deg_stored = deg; break;
      case op_mul: deg += deg_stored; break;
      case op_sqr: deg <<= 1; break;
      // remainder have no effect
      case op_add:
      case op_absx:
      case op_absy:
      case op_negx:
      case op_negy:
        break;
    }
  }
  return deg;
}

std::string print_opcodes(const std::vector<opcode> &ops)
{
  std::ostringstream s;
  bool first_word = true;
  for (const auto & op : ops)
  {
    if (! first_word)
    {
      s << " ";
    }
    first_word = false;
    s << op_string[op];
  }
  return s.str();
}

std::string print_opcodess(const std::vector<std::vector<opcode>> &opss)
{
  std::ostringstream s;
  bool first_line = true;
  for (const auto & ops : opss)
  {
    if (! first_line)
    {
      s << "\n";
    }
    first_line = false;
    s << print_opcodes(ops);
  }
  return s.str();
}

std::vector<int> opcodes_degrees(const std::vector<std::vector<opcode>> &opss)
{
  std::vector<int> result;
  for (const auto & ops : opss)
  {
    result.push_back(opcodes_degree(ops));
  }
  return result;
}

std::vector<std::vector<opcode>> parse_opcodess(const std::string &s)
{
  std::istringstream i(s);
  std::vector<std::vector<opcode>> result;
  std::vector<opcode> current;
  while (true)
  {
    std::string word;
    i >> word;
    if (word == "")
    {
      break;
    }
    int j;
    for (j = 0; j < op_count; ++j)
    {
      if (op_string[j] == word)
      {
        opcode op((opcode(j)));
        current.push_back(op);
        if (op == op_add)
        {
          result.push_back(current);
          current = std::vector<opcode>();
        }
        break;
      }
    }
    if (j == op_count)
    {
      throw std::out_of_range{"unrecognized opcode"};
    }
  }
  if (! current.empty())
  {
    result.push_back(current); // invalid formula, checked later
  }
  return result;
}

std::vector<opcode> parse_opcodes(const std::string &s)
{
  std::vector<std::vector<opcode>> r = parse_opcodess(s);
  if (r.size() == 0)
  {
    return std::vector<opcode>();
  }
  if (r.size() == 1)
  {
    return r[0];
  }
  else
  {
    throw std::out_of_range{"too many opcodes"};
  }
}

bool validate_opcodes(const std::vector<opcode> &ops)
{
  if (ops.empty())
  {
    return false;
  }
  bool have_store = false;
  bool have_add = false;
  for (const auto & op : ops)
  {
    if (have_add)
    {
      return false;
    }
    switch (op)
    {
      case op_add: have_add = true; break;
      case op_store: have_store = true; break;
      case op_mul: if (! have_store) return false; break;
      case op_sqr:
      case op_absx:
      case op_absy:
      case op_negx:
      case op_negy:
        break;
    }
  }
  if (! have_add)
  {
    return false;
  }
  if (! (opcodes_degree(ops) >= 2))
  {
    // return false; // can cause problems in GUI editor // FIXME
  }
  return true;
}

bool validate_opcodess(const std::vector<std::vector<opcode>> &opss)
{
  if (opss.empty())
  {
    return false;
  }
  for (const auto & ops : opss)
  {
    if (! validate_opcodes(ops))
    {
      return false;
    }
  }
  return true;
}

const char * const op_string[op_count] = { "add", "store", "mul", "sqr", "absx", "absy", "negx", "negy" };
