// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <chrono>
#include <thread>
#include <vector>

#include <sys/stat.h>

#include <mpreal.h>
#include <toml.hpp>

#include "bla.h"
#include "colour.h"
#include "display.h"
#include "engine.h"
#include "floatexp.h"
#include "hybrid.h"
#include "map.h"
#include "param.h"
#include "softfloat.h"
#include "stats.h"
#include "types.h"

const char *nt_string[
#ifdef HAVE_FLOAT128
  7
#else
  6
#endif
] = { "none", "float", "double", "long double", "floatexp", "softfloat"
#ifdef HAVE_FLOAT128
  , "float128"
#endif
};

number_type nt_current = nt_none;

std::string pref_path = ""; // current working directory; updated by front-end

#ifdef HAVE_FLOAT128
std::vector<std::vector<complex<float128>>> Zq;
#endif
std::vector<std::vector<complex<softfloat>>> Zsf;
std::vector<std::vector<complex<floatexp>>> Zfe;
std::vector<std::vector<complex<long double>>> Zld;
std::vector<std::vector<complex<double>>> Zd;
std::vector<std::vector<complex<float>>> Zf;

void delete_ref()
{
#ifdef HAVE_FLOAT128
  std::vector<std::vector<complex<float128>>>().swap(Zq);
#endif
  std::vector<std::vector<complex<softfloat>>>().swap(Zsf);
  std::vector<std::vector<complex<floatexp>>>().swap(Zfe);
  std::vector<std::vector<complex<long double>>>().swap(Zld);
  std::vector<std::vector<complex<double>>>().swap(Zd);
  std::vector<std::vector<complex<float>>>().swap(Zf);
}

#ifdef HAVE_FLOAT128
std::vector<blasR2<float128>> Bq;
#endif
std::vector<blasR2<softfloat>> Bsf;
std::vector<blasR2<floatexp>> Bfe;
std::vector<blasR2<long double>> Bld;
std::vector<blasR2<double>> Bd;
std::vector<blasR2<float>> Bf;

void delete_bla()
{
#ifdef HAVE_FLOAT128
  std::vector<blasR2<float128>>().swap(Bq);
#endif
  std::vector<blasR2<softfloat>>().swap(Bsf);
  std::vector<blasR2<floatexp>>().swap(Bfe);
  std::vector<blasR2<long double>>().swap(Bld);
  std::vector<blasR2<double>>().swap(Bd);
  std::vector<blasR2<float>>().swap(Bf);
}

count_t getM(number_type nt, count_t phase)
{
  switch (nt)
  {
    case nt_none:
      return 0;
    case nt_float:
      return Zf[phase].size();
    case nt_double:
      return Zd[phase].size();
    case nt_longdouble:
      return Zld[phase].size();
    case nt_floatexp:
      return Zfe[phase].size();
    case nt_softfloat:
      return Zsf[phase].size();
#ifdef HAVE_FLOAT128
    case nt_float128:
      return Zq[phase].size();
#endif
  }
  return 0;
}

bool number_type_available(const param &par, number_type nt)
{
  return par.p.algorithm.number_types.empty() || std::find(par.p.algorithm.number_types.begin(), par.p.algorithm.number_types.end(), std::string(nt_string[nt])) != par.p.algorithm.number_types.end();
}

struct nt_characteristic
{
  number_type type;
  int mantissa_bits;
  int exponent_bits;
  double iterations_per_second;
};

bool comparing_iterations_per_second(const nt_characteristic &a, const nt_characteristic &b) 
{
  return a.iterations_per_second < b.iterations_per_second;
}

template <typename real>
void compute_characteristic_thread(volatile bool *running, count_t *iterations)
{
  complex<real> C(-1.755, 0.001);
  complex<real> Z(0.0, 0.0);
  complex<real> c(1e-6, 1e-6);
  complex<real> z(0.0, 0.0);
  count_t n = 0;
  while (*running)
  {
    Z = sqr(Z) + C;
    z = (2 * Z + z) * z + c;
    n++;
  }
  *running = z.x > z.y; // ensure work is done
  *iterations = n;
}

template<typename real>
nt_characteristic compute_characteristic(number_type type)
{
  using std::isinf, ::isinf;
  using std::ldexp, ::ldexp;
  int mantissa_bits = 0;
  int exponent_bits = 0;
  const real one = 1;
  for (int bits = 0; ; bits++)
  {
    if (one + ldexp(one, -bits) == one)
    {
      mantissa_bits = bits;
      break;
    }
  }
  for (int bits = 0; ; bits++)
  {
    if (isinf(ldexp(one, 1 << bits)))
    {
      exponent_bits = bits + 1; // signed
      break;
    }
  }
  volatile bool running = true;
  count_t iterations_per_second = 0;
  std::thread bg(compute_characteristic_thread<real>, &running, &iterations_per_second);
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  running = false;
  bg.join();
#ifdef __EMSCRIPTEN__
  std::cerr << "found number type " << type << " " << mantissa_bits << "." << exponent_bits << " " << double(iterations_per_second) << std::endl;
#endif
  return { type, mantissa_bits, exponent_bits, double(iterations_per_second) };
}

std::vector<nt_characteristic> nt_characteristics;

void load_characteristics(const std::string &filename)
{
  nt_characteristics.clear();
  try
  {
    struct stat sbuf;
    if (0 == stat(filename.c_str(), &sbuf))
    {
      std::ifstream ifs(filename, std::ios_base::binary);
      ifs.exceptions(std::ifstream::badbit);
      if (ifs.is_open())
      {
        auto t = toml::parse(ifs, filename);
#define LOAD(type) nt_characteristics.push_back(\
        { type \
        , toml::find<int>(t, nt_string[type], "mantissa") \
        , toml::find<int>(t, nt_string[type], "exponent") \
        , toml::find<double>(t, nt_string[type], "speed") \
        });
        LOAD(nt_float)
        LOAD(nt_double)
        LOAD(nt_longdouble)
        LOAD(nt_floatexp)
        LOAD(nt_softfloat)
#ifdef HAVE_FLOAT128
        LOAD(nt_float128)
#endif
#undef LOAD
      }
    }
  }
  catch (...)
  {
    nt_characteristics.clear();
    std::cerr << "ERROR loading number type characteristics" << std::endl;
  }
}

void save_characteristics(const std::string &filename)
{
  try
  {
    std::ofstream ofs(filename, std::ios_base::binary);
    ofs.exceptions(std::ifstream::badbit);
    toml::value data;
    for (const auto &c : nt_characteristics)
    {
      data[nt_string[c.type]] = toml::value
        { { "mantissa", toml::value(c.mantissa_bits) }
        , { "exponent", toml::value(c.exponent_bits) }
        , { "speed", toml::value(c.iterations_per_second) }
        };
    }
    ofs << std::setprecision(17) << data;
  }
  catch (const std::exception &e)
  {
    std::cerr << "ERROR saving number type characteristics" << std::endl;
    std::cerr << e.what() << std::endl;
  }
}

void compute_characteristics()
{
#ifdef __EMSCRIPTEN__
  std::cerr << "computing characteristics of available number types..." << std::endl;
#endif
  nt_characteristics =
    { compute_characteristic<float>(nt_float)
    , compute_characteristic<double>(nt_double)
    , compute_characteristic<long double>(nt_longdouble)
    , compute_characteristic<floatexp>(nt_floatexp)
    , compute_characteristic<softfloat>(nt_softfloat)
#ifdef HAVE_FLOAT128
    , compute_characteristic<float128>(nt_float128)
#endif
    };
}

void populate_number_type_wisdom(void)
{
  const std::string filename = pref_path + "number-type-wisdom.toml";
  if (nt_characteristics.empty())
  {
    load_characteristics(filename);
  }
  if (nt_characteristics.empty())
  {
    compute_characteristics();
    save_characteristics(filename);
#ifdef __EMSCRIPTEN__
    syncfs();
#endif
  }
}

number_type choose_number_type(const param &par, count_t pixel_spacing_exponent, count_t pixel_spacing_precision)
{
  std::vector<nt_characteristic> candidates;
  for (auto c : nt_characteristics)
  {
    if (pixel_spacing_exponent < (count_t(1) << c.exponent_bits) >> 1 &&
        pixel_spacing_precision < c.mantissa_bits &&
        number_type_available(par, c.type))
    {
      candidates.push_back(c);
    }
  }
  if (candidates.empty())
  {
    return nt_none;
  }
  number_type nt = std::max_element(candidates.begin(), candidates.end(), comparing_iterations_per_second)->type;
  return nt;
} 

void reference_thread(stats &sta, param &par, progress_t *progress, bool *running, bool *ended)
{
  floatexp pixel_spacing = 1 / (par.zoom * par.p.image.height);
  complex<mpreal> offset;
  offset.x.set_prec(par.center.x.get_prec());
  offset.y.set_prec(par.center.y.get_prec());
  offset = par.center - par.reference;
  floatexp pixel_precision = std::max
    ( std::max(abs(offset.x / pixel_spacing), abs(offset.y / pixel_spacing))
    , hypot(floatexp(par.p.image.width), floatexp(par.p.image.height))
    );
  number_type nt = choose_number_type(par, std::max(24, 24 - pixel_spacing.exp), pixel_precision.exp);
  bool have_reference = false;
  bool have_bla = false;
  if (par.p.algorithm.reuse_reference && nt_current != nt_none && nt != nt_none)
  {
    have_reference = nt == nt_current;
  }
  if (par.p.algorithm.reuse_bilinear_approximation && nt != nt_none)
  {
    have_bla = nt == nt_current;
  }
  if (nt != nt_current || nt == nt_none || nt_current == nt_none)
  {
    // will be using a new reference in image center
    nt = choose_number_type(par, std::max(24, 24 - pixel_spacing.exp), hypot(floatexp(par.p.image.width), floatexp(par.p.image.height)).exp);
    have_reference = false;
    have_bla = false;
  }
  count_t maximum_reference_iterations = par.p.bailout.maximum_reference_iterations;
  if (par.p.algorithm.lock_maximum_reference_iterations_to_period)
  {
    maximum_reference_iterations = par.p.reference.period;
  }
  const count_t count = par.p.formula.per.size();
  if (have_reference)
  {
    for (count_t i = 0; i < count; ++i)
    {
      progress[i] = getM(nt, i) / progress_t(maximum_reference_iterations);
    }
  }
  else
  {
    par.reference.x.set_prec(par.center.x.get_prec());
    par.reference.y.set_prec(par.center.y.get_prec());
    par.reference = par.center;
    offset = par.center - par.reference;
    pixel_precision = std::max
      ( std::max(abs(offset.x / pixel_spacing), abs(offset.y / pixel_spacing))
      , hypot(floatexp(par.p.image.width), floatexp((par.p.image.height)))
      );
    delete_ref();
    switch (nt)
    {
#ifdef HAVE_FLOAT128
      case nt_float128:
        Zq.resize(count);
        for (count_t phase = 0; phase < count; ++phase)
        {
          Zq[phase].resize(maximum_reference_iterations);
        }
        hybrid_references(Zq, par.p.formula, maximum_reference_iterations, par.reference, &progress[0], running);
        break;
#endif
      case nt_softfloat:
        Zsf.resize(count);
        for (count_t phase = 0; phase < count; ++phase)
        {
          Zsf[phase].resize(maximum_reference_iterations);
        }
        hybrid_references(Zsf, par.p.formula, maximum_reference_iterations, par.reference, &progress[0], running);
        break;
      case nt_floatexp:
        Zfe.resize(count);
        for (count_t phase = 0; phase < count; ++phase)
        {
          Zfe[phase].resize(maximum_reference_iterations);
        }
        hybrid_references(Zfe, par.p.formula, maximum_reference_iterations, par.reference, &progress[0], running);
        break;
      case nt_longdouble:
        Zld.resize(count);
        for (count_t phase = 0; phase < count; ++phase)
        {
          Zld[phase].resize(maximum_reference_iterations);
        }
        hybrid_references(Zld, par.p.formula, maximum_reference_iterations, par.reference, &progress[0], running);
        break;
      case nt_double:
        Zd.resize(count);
        for (count_t phase = 0; phase < count; ++phase)
        {
          Zd[phase].resize(maximum_reference_iterations);
        }
        hybrid_references(Zd, par.p.formula, maximum_reference_iterations, par.reference, &progress[0], running);
        break;
      case nt_float:
        Zf.resize(count);
        for (count_t phase = 0; phase < count; ++phase)
        {
          Zf[phase].resize(maximum_reference_iterations);
        }
        hybrid_references(Zf, par.p.formula, maximum_reference_iterations, par.reference, &progress[0], running);
        break;
      case nt_none:
        break;
    }
    have_bla = false;
  }
  if (have_bla)
  {
    for (count_t i = 0; i < count; ++i)
    {
      progress[count + i] = 1;
    }
  }
  else
  {
    const count_t height = par.p.image.height;
    const floatexp pixel_spacing = 4 / par.zoom / height;
    const count_t bits = 24; // FIXME
    const float precision = count_t(1) << bits;
    using std::hypot, ::hypot;
    delete_bla();
    switch (nt)
    {
      case nt_none: break;
      case nt_float: hybrid_blas(Bf, Zf, par.p.formula, float(pixel_precision), float(pixel_spacing), float(precision), &progress[count], running); break;
      case nt_double: hybrid_blas(Bd, Zd, par.p.formula, double(pixel_precision), double(pixel_spacing), double(precision), &progress[count], running); break;
      case nt_longdouble: hybrid_blas(Bld, Zld, par.p.formula, (long double)(pixel_precision), (long double)(pixel_spacing), (long double)(precision), &progress[count], running); break;
      case nt_floatexp: hybrid_blas(Bfe, Zfe, par.p.formula, floatexp(pixel_precision), floatexp(pixel_spacing), floatexp(precision), &progress[count], running); break;
      case nt_softfloat: hybrid_blas(Bsf, Zsf, par.p.formula, softfloat(pixel_precision), softfloat(pixel_spacing), softfloat(precision), &progress[count], running); break;
#ifdef HAVE_FLOAT128
      case nt_float128: hybrid_blas(Bq, Zq, par.p.formula, float128(pixel_precision), float128(pixel_spacing), float128(precision), &progress[count], running); break;
#endif
    }
  }
  nt_current = nt;
  reset(sta);
  *ended = true;
}

void subframe_thread(map &out, stats &sta, const param &par, const count_t subframe, progress_t *progress, bool *running, bool *ended)
{
  complex<mpreal> offset;
  offset.x.set_prec(par.center.x.get_prec());
  offset.y.set_prec(par.center.y.get_prec());
  offset = par.center - par.reference;
  switch (nt_current)
  {
    case nt_none: break;
    case nt_float: hybrid_render(out, sta, par.p.formula, Bf, subframe, par, float(par.zoom), complex<float>(float(offset.x), float(offset.y)), Zf, progress, running); break;
    case nt_double: hybrid_render(out, sta, par.p.formula, Bd, subframe, par, double(par.zoom), complex<double>(double(offset.x), double(offset.y)), Zd, progress, running); break;
    case nt_longdouble: hybrid_render(out, sta, par.p.formula, Bld, subframe, par, (long double)(par.zoom), complex<long double>((long double)(offset.x), (long double)(offset.y)), Zld, progress, running); break;
    case nt_floatexp: hybrid_render(out, sta, par.p.formula, Bfe, subframe, par, floatexp(par.zoom), complex<floatexp>(floatexp(offset.x), floatexp(offset.y)), Zfe, progress, running); break;
    case nt_softfloat: hybrid_render(out, sta, par.p.formula, Bsf, subframe, par, softfloat(par.zoom), complex<softfloat>(softfloat(offset.x), softfloat(offset.y)), Zsf, progress, running); break;
#ifdef HAVE_FLOAT128
    case nt_float128: hybrid_render(out, sta, par.p.formula, Bq, subframe, par, float128(par.zoom), complex<float128>(float128(offset.x), float128(offset.y)), Zq, progress, running); break;
#endif
  }
  *ended = true;
}
