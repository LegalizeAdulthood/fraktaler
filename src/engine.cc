// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <chrono>
#include <thread>
#include <vector>

#include <mpreal.h>
#include <toml.hpp>

#include "colour.h"
#include "display.h"
#include "engine.h"
#include "floatexp.h"
#include "formula.h"
#include "map.h"
#include "param.h"
#include "softfloat.h"
#include "stats.h"
#include "types.h"

const char *nt_string[7] = { "none", "float", "double", "long double", "floatexp", "softfloat", "float128" };

number_type nt_current = nt_none;

std::string pref_path = ""; // current working directory; updated by front-end

std::vector<complex<float128>> Zq;
std::vector<complex<softfloat>> Zsf;
std::vector<complex<floatexp>> Zfe;
std::vector<complex<long double>> Zld;
std::vector<complex<double>> Zd;
std::vector<complex<float>> Zf;

blasC<float128> *BCq = nullptr;
blasC<softfloat> *BCsf = nullptr;
blasC<floatexp> *BCfe = nullptr;
blasC<long double> *BCld = nullptr;
blasC<double> *BCd = nullptr;
blasC<float> *BCf = nullptr;

blasR2<float128> *BR2q = nullptr;
blasR2<softfloat> *BR2sf = nullptr;
blasR2<floatexp> *BR2fe = nullptr;
blasR2<long double> *BR2ld = nullptr;
blasR2<double> *BR2d = nullptr;
blasR2<float> *BR2f = nullptr;

void delete_bla()
{
  delete BCq; BCq = nullptr;
  delete BCsf; BCsf = nullptr;
  delete BCfe; BCfe = nullptr;
  delete BCld; BCld = nullptr;
  delete BCd; BCd = nullptr;
  delete BCf; BCf = nullptr;
  delete BR2q; BR2q = nullptr;
  delete BR2sf; BR2sf = nullptr;
  delete BR2fe; BR2fe = nullptr;
  delete BR2ld; BR2ld = nullptr;
  delete BR2d; BR2d = nullptr;
  delete BR2f; BR2f = nullptr;
}

count_t getM(number_type nt)
{
  switch (nt)
  {
    case nt_none:
      return 0;
    case nt_float:
      return Zf.size();
    case nt_double:
      return Zd.size();
    case nt_longdouble:
      return Zld.size();
    case nt_floatexp:
      return Zfe.size();
    case nt_softfloat:
      return Zsf.size();
    case nt_float128:
      return Zq.size();
  }
  return 0;
}

bool convert_reference(const number_type to, const number_type from)
{
  count_t M;
  bool converted = true;
  switch (to)
  {
    case nt_none:
      switch (from)
      {
        case nt_none: break;
        case nt_float:
          Zf.clear();
          M = 0;
          break;
        case nt_double:
          Zd.clear();
          M = 0;
          break;
        case nt_longdouble:
          Zld.clear();
          M = 0;
          break;
        case nt_floatexp:
          Zfe.clear();
          M = 0;
          break;
        case nt_softfloat:
          Zsf.clear();
          M = 0;
          break;
        case nt_float128:
          Zq.clear();
          M = 0;
          break;
      }
      break;

    case nt_float:
      switch (from)
      {
        case nt_none: converted = false; break;
        case nt_float: break;
        case nt_double:
          Zf.resize(Zd.size());
          M = Zf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<double> Z = Zd[m];
            Zf[m] = complex<float>(float(Z.x), float(Z.y));
          }
          Zd.clear();
          break;
        case nt_longdouble:
          Zf.resize(Zld.size());
          M = Zf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<long double> Z = Zld[m];
            Zf[m] = complex<float>(float(Z.x), float(Z.y));
          }
          Zld.clear();
          break;
        case nt_floatexp:
          Zf.resize(Zfe.size());
          M = Zf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<floatexp> Z = Zfe[m];
            Zf[m] = complex<float>(float(Z.x), float(Z.y));
          }
          Zfe.clear();
          break;
        case nt_softfloat:
          Zf.resize(Zsf.size());
          M = Zf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<softfloat> Z = Zsf[m];
            Zf[m] = complex<float>(float(Z.x), float(Z.y));
          }
          Zsf.clear();
          break;
        case nt_float128:
          Zf.resize(Zq.size());
          M = Zf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<float128> Z = Zq[m];
            Zf[m] = complex<float>(float(Z.x), float(Z.y));
          }
          Zq.clear();
          break;
      }
      break;

    case nt_double:
      switch (from)
      {
        case nt_none: converted = false; break;
        case nt_double: break;
        case nt_float:
          Zd.resize(Zf.size());
          M = Zd.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<float> Z = Zf[m];
            Zd[m] = complex<double>(double(Z.x), double(Z.y));
          }
          Zf.clear();
          break;
        case nt_longdouble:
          Zd.resize(Zld.size());
          M = Zd.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<long double> Z = Zld[m];
            Zd[m] = complex<double>(double(Z.x), double(Z.y));
          }
          Zld.clear();
          break;
        case nt_floatexp:
          Zd.resize(Zfe.size());
          M = Zd.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<floatexp> Z = Zfe[m];
            Zd[m] = complex<double>(double(Z.x), double(Z.y));
          }
          Zfe.clear();
          break;
        case nt_softfloat:
          Zd.resize(Zsf.size());
          M = Zd.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<softfloat> Z = Zsf[m];
            Zd[m] = complex<double>(double(Z.x), double(Z.y));
          }
          Zsf.clear();
          break;
        case nt_float128:
          Zd.resize(Zq.size());
          M = Zd.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<float128> Z = Zq[m];
            Zd[m] = complex<double>(double(Z.x), double(Z.y));
          }
          Zq.clear();
          break;
      }
      break;

    case nt_longdouble:
      switch (from)
      {
        case nt_none: converted = false; break;
        case nt_longdouble: break;
        case nt_float:
          Zld.resize(Zf.size());
          M = Zd.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<float> Z = Zf[m];
            Zld[m] = complex<long double>((long double)(Z.x), (long double)(Z.y));
          }
          Zf.clear();
          break;
        case nt_double:
          Zld.resize(Zd.size());
          M = Zd.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<double> Z = Zd[m];
            Zld[m] = complex<long double>((long double)(Z.x), (long double)(Z.y));
          }
          Zd.clear();
          break;
        case nt_floatexp:
          Zld.resize(Zfe.size());
          M = Zld.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<floatexp> Z = Zfe[m];
            Zld[m] = complex<long double>((long double)(Z.x), (long double)(Z.y));
          }
          Zfe.clear();
          break;
        case nt_softfloat:
          Zld.resize(Zsf.size());
          M = Zld.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<softfloat> Z = Zsf[m];
            Zld[m] = complex<long double>((long double)(Z.x), (long double)(Z.y));
          }
          Zsf.clear();
          break;
        case nt_float128:
          Zld.resize(Zq.size());
          M = Zld.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<float128> Z = Zq[m];
            Zld[m] = complex<long double>((long double)(Z.x), (long double)(Z.y));
          }
          Zq.clear();
          break;
      }
      break;

    case nt_floatexp:
      switch (from)
      {
        case nt_none: converted = false; break;
        case nt_floatexp: break;
        case nt_float:
          Zfe.resize(Zf.size());
          M = Zf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<float> Z = Zf[m];
            Zfe[m] = complex<floatexp>(floatexp(Z.x), floatexp(Z.y));
          }
          Zf.clear();
          break;
        case nt_double:
          Zfe.resize(Zd.size());
          M = Zd.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<double> Z = Zd[m];
            Zfe[m] = complex<floatexp>(floatexp(Z.x), floatexp(Z.y));
          }
          Zd.clear();
          break;
        case nt_longdouble:
          Zfe.resize(Zld.size());
          M = Zfe.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<long double> Z = Zld[m];
            Zfe[m] = complex<floatexp>(floatexp(Z.x), floatexp(Z.y));
          }
          Zld.clear();
          break;
        case nt_softfloat:
          Zfe.resize(Zsf.size());
          M = Zfe.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<softfloat> Z = Zsf[m];
            Zfe[m] = complex<floatexp>(floatexp(Z.x), floatexp(Z.y));
          }
          Zsf.clear();
          break;
        case nt_float128:
          Zfe.resize(Zq.size());
          M = Zfe.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<float128> Z = Zq[m];
            Zfe[m] = complex<floatexp>(floatexp(Z.x), floatexp(Z.y));
          }
          Zq.clear();
          break;
      }
      break;

    case nt_softfloat:
      switch (from)
      {
        case nt_none: converted = false; break;
        case nt_softfloat: break;
        case nt_float:
          Zsf.resize(Zf.size());
          M = Zf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<float> Z = Zf[m];
            Zsf[m] = complex<softfloat>(softfloat(Z.x), softfloat(Z.y));
          }
          Zf.clear();
          break;
        case nt_double:
          Zsf.resize(Zd.size());
          M = Zd.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<double> Z = Zd[m];
            Zsf[m] = complex<softfloat>(softfloat(Z.x), softfloat(Z.y));
          }
          Zd.clear();
          break;
        case nt_longdouble:
          Zsf.resize(Zld.size());
          M = Zsf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<long double> Z = Zld[m];
            Zsf[m] = complex<softfloat>(softfloat(Z.x), softfloat(Z.y));
          }
          Zld.clear();
          break;
        case nt_floatexp:
          Zsf.resize(Zfe.size());
          M = Zsf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<floatexp> Z = Zfe[m];
            Zsf[m] = complex<softfloat>(softfloat(Z.x), softfloat(Z.y));
          }
          Zfe.clear();
          break;
        case nt_float128:
          Zsf.resize(Zq.size());
          M = Zsf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<float128> Z = Zq[m];
            Zsf[m] = complex<softfloat>(softfloat(Z.x), softfloat(Z.y));
          }
          Zq.clear();
          break;
      }
      break;

    case nt_float128:
      switch (from)
      {
        case nt_none: converted = false; break;
        case nt_float128: break;
        case nt_float:
          Zq.resize(Zf.size());
          M = Zf.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<float> Z = Zf[m];
            Zq[m] = complex<float128>(float128(Z.x), float128(Z.y));
          }
          Zf.clear();
          break;
        case nt_double:
          Zq.resize(Zd.size());
          M = Zd.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<double> Z = Zd[m];
            Zq[m] = complex<float128>(float128(Z.x), float128(Z.y));
          }
          Zd.clear();
          break;
        case nt_longdouble:
          Zq.resize(Zld.size());
          M = Zq.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<long double> Z = Zld[m];
            Zq[m] = complex<float128>(float128(Z.x), float128(Z.y));
          }
          Zld.clear();
          break;
        case nt_floatexp:
          Zq.resize(Zfe.size());
          M = Zq.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<floatexp> Z = Zfe[m];
            Zq[m] = complex<float128>(float128(Z.x), float128(Z.y));
          }
          Zfe.clear();
          break;
        case nt_softfloat:
          Zq.resize(Zsf.size());
          M = Zq.size();
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<softfloat> Z = Zsf[m];
            Zq[m] = complex<float128>(float128(Z.x), float128(Z.y));
          }
          Zsf.clear();
          break;
      }
      break;
  }
  return converted;
}

bool convert_bla(const number_type to, const number_type from)
{
  (void) to;
  (void) from;
  return false;
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
  bool running = true;
  count_t iterations_per_second = 0;
  std::thread bg(compute_characteristic_thread<real>, &running, &iterations_per_second);
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  running = false;
  bg.join();
  return { type, mantissa_bits, exponent_bits, double(iterations_per_second) };
}

std::vector<nt_characteristic> nt_characteristics;

void load_characteristics(const std::string &filename)
{
  nt_characteristics.clear();
  try
  {
    std::ifstream ifs(filename, std::ios_base::binary);
    ifs.exceptions(std::ifstream::badbit);
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
    LOAD(nt_float128)
#undef LOAD
  }
  catch (std::exception &e)
  {
    nt_characteristics.clear();
    std::cerr << "ERROR loading number type characteristics" << std::endl;
    std::cerr << e.what() << std::endl;
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
  catch (std::exception &e)
  {
    std::cerr << "ERROR saving number type characteristics" << std::endl;
    std::cerr << e.what() << std::endl;
  }
}

void compute_characteristics()
{
  nt_characteristics =
    { compute_characteristic<float>(nt_float)
    , compute_characteristic<double>(nt_double)
    , compute_characteristic<long double>(nt_longdouble)
    , compute_characteristic<floatexp>(nt_floatexp)
    , compute_characteristic<softfloat>(nt_softfloat)
    , compute_characteristic<float128>(nt_float128)
    };
}

number_type choose_number_type(const param &par, count_t pixel_spacing_exponent, count_t pixel_spacing_precision)
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
  }
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
//std::cerr << "choosing none for exponent " << pixel_spacing_exponent << " and mantissa " << pixel_spacing_precision << std::endl;
    return nt_none;
  }
  number_type nt = std::max_element(candidates.begin(), candidates.end(), comparing_iterations_per_second)->type;
//std::cerr << "choosing " << nt_string[nt] << " for exponent " << pixel_spacing_exponent << " and mantissa " << pixel_spacing_precision << std::endl;
  return nt;
} 

void reference_thread(stats &sta, const formula *form, param &par, progress_t *progress, bool *running, bool *ended)
{
  floatexp pixel_spacing = 1 / (par.zoom * par.p.image.height);
  complex<mpreal> offset;
  offset.x.set_prec(par.center.x.get_prec());
  offset.y.set_prec(par.center.y.get_prec());
  offset = par.center - par.reference;
  floatexp pixel_precision = std::max
    ( std::max(abs(offset.x / pixel_spacing), abs(offset.y / pixel_spacing))
    , hypot(floatexp(par.p.image.width), floatexp((par.p.image.height)))
    );
  number_type nt = choose_number_type(par, std::max(24, 24 - pixel_spacing.exp), pixel_precision.exp);
  bool have_reference = false;
  bool have_bla = false;
  if (par.p.algorithm.reuse_reference && nt_current != nt_none && nt != nt_none)
  {
    have_reference = convert_reference(nt, nt_current);
  }
  if (par.p.algorithm.reuse_bilinear_approximation && nt != nt_none)
  {
    have_bla = convert_bla(nt, nt_current);
  }
  if (nt != nt_current || nt == nt_none)
  {
    // will be using a new reference in image center
    nt = choose_number_type(par, std::max(24, 24 - pixel_spacing.exp), hypot(floatexp(par.p.image.width), floatexp((par.p.image.height))).exp);
    have_reference = false;
    have_bla = false;
  }
  count_t maximum_reference_iterations = par.p.bailout.maximum_reference_iterations;
  if (par.p.algorithm.lock_maximum_reference_iterations_to_period)
  {
    maximum_reference_iterations = par.p.reference.period;
  }
  if (have_reference)
  {
    progress[0] = getM(nt) / progress_t(maximum_reference_iterations);
    progress[1] = 1;
  }
  else
  {
    //std::cerr << "computing reference" << std::endl;
    par.reference.x.set_prec(par.center.x.get_prec());
    par.reference.y.set_prec(par.center.y.get_prec());
    par.reference = par.center;
    offset = par.center - par.reference;
    pixel_precision = std::max
      ( std::max(abs(offset.x / pixel_spacing), abs(offset.y / pixel_spacing))
      , hypot(floatexp(par.p.image.width), floatexp((par.p.image.height)))
      );
    count_t M;
    switch (nt)
    {
      case nt_float128:
        Zq.resize(maximum_reference_iterations);
        Zsf.clear();
        Zfe.clear();
        Zld.clear();
        Zd.clear();
        Zf.clear();
        M = form->reference(&Zq[0], maximum_reference_iterations, par.reference, &progress[0], running);
        Zq.resize(M);
        break;
      case nt_softfloat:
        Zq.clear();
        Zsf.resize(maximum_reference_iterations);
        Zfe.clear();
        Zld.clear();
        Zd.clear();
        Zf.clear();
        M = form->reference(&Zsf[0], maximum_reference_iterations, par.reference, &progress[0], running);
        Zsf.resize(M);
        break;
      case nt_floatexp:
        Zq.clear();
        Zsf.clear();
        Zfe.resize(maximum_reference_iterations);
        Zld.clear();
        Zd.clear();
        Zf.clear();
        M = form->reference(&Zfe[0], maximum_reference_iterations, par.reference, &progress[0], running);
        Zfe.resize(M);
        break;
      case nt_longdouble:
        Zq.clear();
        Zsf.clear();
        Zfe.clear();
        Zld.resize(maximum_reference_iterations);
        Zd.clear();
        Zf.clear();
        M = form->reference(&Zld[0], maximum_reference_iterations, par.reference, &progress[0], running);
        Zld.resize(M);
        break;
      case nt_double:
        Zq.clear();
        Zsf.clear();
        Zfe.clear();
        Zld.clear();
        Zd.resize(maximum_reference_iterations);
        Zf.clear();
        M = form->reference(&Zd[0], maximum_reference_iterations, par.reference, &progress[0], running);
        Zd.resize(M);
        break;
      case nt_float:
        Zq.clear();
        Zsf.clear();
        Zfe.clear();
        Zld.clear();
        Zd.clear();
        Zf.resize(maximum_reference_iterations);
        M = form->reference(&Zf[0], maximum_reference_iterations, par.reference, &progress[0], running);
        Zf.resize(M);
        break;
      case nt_none:
        Zq.clear();
        Zsf.clear();
        Zfe.clear();
        Zld.clear();
        Zd.clear();
        Zf.clear();
        M = 0;
        break;
    }
    progress[1] = 1;
  }
  if (have_bla)
  {
    progress[2] = 1;
  }
  else
  {
    //std::cerr << "computing bla" << std::endl;
    const count_t height = par.p.image.height;
    const floatexp pixel_spacing = 4 / par.zoom / height;
    const count_t bits = 24; // FIXME
    const float precision = count_t(1) << bits;
    using std::hypot, ::hypot;
    delete_bla();
    if (form->complex_analytic())
    {
      const formulaCbase *fc = static_cast<const formulaCbase *>(form);
      switch (nt)
      {
        case nt_none: break;
        case nt_float: BCf = fc->bla(&Zf[0], Zf.size(), float(pixel_precision), float(pixel_spacing), float(precision), &progress[2], running); break;
        case nt_double: BCd = fc->bla(&Zd[0], Zd.size(),double(pixel_precision), double(pixel_spacing), double(precision), &progress[2], running); break;
        case nt_longdouble: BCld = fc->bla(&Zld[0], Zld.size(), (long double)(pixel_precision), (long double)(pixel_spacing), (long double)(precision), &progress[2], running); break;
        case nt_floatexp: BCfe = fc->bla(&Zfe[0], Zfe.size(), floatexp(pixel_precision), floatexp(pixel_spacing), floatexp(precision), &progress[2], running); break;
        case nt_softfloat: BCsf = fc->bla(&Zsf[0], Zsf.size(), softfloat(pixel_precision), softfloat(pixel_spacing), softfloat(precision), &progress[2], running); break;
        case nt_float128: BCq = fc->bla(&Zq[0], Zq.size(), float128(pixel_precision), float128(pixel_spacing), float128(precision), &progress[2], running); break;
      }
    }
    else
    {
      const formulaR2base *fr2 = static_cast<const formulaR2base *>(form);
      switch (nt)
      {
        case nt_none: break;
        case nt_float: BR2f = fr2->bla(&Zf[0], Zf.size(), float(pixel_precision), float(pixel_spacing), float(precision), &progress[2], running); break;
        case nt_double: BR2d = fr2->bla(&Zd[0], Zd.size(), double(pixel_precision), double(pixel_spacing), double(precision), &progress[2], running); break;
        case nt_longdouble: BR2ld = fr2->bla(&Zld[0], Zld.size(), (long double)(pixel_precision), (long double)(pixel_spacing), (long double)(precision), &progress[2], running); break;
        case nt_floatexp: BR2fe = fr2->bla(&Zfe[0], Zfe.size(), floatexp(pixel_precision), floatexp(pixel_spacing), floatexp(precision), &progress[2], running); break;
        case nt_softfloat: BR2sf = fr2->bla(&Zsf[0], Zsf.size(), softfloat(pixel_precision), softfloat(pixel_spacing), softfloat(precision), &progress[2], running); break;
        case nt_float128: BR2q = fr2->bla(&Zq[0], Zq.size(), float128(pixel_precision), float128(pixel_spacing), float128(precision), &progress[2], running); break;
      }
    }
  }
  nt_current = nt;
  reset(sta);
  *ended = true;
}

void subframe_thread(map &out, stats &sta, const formula *form, const param &par, const count_t subframe, progress_t *progress, bool *running, bool *ended)
{
  complex<mpreal> offset;
  offset.x.set_prec(par.center.x.get_prec());
  offset.y.set_prec(par.center.y.get_prec());
  offset = par.center - par.reference;
  if (form->complex_analytic())
  {
    const formulaCbase *fc = static_cast<const formulaCbase *>(form);
    switch (nt_current)
    {
      case nt_none: break;
      case nt_float:
        fc->render(out, sta, BCf, subframe, par, float(par.zoom), complex<float>(float(offset.x), float(offset.y)), Zf.size(), &Zf[0], progress, running);
        break;
      case nt_double:
        fc->render(out, sta, BCd, subframe, par, double(par.zoom), complex<double>(double(offset.x), double(offset.y)), Zd.size(), &Zd[0], progress, running);
        break;
      case nt_longdouble:
        fc->render(out, sta, BCld, subframe, par, (long double)(par.zoom), complex<long double>((long double)(offset.x), (long double)(offset.y)), Zld.size(), &Zld[0], progress, running);
        break;
      case nt_floatexp:
        fc->render(out, sta, BCfe, subframe, par, par.zoom, complex<floatexp>(floatexp(offset.x), floatexp(offset.y)), Zfe.size(), &Zfe[0], progress, running);
        break;
      case nt_softfloat:
        fc->render(out, sta, BCsf, subframe, par, softfloat(par.zoom), complex<softfloat>(softfloat(offset.x), softfloat(offset.y)), Zsf.size(), &Zsf[0], progress, running);
        break;
      case nt_float128:
        fc->render(out, sta, BCq, subframe, par, float128(par.zoom), complex<float128>(convert<float128>(offset.x), convert<float128>(offset.y)), Zq.size(), &Zq[0], progress, running);
        break;
    }
  }
  else
  {
    const formulaR2base *fr2 = static_cast<const formulaR2base *>(form);
    switch (nt_current)
    {
      case nt_none: break;
      case nt_float:
        fr2->render(out, sta, BR2f, subframe, par, float(par.zoom), complex<float>(float(offset.x), float(offset.y)), Zf.size(), &Zf[0], progress, running);
        break;
      case nt_double:
        fr2->render(out, sta, BR2d, subframe, par, double(par.zoom), complex<double>(double(offset.x), double(offset.y)), Zd.size(), &Zd[0], progress, running);
        break;
      case nt_longdouble:
        fr2->render(out, sta, BR2ld, subframe, par, (long double)(par.zoom), complex<long double>((long double)(offset.x), (long double)(offset.y)), Zld.size(), &Zld[0], progress, running);
        break;
      case nt_floatexp:
        fr2->render(out, sta, BR2fe, subframe, par, par.zoom, complex<floatexp>(floatexp(offset.x), floatexp(offset.y)), Zfe.size(), &Zfe[0], progress, running);
        break;
      case nt_softfloat:
        fr2->render(out, sta, BR2sf, subframe, par, softfloat(par.zoom), complex<softfloat>(softfloat(offset.x), softfloat(offset.y)), Zsf.size(), &Zsf[0], progress, running);
        break;
      case nt_float128:
        fr2->render(out, sta, BR2q, subframe, par, float128(par.zoom), complex<float128>(convert<float128>(offset.x), convert<float128>(offset.y)), Zq.size(), &Zq[0], progress, running);
        break;
    }
  }
  *ended = true;
}
