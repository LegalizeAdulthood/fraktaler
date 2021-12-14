// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <vector>

#include <mpreal.h>

#include "colour.h"
#include "display.h"
#include "floatexp.h"
#include "formula.h"
#include "main.h"
#include "map.h"
#include "param.h"
#include "stats.h"

enum number_type
{
  nt_none = 0,
  nt_float = 1,
  nt_double = 2,
  nt_longdouble = 3,
  nt_floatexp = 4
};

const char *nt_string[] = { "none", "float", "double", "long double", "floatexp" };

number_type nt_current = nt_none;

std::vector<complex<floatexp>> Zfe;
std::vector<complex<long double>> Zld;
std::vector<complex<double>> Zd;
std::vector<complex<float>> Zf;

blasC<floatexp> *BCfe = nullptr;
blasC<long double> *BCld = nullptr;
blasC<double> *BCd = nullptr;
blasC<float> *BCf = nullptr;

blasR2<floatexp> *BR2fe = nullptr;
blasR2<long double> *BR2ld = nullptr;
blasR2<double> *BR2d = nullptr;
blasR2<float> *BR2f = nullptr;

void delete_bla()
{
  delete BCfe; BCfe = nullptr;
  delete BCld; BCld = nullptr;
  delete BCd; BCd = nullptr;
  delete BCf; BCf = nullptr;
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
  }
  return 0;
}

bool convert_reference(const number_type to, const number_type from)
{
  count_t M;
  bool converted = true;
  switch (to)
  {
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
      }
      break;
  }
  return converted;
}

bool convert_bla(const number_type to, const number_type from)
{
  return false;
}

void reference_thread(stats &sta, const formula *form, const param &par, progress_t *progress, bool *running, bool *ended)
{
  reset(sta);
  floatexp Zoom = par.Zoom;
  number_type nt = nt_none;
  if (Zoom > e10(1, 4900))
  {
    nt = nt_floatexp;
  }
  else if (Zoom > e10(1, 300))
  {
    nt = nt_longdouble; 
  }
  else if (Zoom > e10(1, 30))
  {
    nt = nt_double;
  }
  else
  {
    nt = nt_float;
  }
  bool have_reference = false;
  bool have_bla = false;
  if (par.ReuseBLA)
  {
    have_bla = convert_bla(nt, nt_current);
  }
  if (par.ReuseReference)
  {
    have_reference = convert_reference(nt, nt_current);
  }
  if (have_reference)
  {
    progress[0] = getM(nt) / progress_t(par.MaxRefIters);
    progress[1] = 1;
  }
  else
  {
    count_t M;
    switch (nt)
    {
      case nt_floatexp:
        Zfe.resize(par.MaxRefIters);
        Zld.clear();
        Zd.clear();
        Zf.clear();
        M = form->reference(&Zfe[0], par.MaxRefIters, par.C, &progress[0], running);
        Zfe.resize(M);
        break;
      case nt_longdouble:
        Zfe.clear();
        Zld.resize(par.MaxRefIters);
        Zd.clear();
        Zf.clear();
        M = form->reference(&Zld[0], par.MaxRefIters, par.C, &progress[0], running);
        Zld.resize(M);
        break;
      case nt_double:
        Zfe.clear();
        Zld.clear();
        Zd.resize(par.MaxRefIters);
        Zf.clear();
        M = form->reference(&Zd[0], par.MaxRefIters, par.C, &progress[0], running);
        Zd.resize(M);
        break;
      case nt_float:
        Zfe.clear();
        Zld.clear();
        Zd.clear();
        Zf.resize(par.MaxRefIters);
        M = form->reference(&Zf[0], par.MaxRefIters, par.C, &progress[0], running);
        Zf.resize(M);
    }
    progress[1] = 1;
  }
  if (have_bla)
  {
    progress[2] = 1;
  }
  else
  {
    const count_t width = par.Width;
    const count_t height = par.Height;
    const floatexp pixel_spacing = 4 / par.Zoom / height;
    const count_t bits = 24; // FIXME
    const float precision = count_t(1) << bits;
    using std::hypot;
    delete_bla();
    if (form->complex_analytic())
    {
      const formulaCbase *fc = static_cast<const formulaCbase *>(form);
      switch (nt)
      {
        case nt_float: BCf = fc->bla(&Zf[0], Zf.size(), hypot(float(width), float(height)), float(pixel_spacing), float(precision), &progress[2], running); break;
        case nt_double: BCd = fc->bla(&Zd[0], Zd.size(), hypot(double(width), double(height)), double(pixel_spacing), double(precision), &progress[2], running); break;
        case nt_longdouble: BCld = fc->bla(&Zld[0], Zld.size(), hypot((long double)(width), (long double)(height)), (long double)(pixel_spacing), (long double)(precision), &progress[2], running); break;
        case nt_floatexp: BCfe = fc->bla(&Zfe[0], Zfe.size(), hypot(floatexp(width), floatexp(height)), floatexp(pixel_spacing), floatexp(precision), &progress[2], running); break;
      }
    }
    else
    {
      const formulaR2base *fr2 = static_cast<const formulaR2base *>(form);
      switch (nt)
      {
        case nt_float: BR2f = fr2->bla(&Zf[0], Zf.size(), hypot(float(width), float(height)), float(pixel_spacing), float(precision), &progress[2], running); break;
        case nt_double: BR2d = fr2->bla(&Zd[0], Zd.size(), hypot(double(width), double(height)), double(pixel_spacing), double(precision), &progress[2], running); break;
        case nt_longdouble: BR2ld = fr2->bla(&Zld[0], Zld.size(), hypot((long double)(width), (long double)(height)), (long double)(pixel_spacing), (long double)(precision), &progress[2], running); break;
        case nt_floatexp: BR2fe = fr2->bla(&Zfe[0], Zfe.size(), hypot(floatexp(width), floatexp(height)), floatexp(pixel_spacing), floatexp(precision), &progress[2], running); break;
      }
    }
  }
  nt_current = nt;
  *ended = true;
}

void subframe_thread(map &out, stats &sta, const formula *form, const param &par, const count_t subframe, progress_t *progress, bool *running, bool *ended)
{
  if (form->complex_analytic())
  {
    const formulaCbase *fc = static_cast<const formulaCbase *>(form);
    switch (nt_current)
    {
      case nt_float:
        fc->render(out, sta, BCf, subframe, par, float(par.Zoom), Zf.size(), &Zf[0], progress, running);
        break;
      case nt_double:
        fc->render(out, sta, BCd, subframe, par, double(par.Zoom), Zd.size(), &Zd[0], progress, running);
        break;
      case nt_longdouble:
        fc->render(out, sta, BCld, subframe, par, (long double)(par.Zoom), Zld.size(), &Zld[0], progress, running);
        break;
      case nt_floatexp:
        fc->render(out, sta, BCfe, subframe, par, par.Zoom, Zfe.size(), &Zfe[0], progress, running);
        break;
    }
  }
  else
  {
    const formulaR2base *fr2 = static_cast<const formulaR2base *>(form);
    switch (nt_current)
    {
      case nt_float:
        fr2->render(out, sta, BR2f, subframe, par, float(par.Zoom), Zf.size(), &Zf[0], progress, running);
        break;
      case nt_double:
        fr2->render(out, sta, BR2d, subframe, par, double(par.Zoom), Zd.size(), &Zd[0], progress, running);
        break;
      case nt_longdouble:
        fr2->render(out, sta, BR2ld, subframe, par, (long double)(par.Zoom), Zld.size(), &Zld[0], progress, running);
        break;
      case nt_floatexp:
        fr2->render(out, sta, BR2fe, subframe, par, par.Zoom, Zfe.size(), &Zfe[0], progress, running);
        break;
    }
  }
  *ended = true;
}
