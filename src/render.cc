// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <cmath>
#include <iostream>

#include "bla.h"
#include "complex.h"
#include "dual.h"
#include "floatexp.h"
#include "formula.h"
#include "map.h"
#include "param.h"
#include "render.h"
#include "stats.h"

template <typename real>
void render(map &out, stats &sta, const param &par, const real Zoom, const count_t M, const complex<real> *Zp, const formulaC *form, progress_t *progress, bool *running)
{
  using std::isinf;
  using std::isnan;
  using std::log;
  using std::max;
  using std::min;
  const coord_t width = out.width;
  const coord_t height = out.height;
  const count_t Iterations = par.Iterations;
  const count_t ReferencePeriod = par.ReferencePeriod;
  const count_t PerturbIterations = par.MaxPtbIters;
  // initialize table
  const real ER2 = par.EscapeRadius * par.EscapeRadius;
  const real pixel_spacing = 4 / Zoom / height;
  const real step_count = 1000; // FIXME TODO
  const blasC<real> BLA(M, Zp, form, hypot(width, height), pixel_spacing, step_count, &progress[0], running);
  if (! *running)
  {
    return;
  }
#ifdef VERBOSE
  for (count_t level = 0; level < BLA.L; ++level)
  {
    std::cerr << BLA.b[level][0].l << "\t" << sqrt(BLA.b[level][0].r2) << std::endl;
  }
#endif
  count_t minimum_iterations = sta.minimum_iterations;
  count_t maximum_iterations = sta.maximum_iterations;
  #pragma omp parallel for reduction(min:minimum_iterations) reduction(max:maximum_iterations)
  for (coord_t j = 0; j < height; ++j) if (*running)
  for (coord_t i = 0; i < width; ++i) if (*running)
  {
    count_t bla_steps = 0;
    count_t bla_iterations = 0;
    count_t perturb_iterations = 0;
    count_t rebases = 0;
    // FIXME TODO ExponentialMap
    const real cx = real(((i + 0.5) / width - 0.5) * width) * pixel_spacing;
    const real cy = real((0.5 - (j + 0.5) / height) * height) * pixel_spacing;
    const complex<real> C (Zp[1]);
    dual<1, complex<real>> c (complex<real>(cx, cy));
    c.dx[0] = complex<real>(pixel_spacing);
    count_t m = 0;
    count_t n = 0;
    complex<real> Z (Zp[0]);
    dual<1, complex<real>> z (0);
    real z2 (norm(z.x));
    dual<1, complex<real>> Zz (Z + z);
    real Zz2 (norm(Zz.x));

    while (n < Iterations && Zz2 < ER2 && perturb_iterations < PerturbIterations)
    {
      // bla steps
      const blaC<real> *b = 0;
      while (n < Iterations && Zz2 < ER2 && (b = BLA.lookup(m, z2)))
      {
        const complex<real> A = b->A;
        const complex<real> B = b->B;
        count_t l = b->l;
        z = A * z + B * c;
        z2 = norm(z.x);
        n += l;
        m += l;
        bla_steps++;
        bla_iterations += l;
        if (ReferencePeriod)
        {
          while (m >= ReferencePeriod)
          {
            m -= ReferencePeriod;
          }
        }

        // rebase
        if (! (n < Iterations && Zz2 < ER2 && perturb_iterations < PerturbIterations))
        {
          break;
        }
        if (! (m < M))
        {
          break;
        }
        complex<real> Z = Zp[m];
        Zz = Z + z;
        Zz2 = norm(Zz.x);
        if (Zz2 < z2 || (ReferencePeriod == 0 && m == M - 1))
        {
          z = Zz;
          m = 0;
          rebases++;
        }
      }

      // perturbation iteration
      {
        if (! (n < Iterations && Zz2 < ER2 && perturb_iterations < PerturbIterations))
        {
          break;
        }
        if (! (m < M))
        {
          break;
        }
        complex<real> Z = Zp[m];
        Zz = Z + z;
        // z = (2 Z + z) z + c
        z = form->perturb(C, Zp[m], c, z);
        z2 = norm(z.x);
        n++;
        m++;
        perturb_iterations++;
        if (ReferencePeriod)
        {
          while (m >= ReferencePeriod)
          {
            m -= ReferencePeriod;
          }
        }
      }

      {
        // rebase
        if (! (n < Iterations && Zz2 < ER2 && perturb_iterations < PerturbIterations))
        {
          break;
        }
        if (! (m < M))
        {
          break;
        }
        complex<real> Z = Zp[m];
        Zz = Z + z;
        Zz2 = norm(Zz.x);
        if (Zz2 < z2 || (ReferencePeriod == 0 && m == M - 1))
        {
          z = Zz;
          m = 0;
          rebases++;
        }
      }
    }

    // compute output
    complex<float> Z1 = complex<float>(float(Zz.x.x), float(Zz.x.y));
    complex<float> dC = complex<float>(float(Zz.dx[0].x), float(Zz.dx[0].y));
    complex<float> de = abs(Z1) * log(abs(Z1)) / dC;
    float nf = 0; // FIXME TODO
    float t = arg(Z1) / (2 * M_PI);
    t -= floor(t);
    if (Zz2 < ER2 || isnan(de.x) || isinf(de.x) || isnan(de.y) || isinf(de.y))
    {
      n = Iterations;
      nf = 0;
      t = 0;
      de = 0;
    }
    out.setN(i, j, n);
    out.setNF(i, j, nf);
    out.setT(i, j, t);
    out.setDE(i, j, de);

    // accumulate statistics
    if (n < Iterations)
    {
      maximum_iterations = maximum_iterations > n ? maximum_iterations : n;
    }
    minimum_iterations = minimum_iterations < n ? minimum_iterations : n;
    count_t count;
    #pragma omp atomic capture
    count = ++sta.pixels;
    #pragma omp atomic
    sta.bla_iterations += bla_iterations;
    #pragma omp atomic
    sta.bla_steps += bla_steps;
    #pragma omp atomic
    sta.iterations += n;
    #pragma omp atomic
    sta.perturb_iterations += perturb_iterations;
    #pragma omp atomic
    sta.rebases += rebases;
    progress[1] = count / progress_t(width * height);
  }
  sta.minimum_iterations = minimum_iterations;
  sta.maximum_iterations = maximum_iterations;
}

template void render(map &out, stats &sta, const param &par, const float Zoom, const count_t M, const complex<float> *Zp, const formulaC *form, progress_t *progress, bool *running);
template void render(map &out, stats &sta, const param &par, const double Zoom, const count_t M, const complex<double> *Zp, const formulaC *form, progress_t *progress, bool *running);
template void render(map &out, stats &sta, const param &par, const long double Zoom, const count_t M, const complex<long double> *Zp, const formulaC *form, progress_t *progress, bool *running);
template void render(map &out, stats &sta, const param &par, const floatexp Zoom, const count_t M, const complex<floatexp> *Zp, const formulaC *form, progress_t *progress, bool *running);

template <typename real>
void render(map &out, stats &sta, const param &par, const real Zoom, const count_t M, const complex<real> *Zp, const formulaR2 *form, progress_t *progress, bool *running)
{
  using std::isinf;
  using std::isnan;
  using std::log;
  using std::max;
  using std::min;
  const coord_t width = out.width;
  const coord_t height = out.height;
  const count_t Iterations = par.Iterations;
  const count_t ReferencePeriod = par.ReferencePeriod;
  const count_t PerturbIterations = par.MaxPtbIters;
  // initialize table
  const real ER2 = par.EscapeRadius * par.EscapeRadius;
  const real pixel_spacing = 4 / Zoom / height;
  const real step_count = 1000; // FIXME TODO
  const blasR2<real> BLA(M, Zp, form, hypot(width, height), pixel_spacing, step_count, &progress[0], running);
  if (! *running)
  {
    return;
  }
#ifdef VERBOSE
  for (count_t level = 0; level < BLA.L; ++level)
  {
    std::cerr << BLA.b[level][0].l << "\t" << sqrt(BLA.b[level][0].r2) << std::endl;
  }
#endif
  count_t minimum_iterations = sta.minimum_iterations;
  count_t maximum_iterations = sta.maximum_iterations;
  #pragma omp parallel for reduction(min:minimum_iterations) reduction(max:maximum_iterations)
  for (coord_t j = 0; j < height; ++j) if (*running)
  for (coord_t i = 0; i < width; ++i) if (*running)
  {
    count_t bla_steps = 0;
    count_t bla_iterations = 0;
    count_t perturb_iterations = 0;
    count_t rebases = 0;
    // FIXME TODO ExponentialMap
    dual<2, real> cx (real(((i + 0.5) / width - 0.5) * width) * pixel_spacing);
    cx.dx[0] = pixel_spacing;
    dual<2, real> cy (real((0.5 - (j + 0.5) / height) * height) * pixel_spacing);
    cy.dx[1] = pixel_spacing;
    const complex<real> C (Zp[1]);
    const complex<dual<2, real>> c (cx, cy);
    count_t m = 0;
    count_t n = 0;
    complex<real> Z (Zp[0]);
    complex<dual<2, real>> z (0);
    real z2 (norm(z).x);
    complex<dual<2, real>> Zz (Z + z);
    real Zz2 (norm(Zz).x);

    while (n < Iterations && Zz2 < ER2 && perturb_iterations < PerturbIterations)
    {
      // bla steps
      const blaR2<real> *b = 0;
      while (n < Iterations && Zz2 < ER2 && (b = BLA.lookup(m, z2)))
      {
        const mat2<real> A = b->A;
        const mat2<real> B = b->B;
        count_t l = b->l;
        z = A * z + B * c;
        z2 = norm(z).x;
        n += l;
        m += l;
        bla_steps++;
        bla_iterations += l;
        if (ReferencePeriod)
        {
          while (m >= ReferencePeriod)
          {
            m -= ReferencePeriod;
          }
        }

        // rebase
        if (! (n < Iterations && Zz2 < ER2 && perturb_iterations < PerturbIterations))
        {
          break;
        }
        if (! (m < M))
        {
          break;
        }
        complex<real> Z = Zp[m];
        Zz = Z + z;
        Zz2 = norm(Zz).x;
        if (Zz2 < z2 || (ReferencePeriod == 0 && m == M - 1))
        {
          z = Zz;
          m = 0;
          rebases++;
        }
      }

      // perturbation iteration
      {
        if (! (n < Iterations && Zz2 < ER2 && perturb_iterations < PerturbIterations))
        {
          break;
        }
        if (! (m < M))
        {
          break;
        }
        complex<real> Z = Zp[m];
        Zz = Z + z;
        // z = (2 Z + z) z + c
        z = form->perturb(C, Zp[m], c, z);
        z2 = norm(z).x;
        n++;
        m++;
        perturb_iterations++;
        if (ReferencePeriod)
        {
          while (m >= ReferencePeriod)
          {
            m -= ReferencePeriod;
          }
        }
      }

      {
        // rebase
        if (! (n < Iterations && Zz2 < ER2 && perturb_iterations < PerturbIterations))
        {
          break;
        }
        if (! (m < M))
        {
          break;
        }
        complex<real> Z = Zp[m];
        Zz = Z + z;
        Zz2 = norm(Zz).x;
        if (Zz2 < z2 || (ReferencePeriod == 0 && m == M - 1))
        {
          z = Zz;
          m = 0;
          rebases++;
        }
      }
    }

    // compute output
    complex<float> Z1 = complex<float>(float(Zz.x.x), float(Zz.y.x));
    complex<float> dC = complex<float>(float(Zz.x.dx[0]), float(Zz.y.dx[1])); // FIXME
    complex<float> de = abs(Z1) * log(abs(Z1)) / dC;
    float nf = 0; // FIXME TODO
    float t = arg(Z1) / (2 * M_PI);
    t -= floor(t);
    if (Zz2 < ER2 || isnan(de.x) || isinf(de.x) || isnan(de.y) || isinf(de.y))
    {
      n = Iterations;
      nf = 0;
      t = 0;
      de = 0;
    }
    out.setN(i, j, n);
    out.setNF(i, j, nf);
    out.setT(i, j, t);
    out.setDE(i, j, de);

    // accumulate statistics
    if (n < Iterations)
    {
      maximum_iterations = maximum_iterations > n ? maximum_iterations : n;
    }
    minimum_iterations = minimum_iterations < n ? minimum_iterations : n;
    count_t count;
    #pragma omp atomic capture
    count = ++sta.pixels;
    #pragma omp atomic
    sta.bla_iterations += bla_iterations;
    #pragma omp atomic
    sta.bla_steps += bla_steps;
    #pragma omp atomic
    sta.iterations += n;
    #pragma omp atomic
    sta.perturb_iterations += perturb_iterations;
    #pragma omp atomic
    sta.rebases += rebases;
    progress[1] = count / progress_t(width * height);
  }
  sta.minimum_iterations = minimum_iterations;
  sta.maximum_iterations = maximum_iterations;
}

template void render(map &out, stats &sta, const param &par, const float Zoom, const count_t M, const complex<float> *Zp, const formulaR2 *form, progress_t *progress, bool *running);
template void render(map &out, stats &sta, const param &par, const double Zoom, const count_t M, const complex<double> *Zp, const formulaR2 *form, progress_t *progress, bool *running);
template void render(map &out, stats &sta, const param &par, const long double Zoom, const count_t M, const complex<long double> *Zp, const formulaR2 *form, progress_t *progress, bool *running);
template void render(map &out, stats &sta, const param &par, const floatexp Zoom, const count_t M, const complex<floatexp> *Zp, const formulaR2 *form, progress_t *progress, bool *running);

template <typename real>
void render(map &out, stats &sta, const param &par, const real Zoom, const count_t M, const complex<real> *Zp, const formula *form, progress_t *progress, bool *running)
{
  if (form->complex_analytic())
  {
    render(out, sta, par, Zoom, M, Zp, dynamic_cast<const formulaC *>(form), progress, running);
  }
  else
  {
    render(out, sta, par, Zoom, M, Zp, dynamic_cast<const formulaR2 *>(form), progress, running);
  }
}

template void render(map &out, stats &sta, const param &par, const float Zoom, const count_t M, const complex<float> *Zp, const formula *form, progress_t *progress, bool *running);
template void render(map &out, stats &sta, const param &par, const double Zoom, const count_t M, const complex<double> *Zp, const formula *form, progress_t *progress, bool *running);
template void render(map &out, stats &sta, const param &par, const long double Zoom, const count_t M, const complex<long double> *Zp, const formula *form, progress_t *progress, bool *running);
template void render(map &out, stats &sta, const param &par, const floatexp Zoom, const count_t M, const complex<floatexp> *Zp, const formula *form, progress_t *progress, bool *running);
