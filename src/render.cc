// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <cmath>
#include <iostream>

#include "bla.h"
#include "map.h"
#include "param.h"
#include "render.h"

template <typename real>
void render(map &out, const param &par, const real Zoom, const count_t M, const complex<real> *Zp, progress_t *progress, bool *running)
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
  const count_t PerturbIterations = par.PerturbIterations;
  // initialize table
  const real ER2 = 65536.0 * 65536.0;
  const real pixel_spacing = 4 / Zoom / height;
  const real step_count = 1000; // FIXME TODO
  const blas<real> BLA(M, Zp, hypot(width, height), pixel_spacing, step_count, &progress[0], running);
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
  count_t total_bla_steps = 0;
  count_t total_bla_iterations = 0;
  count_t total_perturb_iterations = 0;
  count_t total_rebases = 0;
  count_t total_iterations = 0;
  count_t total_pixels = 0;
  count_t minimum_iterations = 0x7fffFFFFffffFFFFLL;;
  count_t maximum_iterations = 0;
  #pragma omp parallel for reduction(+:total_bla_steps) reduction(+:total_bla_iterations) reduction(+:total_perturb_iterations) reduction(+:total_iterations) reduction(+:total_rebases) reduction(min:minimum_iterations) reduction(max:maximum_iterations)
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
    const complex<real> c (cx, cy);

    count_t m = 0;
    count_t n = 0;
    complex<real> Z (Zp[0]);
    complex<real> z (0);
    real z2 (norm(z));
    complex<real> dZdC (0);
    complex<real> Zz (Z + z);
    real Zz2 (norm(Zz));
#ifdef VERBOSE
    if (i == 0 && j == 0)
    {
      std::cerr << (2 / Zoom) << " ";
    }
#endif

    while (n < Iterations && Zz2 < ER2 && perturb_iterations < PerturbIterations)
    {
      // bla steps
      const bla<real> *b = 0;
      while (n < Iterations && Zz2 < ER2 && (b = BLA.lookup(m, z2)))
      {
        const complex<real> A = b->A;
        const complex<real> B = b->B;
        count_t l = b->l;
        const complex<real> zn = A * z + B * c;
        const complex<real> dZdCn = A * dZdC + B * pixel_spacing;
        z = zn;
        z2 = norm(z);
        dZdC = dZdCn;
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
#ifdef VERBOSE
        if (i == 0 && j == 0)
        {
          std::cerr << l << " ";
        }
#endif

        // rebase
        assert(m < M);
        complex<real> Z = Zp[m];
        Zz = Z + z;
        Zz2 = norm(Zz);
        if (Zz2 < z2 || (ReferencePeriod == 0 && m == M - 1))
        {
          z = Zz;
          m = 0;
          rebases++;
#ifdef VERBOSE
          if (i == 0 && j == 0)
          {
            std::cerr << "B ";
          }
#endif
        }
      }

      // perturbation iteration
      {
        assert(m < M);
        complex<real> Z = Zp[m];
        Zz = Z + z;
        // z = (2 Z + z) z + c
        complex<real> ZZz = 2 * Z + z;
        complex<real> zn = ZZz * z + c;
        complex<real> dZdCn = 2 * dZdC * Zz + pixel_spacing;
        z = zn;
        z2 = norm(z);
        dZdC = dZdCn;
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
#ifdef VERBOSE
        if (i == 0 && j == 0)
        {
          std::cerr << "p ";
        }
#endif
      }

      {
        // rebase
        assert(m < M);
        complex<real> Z = Zp[m];
        Zz = Z + z;
        Zz2 = norm(Zz);
        if (Zz2 < z2 || (ReferencePeriod == 0 && m == M - 1))
        {
          z = Zz;
          m = 0;
          rebases++;
#ifdef VERBOSE
          if (i == 0 && j == 0)
          {
            std::cerr << "b ";
          }
#endif
        }
      }
    }

#ifdef VERBOSE
    // output
    if (i == 0 && j == 0)
    {
      std::cerr << n << " ]" << std::endl;
      std::cerr << "estimated speedup: " << (n * 1.0 / (n - bla_iterations)) << " (" << rebases << ")" << std::endl;
    }
#endif

    // compute colour
    complex<float> Z1 = complex<float>(float(Zz.x), float(Zz.y));
    complex<float> dC = complex<float>(float(dZdC.x), float(dZdC.y));
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
    out.setDE(i, j, complex<float>(float(de.x), float(de.y)));

    // accumulate statistics
    maximum_iterations = maximum_iterations > n ? maximum_iterations : n;
    minimum_iterations = minimum_iterations < n ? minimum_iterations : n;
    total_bla_iterations += bla_iterations;
    total_bla_steps += bla_steps;
    total_iterations += n;
    total_perturb_iterations += perturb_iterations;
    total_rebases += rebases;
    count_t count;
    #pragma omp atomic capture
    count = ++total_pixels;
    progress[1] = count / progress_t(width * height);
  }

#ifdef VERBOSE
  // output statistics
  std::cerr << minimum_iterations << " minimum iterations" << std::endl;
  std::cerr << maximum_iterations << " minimum iterations" << std::endl;
  std::cerr << (total_iterations / (double) total_pixels) << " average iterations" << std::endl;
  std::cerr << (total_bla_iterations / (double) total_pixels) << " average bla iterations" << std::endl;
  std::cerr << (total_perturb_iterations / (double) total_pixels) << " average ptb iterations" << std::endl;
  std::cerr << ((total_perturb_iterations + total_bla_steps) / (double) total_pixels) << " average steps" << std::endl;
  std::cerr << (total_bla_steps / (double) total_pixels) << " average bla steps" << std::endl;
  std::cerr << (total_perturb_iterations / (double) total_pixels) << " average ptb steps" << std::endl;
  std::cerr << (total_bla_iterations / (double) total_bla_steps) << " iterations per bla" << std::endl;
  std::cerr << (total_rebases / (double) total_pixels) << " rebases" << std::endl;
  std::cerr << "speedup: " << (total_iterations / (double) (total_perturb_iterations + total_bla_steps)) << "x" << std::endl;
#endif
}

template void render(map &out, const param &par, const float Zoom, const count_t M, const complex<float> *Zp, progress_t *progress, bool *running);
template void render(map &out, const param &par, const double Zoom, const count_t M, const complex<double> *Zp, progress_t *progress, bool *running);
template void render(map &out, const param &par, const long double Zoom, const count_t M, const complex<long double> *Zp, progress_t *progress, bool *running);
template void render(map &out, const param &par, const floatexp Zoom, const count_t M, const complex<floatexp> *Zp, progress_t *progress, bool *running);
