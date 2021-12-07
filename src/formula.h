// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <vector>

#include <mpreal.h>

using mpreal = mpfr::mpreal;

#include "complex.h"
#include "dual.h"
#include "floatexp.h"
#include "map.h"
#include "matrix.h"
#include "param.h"
#include "stats.h"
#include "types.h"

template
  < typename t
  , complex<mpreal> PLAIN(const complex<mpreal>&, const complex<mpreal>&)
  >
count_t reference(complex<t> *Zp, const count_t &MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running)
{
  complex<mpreal> Z (0);
  count_t M = MaxRefIters - 1;
  // calculate reference in high precision
  for (count_t i = 0; i < MaxRefIters; ++i)
  {
    // store low precision orbit
    Zp[i] = complex<t>(t(Z.x), t(Z.y));
    // escape check
    if (norm(Zp[i]) > 4 || ! *running) // FIXME escape radius
    {
      M = i;
      break;
    }
    // step
    Z = PLAIN(C, Z);
    *progress = (i + 1) / progress_t(MaxRefIters);
  }
  return M;
}

template
  < typename t
  , dual<1, complex<t>> PERTURB(const complex<t>&, const complex<t>&, const dual<1, complex<t>>&, const dual<1, complex<t>>&) noexcept
  >
bool periodC(count_t &period, const complex<t> *Zp, const count_t M, const complex<t> c0, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) noexcept
{
  const complex<t> C (Zp[1]); // FIXME
  dual<1, complex<t>> c (c0);
  c.dx[0] = 1;
  dual<1, complex<t>> z (0);
  mat2<double> K1(inverse(K));
  t z2 = 0;
  t r2 = 1e16;
  bool p = true;
  count_t i = 0;
  count_t m = 0;
  while (i < N && z2 < r2 && p && *running)
  {
    // progress
    progress[0] = i / progress_t(N);
    // formula
    z = PERTURB(C, Zp[m], c, z);
    m++;
    // rebase
    const complex<t> Z = Zp[m];
    const dual<1, complex<t>> Zz = Z + z;
    const t Zz2 = norm(Zz.x);
    if (Zz2 < z2 || m == M - 1)
    {
      z = Zz;
      m = 0;
    }
    // (u v) = s^{-1} K^{-1} J^{-1} (x y)
    complex<t> w = (K1 * (z.x / z.dx[0]));
    p = 1 <= floatexp(norm(w)) / (s * s);
    ++i;
    z2 = norm(Zp[m] + z.x);
  }
  if (i == N || r2 <= z2 || p || ! *running)
  {
    return false;
  }
  period = i;
  return true;
}

template
  < typename t
  , complex<dual<2, t>> PERTURB(const complex<t>&, const complex<t>&, const complex<dual<2, t>>&, const complex<dual<2, t>>&) noexcept
  >
bool periodR2(count_t &period, const complex<t> *Zp, const count_t &M, const complex<t> &c0, const count_t &N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) noexcept
{
  complex<t> C (Zp[1]); // FIXME
  dual<2, t> cx (c0.x); cx.dx[0] = 1;
  dual<2, t> cy (c0.y); cy.dx[1] = 1;
  const complex<dual<2, t>> c (cx, cy);
  complex<dual<2, t>> z (0);
  mat2<double> K1 (inverse(K));
  t z2 = 0;
  t r2 = 1e16;
  bool p = true;
  count_t i = 0;
  count_t m = 0;
  while (i < N && z2 < r2 && p && *running)
  {
    // progress
    progress[0] = i / progress_t(N);
    // formula
    z = PERTURB(C, Zp[m], c, z);
    m++;
    // rebase
    const complex<t> Z = Zp[m];
    const complex<dual<2, t>> Zz = Z + z;
    const t Zz2 = norm(complex<t>(Zz.x.x, Zz.y.x));
    if (Zz2 < z2 || m == M - 1)
    {
      z = Zz;
      m = 0;
    }
    // (u1 v1) = s^{-1} K^{-1} J^{-1} (u0 v0)
    const mat2<t> J(z.x.dx[0], z.x.dx[1], z.y.dx[0], z.y.dx[1]);
    complex<t> w = (K1 * (inverse(J) * complex<t>(z.x.x, z.y.x)));
    p = 1 <= floatexp(norm(w)) / (s * s);
    ++i;
    z2 = norm(Zp[m] + complex<t>(z.x.x, z.y.x));
  }
  if (i == N || r2 <= z2 || p || ! *running)
  {
    return false;
  }
  period = i;
  return true;
}

template < dual<1, complex<mpreal>> PLAIN(const dual<1, complex<mpreal>>&, const dual<1, complex<mpreal>>&) >
bool centerC(complex<mpreal> &C0, const count_t period, progress_t *progress, bool *running)
{
  mpfr_prec_t prec = std::max(mpfr_get_prec(C0.x.mpfr_srcptr()), mpfr_get_prec(C0.y.mpfr_srcptr()));
  const floatexp epsilon2 = floatexp(1, 16 - 2 * prec);
  double lepsilon2 = double(log(epsilon2));
  double ldelta0 = 0;
  double ldelta1 = 0;
  progress_t eta = 0;
  bool converged = false;
  const count_t maxsteps = 64;
  for (count_t j = 0; j < maxsteps && *running && ! converged; ++j)
  {
    progress[0] = j / eta;
    progress[1] = 0;
    dual<1, complex<mpreal>> c(C0); c.dx[0] = 1;
    dual<1, complex<mpreal>> z(0);
    // iteration
    for (count_t i = 0; i < period && *running; ++i)
    {
      progress[1] = i / progress_t(period);
      z = PLAIN(c, z);
    }
    if (*running)
    {
      // Newton step
      const complex<mpreal> u = - z.x / z.dx[0];
      C0 += u;
      // check convergence
      floatexp uf = floatexp(u.x);
      floatexp vf = floatexp(u.y);
      floatexp delta = sqr(uf) + sqr(vf);
      converged = delta < epsilon2;
      ldelta0 = ldelta1;
      ldelta1 = double(log(delta));
      eta = log2((lepsilon2 - ldelta0) / (ldelta1 - ldelta0));
    }
  }
  return converged;
}

template < complex<dual<2, mpreal>> PLAIN(const complex<dual<2, mpreal>>&, const complex<dual<2, mpreal>>&) >
bool centerR2(complex<mpreal> &C0, const count_t period, progress_t *progress, bool *running)
{
  mpfr_prec_t prec = std::max(mpfr_get_prec(C0.x.mpfr_srcptr()), mpfr_get_prec(C0.y.mpfr_srcptr()));
  const floatexp epsilon2 = floatexp(1, 16 - 2 * prec);
  double lepsilon2 = double(log(epsilon2));
  double ldelta0 = 0;
  double ldelta1 = 0;
  progress_t eta = 0;
  bool converged = false;
  const count_t maxsteps = 64;
  for (count_t j = 0; j < maxsteps && *running && ! converged; ++j)
  {
    progress[0] = j / eta;
    progress[1] = 0;
    dual<2, mpreal> cx(C0.x); cx.dx[0] = 1;
    dual<2, mpreal> cy(C0.y); cy.dx[1] = 1;
    complex<dual<2, mpreal>> c(cx, cy);
    complex<dual<2, mpreal>> z(0, 0);
    // iteration
    for (count_t i = 0; i < period && *running; ++i)
    {
      progress[1] = i / progress_t(period);
      z = PLAIN(c, z);
    }
    if (*running)
    {
      const mpreal &x = z.x.x;
      const mpreal &y = z.y.x;
      const mpreal &dxa = z.x.dx[0];
      const mpreal &dxb = z.x.dx[1];
      const mpreal &dya = z.y.dx[0];
      const mpreal &dyb = z.y.dx[1];
      // Newton step
      const mpreal det = dxa * dyb - dxb * dya;
      const mpreal u = -( dyb * x - dxb * y) / det;
      const mpreal v = -(-dya * x + dxa * y) / det;
      C0.x += u;
      C0.y += v;
      // check convergence
      floatexp uf = floatexp(u);
      floatexp vf = floatexp(v);
      floatexp delta = sqr(uf) + sqr(vf);
      converged = delta < epsilon2;
      ldelta0 = ldelta1;
      ldelta1 = double(log(delta));
      eta = log2((lepsilon2 - ldelta0) / (ldelta1 - ldelta0));
    }
  }
  return converged;
}

template
  < typename t
  , dual<1, complex<t>> PERTURB(const complex<t>&, const complex<t>&, const dual<1, complex<t>>&, const dual<1, complex<t>>&) noexcept
  >
bool sizeC(floatexp &s, mat2<double> &K, const complex<t> *Zp, count_t period, const complex<t> &c0, progress_t *progress, bool *running) noexcept
{
  using std::abs;
  using std::exp;
  using std::log;
  const double degree = 2; // FIXME
  const complex<t> C (Zp[1]); // FIXME
  dual<1, complex<t>> c (c0);
  dual<1, complex<t>> z (c);
  z.dx[0] = 1;
  complex<t> b (1);
  count_t j = 1;
  count_t m = 1;
  if (m == period)
  {
    m = 0;
  }
  while (j < period && *running)
  {
    progress[0] = j / progress_t(period);
    z = PERTURB(C, Zp[m], c, z);
    m++;
    if (m == period)
    {
      m = 0;
    }
    // rebase
    const complex<t> Z = Zp[m];
    const dual<1, complex<t>> Zz = Z + z;
    const t z2 = norm(z.x);
    const t Zz2 = norm(Zz.x);
    if (Zz2 < z2)
    {
      z = Zz;
      m = 0;
    }
    b += 1 / z.dx[0];
    ++j;
  }
  // l^d b
  if (*running)
  {
    double d = degree / (degree - 1);
    const t l = abs(z.dx[0]);
    const t beta = abs(b);
    const t llb = exp(log(l) * d) * beta;
    s = floatexp(1 / llb);
    K = mat2<double>(1);
    return true;
  }
  return false;
}

template
  < typename t
  , complex<dual<2, t>> PERTURB(const complex<t>&, const complex<t>&, const complex<dual<2, t>>&, const complex<dual<2, t>>&) noexcept
  >
bool sizeR2(floatexp &s, mat2<double> &K, const complex<t> *Zp, count_t period, const complex<t> &c0, progress_t *progress, bool *running) noexcept
{
  using std::abs;
  using std::exp;
  using std::log;
  using std::sqrt;
  const double degree = 2; // FIXME
  const complex<t> C (Zp[1]); // FIXME
  const complex<dual<2, t>> c (dual<2, t>(c0.x), dual<2, t>(c0.y));
  complex<dual<2, t>> z(c);
  z.x.dx[0] = 1;
  z.y.dx[1] = 1;
  mat2<t> b (1);
  count_t j = 1;
  count_t m = 1;
  if (m == period)
  {
    m = 0;
  }
  while (j < period && *running)
  {
    progress[0] = j / progress_t(period);
    z = PERTURB(C, Zp[m], c, z);
    m++;
    if (m == period)
    {
      m = 0;
    }
    // rebase
    const complex<t> Z = Zp[m];
    const complex<dual<2, t>> Zz = Z + z;
    const t z2 = norm(complex<t>(z.x.x, z.y.x));
    const t Zz2 = norm(complex<t>(Zz.x.x, Zz.y.x));
    if (Zz2 < z2)
    {
      z = Zz;
      m = 0;
    }
    mat2<t> l (z.x.dx[0], z.x.dx[1], z.y.dx[0], z.y.dx[1]);
    b += inverse(l);
    ++j;
  }
  // l^d b
  if (*running)
  {
    double d = degree / (degree - 1);
    mat2<t> l (z.x.dx[0], z.x.dx[1], z.y.dx[0], z.y.dx[1]);
    const t lambda = sqrt(abs(determinant(l)));
    const t beta = sqrt(abs(determinant(b)));
    const t llb = exp(log(lambda) * d) * beta;
    s = floatexp(1 / llb);
    b = inverse(transpose(b)) / beta;
    K = mat2<double>(double(b.x[0][0]), double(b.x[0][1]), double(b.x[1][0]), double(b.x[1][1]));
    return true;
  }
  return false;
}

template
  < typename t
  , dual<1, complex<t>> PERTURB(const complex<t>&, const complex<t>&, const dual<1, complex<t>>&, const dual<1, complex<t>>&) noexcept
  >
bool domain_sizeC(floatexp &s, const complex<t> *Zp, count_t period, const complex<t> &c0, progress_t *progress, bool *running) noexcept
{
  const complex<t> C (Zp[1]); // FIXME
  dual<1, complex<t>> c (c0);
  c.dx[0] = 1;
  dual<1, complex<t>> z (c);
  count_t j = 2;
  count_t m = 1;
  if (m == period)
  {
    m = 0;
  }
  t zq2 = norm(z.x);
  while (j <= period && *running)
  {
    // progress
    progress[0] = j / progress_t(period);
    // formula
    z = PERTURB(C, Zp[m], c, z);
    m++;
    if (m == period)
    {
      m = 0;
    }
    // rebase
    const complex<t> Z = Zp[m];
    const dual<1, complex<t>> Zz = Z + z;
    const t z2 = norm(z.x);
    const t Zz2 = norm(Zz.x);
    if (Zz2 < z2)
    {
      z = Zz;
      m = 0;
    }
    // capture penultimate minimum |z|
    t zp2 = norm(z.x);
    if (j < period && zp2 < zq2)
    {
      zq2 = zp2;
    }
    ++j;
  }
  if (*running)
  {
    s = sqrt(zq2) / abs(z.dx[0]);
    return true;
  }
  return false;
}

template
  < typename t
  , complex<dual<2, t>> PERTURB(const complex<t>&, const complex<t>&, const complex<dual<2, t>>&, const complex<dual<2, t>>&) noexcept
  >
bool domain_sizeR2(floatexp &s, const complex<t> *Zp, count_t period, const complex<t> &c0, progress_t *progress, bool *running) noexcept
{
  const complex<t> C (Zp[1]); // FIXME
  dual<2, t> cx (c0.x);
  cx.dx[0] = 1;
  dual<2, t> cy (c0.y);
  cy.dx[1] = 1;
  const complex<dual<2, t>> c (cx, cy);
  complex<dual<2, t>> z(c);
  count_t j = 2;
  count_t m = 1;
  if (m == period)
  {
    m = 0;
  }
  t zq2 = norm(Zp[m] + complex<t>(z.x.x, z.y.x));
  while (j <= period && *running)
  {
    // progress
    progress[0] = j / progress_t(period);
    // formula
    z = PERTURB(C, Zp[m], c, z);
    m++;
    if (m == period)
    {
      m = 0;
    }
    // rebase
    const complex<t> Z = Zp[m];
    const complex<dual<2, t>> Zz = Z + z;
    const t z2 = norm(complex<t>(z.x.x, z.y.x));
    const t Zz2 = norm(complex<t>(Zz.x.x, Zz.y.x));
    if (Zz2 < z2)
    {
      z = Zz;
      m = 0;
    }
    // capture penultimate minimum |z|
    t zp2 = norm(Zp[m] + complex<t>(z.x.x, z.y.x));
    if (j < period && zp2 < zq2)
    {
      zq2 = zp2;
    }
    ++j;
  }
  if (*running)
  {
    mat2<t> L (z.x.dx[0], z.x.dx[1], z.y.dx[0], z.y.dx[1]);
    s = sqrt(zq2) / sqrt(abs(determinant(L)));
    return true;
  }
  return false;
}

template
  < typename t
  , blaC<t> BLA(const t &h, const t &k, const t &L, const complex<t> &Z) noexcept
  >
void blas_init1C(blasC<t> *Bp, const complex<t> *Zp, const t h, const t k, t L, progress_t *progress, bool *running) noexcept
{
  using std::max;
  const count_t M = Bp->M;
  count_t total = 0;
  #pragma omp parallel for
  for (count_t m = 1; m < M; ++m) if (running)
  {
    Bp->b[0][m - 1] = BLA(h, k, L, Zp[m]);
    count_t done;
    #pragma omp atomic capture
    done = total++;
    progress[0] = done / progress_t(2 * M);
  }
}

template
  < typename t
  , blaR2<t> BLA(const t &h, const t &k, const t &L, const complex<t> &Z) noexcept
  >
void blas_init1R2(blasR2<t> *Bp, const complex<t> *Zp, const t h, const t k, t L, progress_t *progress, bool *running) noexcept
{
  using std::max;
  const count_t M = Bp->M;
  count_t total = 0;
  #pragma omp parallel for
  for (count_t m = 1; m < M; ++m) if (running)
  {
    Bp->b[0][m - 1] = BLA(h, k, L, Zp[m]);
    count_t done;
    #pragma omp atomic capture
    done = total++;
    progress[0] = done / progress_t(2 * M);
  }
}

template
  < typename real
  , dual<1, complex<real>> PERTURB(const complex<real> &, const complex<real> &, const dual<1, complex<real>> &, const dual<1, complex<real>> &)
  >
void renderC(map &out, stats &sta, const param &par, const real Zoom, const count_t M, const complex<real> *Zp, const formulaCbase *form, progress_t *progress, bool *running)
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
  const count_t precision = 24; // FIXME TODO
  const real step_count = count_t(1) << precision;
  const blasC<real> BLA (M, Zp, form, hypot(width, height), pixel_spacing, step_count, &progress[0], running);
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
  const mat2<float> K (1); // FIXME
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
        z = PERTURB(C, Zp[m], c, z);
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
    complex<float> J = complex<float>(float(Zz.dx[0].x), float(Zz.dx[0].y));
    complex<float> dC = Z1 * J * K;
    complex<float> de = norm(Z1) * log(abs(Z1)) / dC;
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

template
  < typename real
  , complex<dual<2, real>> PERTURB(const complex<real> &, const complex<real> &, const complex<dual<2, real>> &, const complex<dual<2, real>> &)
  >
void renderR2(map &out, stats &sta, const param &par, const real Zoom, const count_t M, const complex<real> *Zp, const formulaR2base *form, progress_t *progress, bool *running)
{
#define normx(w) norm(complex<real>((w).x.x, (w).y.x))
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
  const count_t precision = 24; // FIXME TODO
  const real step_count = count_t(1) << precision;
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
  const mat2<float> K (1); // FIXME
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
    real z2 (normx(z));
    complex<dual<2, real>> Zz (Z + z);
    real Zz2 (normx(Zz));

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
        z2 = normx(z);
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
        Zz2 = normx(Zz);
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
        // z = (2 Z + z) z + c
        z = PERTURB(C, Zp[m], c, z);
        z2 = normx(z);
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
        Zz2 = normx(Zz);
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
    mat2<float> J (float(Zz.x.dx[0]), float(Zz.x.dx[1]), float(Zz.y.dx[0]), float(Zz.y.dx[1]));
    complex<float> dC = Z1 * J * K;
    complex<float> de = norm(Z1) * log(abs(Z1)) / dC;
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
#undef normx
}

struct formula
{
  formula() { }
  virtual ~formula() { }
  virtual bool complex_analytic() const noexcept = 0;
  virtual std::string name() const = 0;

  virtual count_t reference(complex<float      > *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const = 0;
  virtual count_t reference(complex<double     > *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const = 0;
  virtual count_t reference(complex<long double> *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const = 0;
  virtual count_t reference(complex<floatexp   > *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const = 0;

  virtual bool period(count_t &period, const complex<float      > *Zp, const count_t M, const complex<float      > c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept = 0;
  virtual bool period(count_t &period, const complex<double     > *Zp, const count_t M, const complex<double     > c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept = 0;
  virtual bool period(count_t &period, const complex<long double> *Zp, const count_t M, const complex<long double> c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept = 0;
  virtual bool period(count_t &period, const complex<floatexp   > *Zp, const count_t M, const complex<floatexp   > c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept = 0;

  virtual bool center(complex<mpreal> &C, const count_t period, progress_t *progress, bool *running) const = 0;

  virtual bool size(floatexp &s, mat2<double> &K, const complex<float      > *Zp, count_t period, const complex<float      > &c0, progress_t *progress, bool *running) const = 0;
  virtual bool size(floatexp &s, mat2<double> &K, const complex<double     > *Zp, count_t period, const complex<double     > &c0, progress_t *progress, bool *running) const = 0;
  virtual bool size(floatexp &s, mat2<double> &K, const complex<long double> *Zp, count_t period, const complex<long double> &c0, progress_t *progress, bool *running) const = 0;
  virtual bool size(floatexp &s, mat2<double> &K, const complex<floatexp   > *Zp, count_t period, const complex<floatexp   > &c0, progress_t *progress, bool *running) const = 0;

  virtual bool domain_size(floatexp &s, const complex<float      > *Zp, count_t period, const complex<float      > &c0, progress_t *progress, bool *running) const = 0;
  virtual bool domain_size(floatexp &s, const complex<double     > *Zp, count_t period, const complex<double     > &c0, progress_t *progress, bool *running) const = 0;
  virtual bool domain_size(floatexp &s, const complex<long double> *Zp, count_t period, const complex<long double> &c0, progress_t *progress, bool *running) const = 0;
  virtual bool domain_size(floatexp &s, const complex<floatexp   > *Zp, count_t period, const complex<floatexp   > &c0, progress_t *progress, bool *running) const = 0;

  virtual void render(map &out, stats &sta, const param &par, const float Zoom, const count_t M, const complex<float> *Zp, progress_t *progress, bool *running) const = 0;
  virtual void render(map &out, stats &sta, const param &par, const double Zoom, const count_t M, const complex<double> *Zp, progress_t *progress, bool *running) const = 0;
  virtual void render(map &out, stats &sta, const param &par, const long double Zoom, const count_t M, const complex<long double> *Zp, progress_t *progress, bool *running) const = 0;
  virtual void render(map &out, stats &sta, const param &par, const floatexp Zoom, const count_t M, const complex<floatexp> *Zp, progress_t *progress, bool *running) const = 0;
};

struct formulaCbase : public formula
{
  formulaCbase() { }
  virtual ~formulaCbase() { }
  virtual bool complex_analytic() const noexcept { return true; }
  virtual void blas_init1(blasC<float> *Bp, const complex<float> *Zp, const float h, const float k, float L, progress_t *progress, bool *running) const noexcept = 0;
  virtual void blas_init1(blasC<double> *Bp, const complex<double> *Zp, const double h, const double k, double L, progress_t *progress, bool *running) const noexcept = 0;
  virtual void blas_init1(blasC<long double> *Bp, const complex<long double> *Zp, const long double h, const long double k, long double L, progress_t *progress, bool *running) const noexcept = 0;
  virtual void blas_init1(blasC<floatexp> *Bp, const complex<floatexp> *Zp, const floatexp h, const floatexp k, floatexp L, progress_t *progress, bool *running) const noexcept = 0;
};

struct formulaR2base : public formula
{
  formulaR2base() { }
  virtual ~formulaR2base() { }
  virtual bool complex_analytic() const noexcept { return false; }
  virtual void blas_init1(blasR2<float> *Bp, const complex<float> *Zp, const float h, const float k, float L, progress_t *progress, bool *running) const noexcept = 0;
  virtual void blas_init1(blasR2<double> *Bp, const complex<double> *Zp, const double h, const double k, double L, progress_t *progress, bool *running) const noexcept = 0;
  virtual void blas_init1(blasR2<long double> *Bp, const complex<long double> *Zp, const long double h, const long double k, long double L, progress_t *progress, bool *running) const noexcept = 0;
  virtual void blas_init1(blasR2<floatexp> *Bp, const complex<floatexp> *Zp, const floatexp h, const floatexp k, floatexp L, progress_t *progress, bool *running) const noexcept = 0;
};

template
  < const char * NAME
  , complex<mpreal> PLAIN_mpreal(const complex<mpreal> &, const complex<mpreal> &)
  , dual<1, complex<mpreal>> PLAIN_dual_mpreal(const dual<1, complex<mpreal>> &, const dual<1, complex<mpreal>> &)
  , dual<1, complex<float>> PERTURB_dual_float(const complex<float>&, const complex<float>&, const dual<1, complex<float>>&, const dual<1, complex<float>>&) noexcept
  , dual<1, complex<double>> PERTURB_dual_double(const complex<double>&, const complex<double>&, const dual<1, complex<double>>&, const dual<1, complex<double>>&) noexcept
  , dual<1, complex<long double>> PERTURB_dual_longdouble(const complex<long double>&, const complex<long double>&, const dual<1, complex<long double>>&, const dual<1, complex<long double>>&) noexcept
  , dual<1, complex<floatexp>> PERTURB_dual_floatexp(const complex<floatexp>&, const complex<floatexp>&, const dual<1, complex<floatexp>>&, const dual<1, complex<floatexp>>&) noexcept
  , blaC<float> BLA_float(const float &h, const float &k, const float &L, const complex<float> &Z) noexcept
  , blaC<double> BLA_double(const double &h, const double &k, const double &L, const complex<double> &Z) noexcept
  , blaC<long double> BLA_longdouble(const long double &h, const long double &k, const long double &L, const complex<long double> &Z) noexcept
  , blaC<floatexp> BLA_floatexp(const floatexp &h, const floatexp &k, const floatexp &L, const complex<floatexp> &Z) noexcept
  >
struct formulaC : public formulaCbase
{
  formulaC() { }
  virtual ~formulaC() { }
  virtual std::string name() const { return NAME; }

  virtual count_t reference(complex<float      > *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const
  {
    return ::reference<float, PLAIN_mpreal>(Zp, MaxRefIters, C, progress, running);
  }
  virtual count_t reference(complex<double     > *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const
  {
    return ::reference<double, PLAIN_mpreal>(Zp, MaxRefIters, C, progress, running);
  }
  virtual count_t reference(complex<long double> *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const
  {
    return ::reference<long double, PLAIN_mpreal>(Zp, MaxRefIters, C, progress, running);
  }
  virtual count_t reference(complex<floatexp   > *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const
  {
    return ::reference<floatexp, PLAIN_mpreal>(Zp, MaxRefIters, C, progress, running);
  }

  virtual bool period(count_t &period, const complex<float      > *Zp, const count_t M, const complex<float      > c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept
  {
    return periodC<float, PERTURB_dual_float>(period, Zp, M, c, N, s, K, progress, running);
  }
  virtual bool period(count_t &period, const complex<double     > *Zp, const count_t M, const complex<double     > c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept
  {
    return periodC<double, PERTURB_dual_double>(period, Zp, M, c, N, s, K, progress, running);
  }
  virtual bool period(count_t &period, const complex<long double> *Zp, const count_t M, const complex<long double> c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept
  {
    return periodC<long double, PERTURB_dual_longdouble>(period, Zp, M, c, N, s, K, progress, running);
  }
  virtual bool period(count_t &period, const complex<floatexp   > *Zp, const count_t M, const complex<floatexp   > c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept
  {
    return periodC<floatexp, PERTURB_dual_floatexp>(period, Zp, M, c, N, s, K, progress, running);
  }

  virtual bool center(complex<mpreal> &C, const count_t period, progress_t *progress, bool *running) const
  {
    return centerC<PLAIN_dual_mpreal>(C, period, progress, running);
  }

  virtual bool size(floatexp &s, mat2<double> &K, const complex<float      > *Zp, count_t period, const complex<float      > &c0, progress_t *progress, bool *running) const noexcept
  {
    return sizeC<float, PERTURB_dual_float>(s, K, Zp, period, c0, progress, running);
  }
  virtual bool size(floatexp &s, mat2<double> &K, const complex<double     > *Zp, count_t period, const complex<double     > &c0, progress_t *progress, bool *running) const noexcept
  {
    return sizeC<double, PERTURB_dual_double>(s, K, Zp, period, c0, progress, running);
  }
  virtual bool size(floatexp &s, mat2<double> &K, const complex<long double> *Zp, count_t period, const complex<long double> &c0, progress_t *progress, bool *running) const noexcept
  {
    return sizeC<long double, PERTURB_dual_longdouble>(s, K, Zp, period, c0, progress, running);
  }
  virtual bool size(floatexp &s, mat2<double> &K, const complex<floatexp   > *Zp, count_t period, const complex<floatexp   > &c0, progress_t *progress, bool *running) const noexcept
  {
    return sizeC<floatexp, PERTURB_dual_floatexp>(s, K, Zp, period, c0, progress, running);
  }

  virtual bool domain_size(floatexp &s, const complex<float      > *Zp, count_t period, const complex<float      > &c0, progress_t *progress, bool *running) const noexcept
  {
    return domain_sizeC<float, PERTURB_dual_float>(s, Zp, period, c0, progress, running);
  }
  virtual bool domain_size(floatexp &s, const complex<double     > *Zp, count_t period, const complex<double     > &c0, progress_t *progress, bool *running) const noexcept
  {
    return domain_sizeC<double, PERTURB_dual_double>(s, Zp, period, c0, progress, running);
  }
  virtual bool domain_size(floatexp &s, const complex<long double> *Zp, count_t period, const complex<long double> &c0, progress_t *progress, bool *running) const noexcept
  {
    return domain_sizeC<long double, PERTURB_dual_longdouble>(s, Zp, period, c0, progress, running);
  }
  virtual bool domain_size(floatexp &s, const complex<floatexp   > *Zp, count_t period, const complex<floatexp   > &c0, progress_t *progress, bool *running) const noexcept
  {
    return domain_sizeC<floatexp, PERTURB_dual_floatexp>(s, Zp, period, c0, progress, running);
  }

  virtual void blas_init1(blasC<float> *Bp, const complex<float> *Zp, const float h, const float k, float L, progress_t *progress, bool *running) const noexcept
  {
    return blas_init1C<float, BLA_float>(Bp, Zp, h, k, L, progress, running);
  }
  virtual void blas_init1(blasC<double> *Bp, const complex<double> *Zp, const double h, const double k, double L, progress_t *progress, bool *running) const noexcept
  {
    return blas_init1C<double, BLA_double>(Bp, Zp, h, k, L, progress, running);
  }
  virtual void blas_init1(blasC<long double> *Bp, const complex<long double> *Zp, const long double h, const long double k, long double L, progress_t *progress, bool *running) const noexcept
  {
    return blas_init1C<long double, BLA_longdouble>(Bp, Zp, h, k, L, progress, running);
  }
  virtual void blas_init1(blasC<floatexp> *Bp, const complex<floatexp> *Zp, const floatexp h, const floatexp k, floatexp L, progress_t *progress, bool *running) const noexcept
  {
    return blas_init1C<floatexp, BLA_floatexp>(Bp, Zp, h, k, L, progress, running);
  }

  virtual void render(map &out, stats &sta, const param &par, const float Zoom, const count_t M, const complex<float> *Zp, progress_t *progress, bool *running) const
  {
    return renderC<float, PERTURB_dual_float>(out, sta, par, Zoom, M, Zp, dynamic_cast<const formulaCbase *>(this), progress, running);
  }
  virtual void render(map &out, stats &sta, const param &par, const double Zoom, const count_t M, const complex<double> *Zp, progress_t *progress, bool *running) const
  {
    return renderC<double, PERTURB_dual_double>(out, sta, par, Zoom, M, Zp, dynamic_cast<const formulaCbase *>(this), progress, running);
  }
  virtual void render(map &out, stats &sta, const param &par, const long double Zoom, const count_t M, const complex<long double> *Zp, progress_t *progress, bool *running) const
  {
    return renderC<long double, PERTURB_dual_longdouble>(out, sta, par, Zoom, M, Zp, dynamic_cast<const formulaCbase *>(this), progress, running);
  }
  virtual void render(map &out, stats &sta, const param &par, const floatexp Zoom, const count_t M, const complex<floatexp> *Zp, progress_t *progress, bool *running) const
  {
    return renderC<floatexp, PERTURB_dual_floatexp>(out, sta, par, Zoom, M, Zp, dynamic_cast<const formulaCbase *>(this), progress, running);
  }
};

template
  < const char * NAME
  , complex<mpreal> PLAIN_mpreal(const complex<mpreal> &, const complex<mpreal> &)
  , complex<dual<2, mpreal>> PLAIN_dual_mpreal(const complex<dual<2, mpreal>> &, const complex<dual<2, mpreal>> &)
  , complex<dual<2, float>> PERTURB_dual_float(const complex<float>&, const complex<float>&, const complex<dual<2, float>>&, const complex<dual<2, float>>&) noexcept
  , complex<dual<2, double>> PERTURB_dual_double(const complex<double>&, const complex<double>&, const complex<dual<2, double>>&, const complex<dual<2, double>>&) noexcept
  , complex<dual<2, long double>> PERTURB_dual_longdouble(const complex<long double>&, const complex<long double>&, const complex<dual<2, long double>>&, const complex<dual<2, long double>>&) noexcept
  , complex<dual<2, floatexp>> PERTURB_dual_floatexp(const complex<floatexp>&, const complex<floatexp>&, const complex<dual<2, floatexp>>&, const complex<dual<2, floatexp>>&) noexcept
  , blaR2<float> BLA_float(const float &h, const float &k, const float &L, const complex<float> &Z) noexcept
  , blaR2<double> BLA_double(const double &h, const double &k, const double &L, const complex<double> &Z) noexcept
  , blaR2<long double> BLA_longdouble(const long double &h, const long double &k, const long double &L, const complex<long double> &Z) noexcept
  , blaR2<floatexp> BLA_floatexp(const floatexp &h, const floatexp &k, const floatexp &L, const complex<floatexp> &Z) noexcept
  >
struct formulaR2 : public formulaR2base
{
  formulaR2() { }
  virtual ~formulaR2() { }
  virtual std::string name() const { return NAME; }
  virtual bool complex_analytic() const noexcept { return false; }

  virtual count_t reference(complex<float      > *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const
  {
    return ::reference<float, PLAIN_mpreal>(Zp, MaxRefIters, C, progress, running);
  }
  virtual count_t reference(complex<double     > *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const
  {
    return ::reference<double, PLAIN_mpreal>(Zp, MaxRefIters, C, progress, running);
  }
  virtual count_t reference(complex<long double> *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const
  {
    return ::reference<long double, PLAIN_mpreal>(Zp, MaxRefIters, C, progress, running);
  }
  virtual count_t reference(complex<floatexp   > *Zp, const count_t MaxRefIters, const complex<mpreal> &C, progress_t *progress, bool *running) const
  {
    return ::reference<floatexp, PLAIN_mpreal>(Zp, MaxRefIters, C, progress, running);
  }

  virtual bool period(count_t &period, const complex<float      > *Zp, const count_t M, const complex<float      > c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept
  {
    return periodR2<float, PERTURB_dual_float>(period, Zp, M, c, N, s, K, progress, running);
  }
  virtual bool period(count_t &period, const complex<double     > *Zp, const count_t M, const complex<double     > c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept
  {
    return periodR2<double, PERTURB_dual_double>(period, Zp, M, c, N, s, K, progress, running);
  }
  virtual bool period(count_t &period, const complex<long double> *Zp, const count_t M, const complex<long double> c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept
  {
    return periodR2<long double, PERTURB_dual_longdouble>(period, Zp, M, c, N, s, K, progress, running);
  }
  virtual bool period(count_t &period, const complex<floatexp   > *Zp, const count_t M, const complex<floatexp   > c, const count_t N, const floatexp &s, const mat2<double> &K, progress_t *progress, bool *running) const noexcept
  {
    return periodR2<floatexp, PERTURB_dual_floatexp>(period, Zp, M, c, N, s, K, progress, running);
  }

  virtual bool center(complex<mpreal> &C, const count_t period, progress_t *progress, bool *running) const
  {
    return centerR2<PLAIN_dual_mpreal>(C, period, progress, running);
  }

  virtual bool size(floatexp &s, mat2<double> &K, const complex<float      > *Zp, count_t period, const complex<float      > &c0, progress_t *progress, bool *running) const noexcept
  {
    return sizeR2<float, PERTURB_dual_float>(s, K, Zp, period, c0, progress, running);
  }
  virtual bool size(floatexp &s, mat2<double> &K, const complex<double     > *Zp, count_t period, const complex<double     > &c0, progress_t *progress, bool *running) const noexcept
  {
    return sizeR2<double, PERTURB_dual_double>(s, K, Zp, period, c0, progress, running);
  }
  virtual bool size(floatexp &s, mat2<double> &K, const complex<long double> *Zp, count_t period, const complex<long double> &c0, progress_t *progress, bool *running) const noexcept
  {
    return sizeR2<long double, PERTURB_dual_longdouble>(s, K, Zp, period, c0, progress, running);
  }
  virtual bool size(floatexp &s, mat2<double> &K, const complex<floatexp   > *Zp, count_t period, const complex<floatexp   > &c0, progress_t *progress, bool *running) const noexcept
  {
    return sizeR2<floatexp, PERTURB_dual_floatexp>(s, K, Zp, period, c0, progress, running);
  }

  virtual bool domain_size(floatexp &s, const complex<float      > *Zp, count_t period, const complex<float      > &c0, progress_t *progress, bool *running) const noexcept
  {
    return domain_sizeR2<float, PERTURB_dual_float>(s, Zp, period, c0, progress, running);
  }
  virtual bool domain_size(floatexp &s, const complex<double     > *Zp, count_t period, const complex<double     > &c0, progress_t *progress, bool *running) const noexcept
  {
    return domain_sizeR2<double, PERTURB_dual_double>(s, Zp, period, c0, progress, running);
  }
  virtual bool domain_size(floatexp &s, const complex<long double> *Zp, count_t period, const complex<long double> &c0, progress_t *progress, bool *running) const noexcept
  {
    return domain_sizeR2<long double, PERTURB_dual_longdouble>(s, Zp, period, c0, progress, running);
  }
  virtual bool domain_size(floatexp &s, const complex<floatexp   > *Zp, count_t period, const complex<floatexp   > &c0, progress_t *progress, bool *running) const noexcept
  {
    return domain_sizeR2<floatexp, PERTURB_dual_floatexp>(s, Zp, period, c0, progress, running);
  }

  virtual void blas_init1(blasR2<float> *Bp, const complex<float> *Zp, const float h, const float k, float L, progress_t *progress, bool *running) const noexcept
  {
    return blas_init1R2<float, BLA_float>(Bp, Zp, h, k, L, progress, running);
  }
  virtual void blas_init1(blasR2<double> *Bp, const complex<double> *Zp, const double h, const double k, double L, progress_t *progress, bool *running) const noexcept
  {
    return blas_init1R2<double, BLA_double>(Bp, Zp, h, k, L, progress, running);
  }
  virtual void blas_init1(blasR2<long double> *Bp, const complex<long double> *Zp, const long double h, const long double k, long double L, progress_t *progress, bool *running) const noexcept
  {
    return blas_init1R2<long double, BLA_longdouble>(Bp, Zp, h, k, L, progress, running);
  }
  virtual void blas_init1(blasR2<floatexp> *Bp, const complex<floatexp> *Zp, const floatexp h, const floatexp k, floatexp L, progress_t *progress, bool *running) const noexcept
  {
    return blas_init1R2<floatexp, BLA_floatexp>(Bp, Zp, h, k, L, progress, running);
  }

  virtual void render(map &out, stats &sta, const param &par, const float Zoom, const count_t M, const complex<float> *Zp, progress_t *progress, bool *running) const
  {
    return renderR2<float, PERTURB_dual_float>(out, sta, par, Zoom, M, Zp, dynamic_cast<const formulaR2base *>(this), progress, running);
  }
  virtual void render(map &out, stats &sta, const param &par, const double Zoom, const count_t M, const complex<double> *Zp, progress_t *progress, bool *running) const
  {
    return renderR2<double, PERTURB_dual_double>(out, sta, par, Zoom, M, Zp, dynamic_cast<const formulaR2base *>(this), progress, running);
  }
  virtual void render(map &out, stats &sta, const param &par, const long double Zoom, const count_t M, const complex<long double> *Zp, progress_t *progress, bool *running) const
  {
    return renderR2<long double, PERTURB_dual_longdouble>(out, sta, par, Zoom, M, Zp, dynamic_cast<const formulaR2base *>(this), progress, running);
  }
  virtual void render(map &out, stats &sta, const param &par, const floatexp Zoom, const count_t M, const complex<floatexp> *Zp, progress_t *progress, bool *running) const
  {
    return renderR2<floatexp, PERTURB_dual_floatexp>(out, sta, par, Zoom, M, Zp, dynamic_cast<const formulaR2base *>(this), progress, running);
  }
};

extern std::vector<formula *> formulas;
void formulas_init();
