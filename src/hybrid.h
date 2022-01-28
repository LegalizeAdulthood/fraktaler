// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <thread>
#include <vector>

#include "map.h"
#include "param.h"
#include "stats.h"

// http://www.burtleburtle.net/bob/hash/integer.html
inline CONSTEXPR uint32_t burtle_hash(uint32_t a) noexcept
{
  a = (a+0x7ed55d16) + (a<<12);
  a = (a^0xc761c23c) ^ (a>>19);
  a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);
  a = (a+0xfd7046c5) + (a<<3);
  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}

inline double radical_inverse(coord_t a, const coord_t base) noexcept
{
  constexpr double one_minus_epsilon = 0.99999999999999989;
  const double base1 = 1.0 / base;
  coord_t reversed = 0;
  double base1n = 1;
  while (a)
  {
    const coord_t next  = a / base;
    const coord_t digit = a - base * next;
    reversed = reversed * base + digit;
    base1n *= base1;
    a = next;
  }
  return std::min(reversed * base1n, one_minus_epsilon);
}

inline double wrap(const double v) noexcept
{
  return v - std::floor(v);
}

inline double triangle(const double a) noexcept
{
  const double b = a * 2 - 1;
  const double c = std::sqrt(std::abs(b));
  const double e = b > 0 ? c - 1 : 1 - c;
  return e;
}

inline void jitter(const coord_t width, const coord_t height, const coord_t i, const coord_t j, const coord_t k, double &x, double &y) noexcept
{
  coord_t ix = (k * height + j) * width + i;
  double h = burtle_hash(ix) / double(0x100000000LL);
  x = triangle(wrap(radical_inverse(k, 2) + h));
  y = triangle(wrap(radical_inverse(k, 3) + h));
}

template <typename T>
inline complex<T> hybrid_plain(const struct phybrid1 &H, const complex<T> &C, const complex<T> &Z)
{
  complex<T> W = Z;
  if (H.abs_x) W.x = abs(W.x);
  if (H.abs_y) W.y = abs(W.y);
  if (H.neg_x) W.x = -W.x;
  if (H.neg_y) W.y = -W.y;
  return pow(W, H.power) + C;
}

template <typename T, typename t>
inline constexpr complex<t> hybrid_perturb(const struct phybrid1 &H, const complex<T> &C, const complex<T> &Z, const complex<t> &c, const complex<t> &z) noexcept
{
  T X = Z.x;
  T Y = Z.y;
  t x = z.x;
  t y = z.y;
  complex<t> W = Z + z;
  complex<T> B = Z;
  if (H.abs_x)
  {
    x = diffabs(X, x);
    W.x = abs(W.x);
    B.x = abs(B.x);
  }
  if (H.abs_y)
  {
    y = diffabs(Y, y);
    W.y = abs(W.y);
    B.y = abs(B.y);
  }
  if (H.neg_x)
  {
    x = -x;
    W.x = -W.x;
    B.x = -B.x;
  }
  if (H.neg_y)
  {
    y = -y;
    W.y = -W.y;
    B.y = -B.y;
  }
  complex<t> P(x, y);
  complex<t> S(0);
  for (int i = 0; i <= H.power - 1; ++i)
  {
    int j = H.power - 1 - i;
    S += pow(W, i) * pow(B, j);
  }
  return P * S + c;
}

template <typename real>
inline constexpr blaR2<real> hybrid_bla(const struct phybrid1 &H, const real &h, const real &k, const real &L, const complex<real> &Z) noexcept
{
  using std::abs, ::abs;
  using std::min;
  using std::max;
  dual<2, real> x(Z.x); x.dx[0] = 1;
  dual<2, real> y(Z.y); y.dx[1] = 1;
  complex<dual<2, real>> W(x, y);
  complex<dual<2, real>> C(0, 0);
  W = hybrid_plain(H, C, W);
  const mat2<real> A(W.x.dx[0], W.x.dx[1], W.y.dx[0], W.y.dx[1]);
  const mat2<real> B(1);
  const real mZ = min(abs(Z.x), abs(Z.y)) / 2; // FIXME arbitrary factor
  const real mA = abs(A);
  const real mB = abs(B);
  const real r = max(real(0), (mZ - mB * h * k) / (mA + 1)) / L;
  const real r2 = r * r;
  const count_t l = 1;
  blaR2<real> b = { A, B, r2, l };
  return b;
}

template <typename t>
void hybrid_blas(std::vector<blasR2<t>> &B, const std::vector<std::vector<complex<t>>> &Z, const phybrid &H, t h, t k, t L, volatile progress_t *progress, volatile bool *running)
{
  count_t count = H.per.size();
  for (count_t phase = 0; phase < count; ++phase)
  {
    B.push_back(blasR2(Z[phase], H, phase, h, k, L, &progress[1], running));
    progress[0] = progress_t(phase + 1) / progress_t(count);
  }
}


template <typename t>
count_t hybrid_reference(complex<t> *Zp, const struct phybrid &H, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running)
{
  complex<mpreal> Z (0);
  count_t M = MaxRefIters;
  // calculate reference in high precision
  for (count_t i = 0; i < MaxRefIters; ++i)
  {
    // store low precision orbit
    Zp[i] = complex<t>(convert<t>(Z.x), convert<t>(Z.y));
    // escape check
    if (norm(Zp[i]) > 4 || ! *running) // FIXME escape radius
    {
      M = i;
      break;
    }
    // step
    Z = hybrid_plain(H.per[(phase + i) % H.per.size()], C, Z);
    *progress = (i + 1) / progress_t(MaxRefIters);
  }
  return M;
}

template <typename t>
void hybrid_references(std::vector<std::vector<complex<t>>> &Zp, const struct phybrid &H, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running)
{
  parallel1d(std::thread::hardware_concurrency(), 0, H.per.size(), 1, running, [&](count_t phase)
  {
    count_t M = hybrid_reference(&Zp[phase][0], H, phase, MaxRefIters, C, &progress[phase], running);
    Zp[phase].resize(M);
  });
}

template <typename real, bool gather_statistics>
void hybrid_render_stats(map &out, stats &sta, const phybrid &H, const std::vector<blasR2<real>> &bla, const count_t subframe, const param &par, const real Zoom, const complex<real> offset, const std::vector<std::vector<complex<real>>> &Zp, volatile progress_t *progress, volatile bool *running)
{
#define normx(w) norm(complex<real>((w).x.x, (w).y.x))
  using std::isinf;
  using std::isnan;
  using std::log;
  using std::max;
  using std::min;
  const coord_t width = out.width;
  const coord_t height = out.height;
  const count_t Iterations = par.p.bailout.iterations;
  const count_t ReferencePeriod = par.p.reference.period;
  const count_t PerturbIterations = par.p.bailout.maximum_perturb_iterations;
  const real ER2 = par.p.bailout.escape_radius * par.p.bailout.escape_radius;
  const real pixel_spacing = 4 / Zoom / height;
  const mat2<real> K (real(par.transform.x[0][0]), real(par.transform.x[0][1]), real(par.transform.x[1][0]), real(par.transform.x[1][1]));
  const mat2<float> Kf (float(par.transform.x[0][0]), float(par.transform.x[0][1]), float(par.transform.x[1][0]), float(par.transform.x[1][1]));
  const float degree (2); // FIXME
  std::atomic<count_t> pixels = 0;
  sta += parallel2dr<stats>(std::thread::hardware_concurrency(), 0, width, 32, 0, height, 32, running, [&](coord_t i, coord_t j) -> stats
  {
    // statistics
    count_t iters_ptb = 0;
    count_t iters_bla = 0;
    count_t steps_ptb = 0;
    count_t steps_bla = 0;
    count_t rebases_small = 0;
    count_t rebases_noref = 0;
    count_t iters_ref = 0;
    double di, dj;
    jitter(width, height, i, j, subframe, di, dj);
	  dual<2, real> u0(real(i + di)); u0.dx[0] = real(1);
	  dual<2, real> v0(real(j + dj)); v0.dx[1] = real(1);
    if (par.p.transform.exponential_map)
    {
      auto re = (-0.6931471805599453 / height) * v0; // log 2
      auto im = (6.283185307179586 / width) * u0; // 2 pi
      auto R = 0.5 * std::hypot(width, height);
      auto r = exp(re);
      auto c = cos(im);
      auto s = sin(im);
      u0 = R * r * c;
      v0 = R * r * s;
    }
    else
    {
      u0 -= real(width / 2.0);
      v0 -= real(height / 2.0);
    }
    // FIXME should K multiply offset?
    dual<2, real> cx (u0 * pixel_spacing + offset.x);
    cx.dx[0] = pixel_spacing;
    dual<2, real> cy (v0 * pixel_spacing + offset.y);
    cy.dx[1] = pixel_spacing;
    const complex<real> C (Zp[0][1]); // FIXME
    if constexpr (gather_statistics)
    {
      iters_ref = 1;
    }
    complex<dual<2, real>> c (cx, cy);
    c = K * c;
    count_t phase = 0;
    count_t m = 0;
    count_t n = 0;
    complex<real> Z (Zp[phase][0]);
    complex<dual<2, real>> z (0);
    real z2 (normx(z));
    complex<dual<2, real>> Zz (Z + z);
    real Zz2 (normx(Zz));

    while (n < Iterations && Zz2 < ER2 && iters_ptb < PerturbIterations)
    {
      // bla steps
      const blaR2<real> *b = 0;
      while (n < Iterations && Zz2 < ER2 && (b = bla[phase].lookup(m, z2)))
      {
        const mat2<real> A = b->A;
        const mat2<real> B = b->B;
        count_t l = b->l;
        z = A * z + B * c;
        z2 = normx(z);
        n += l;
        m += l;
        if constexpr (gather_statistics)
        {
          steps_bla++;
          iters_bla += l;
        }
        if (ReferencePeriod)
        {
          while (m >= ReferencePeriod)
          {
            phase = (phase + m - ReferencePeriod) % Zp.size();
            m -= ReferencePeriod;
          }
        }

        // rebase
        if (! (n < Iterations && Zz2 < ER2 && iters_ptb < PerturbIterations))
        {
          break;
        }
        if (! (m < Zp[phase].size()))
        {
          break;
        }
        complex<real> Z = Zp[phase][m];
        if constexpr (gather_statistics)
        {
          iters_ref = iters_ref > m ? iters_ref : m;
        }
        Zz = Z + z;
        Zz2 = normx(Zz);
        if (Zz2 < z2 || (ReferencePeriod == 0 && m + 1 == Zp[phase].size()))
        {
          z = Zz;
          phase = (phase + m) % Zp.size();
          m = 0;
          if constexpr (gather_statistics)
          {
            if (Zz2 < z2)
            {
              rebases_small++;
            }
            else
            {
              rebases_noref++;
            }
          }
        }
      }

      // perturbation iteration
      {
        if (! (n < Iterations && Zz2 < ER2 && iters_ptb < PerturbIterations))
        {
          break;
        }
        if (! (m < Zp[phase].size()))
        {
          break;
        }
        // z = (2 Z + z) z + c
        z = hybrid_perturb(H.per[phase], C, Zp[phase][m], c, z);
        iters_ref = iters_ref > m ? iters_ref : m;
        z2 = normx(z);
        n++;
        m++;
        if constexpr (gather_statistics)
        {
          steps_ptb++;
        }
        iters_ptb++;
        if (ReferencePeriod)
        {
          while (m >= ReferencePeriod)
          {
            phase = (phase + m - ReferencePeriod) % Zp.size();
            m -= ReferencePeriod;
          }
        }
      }

      {
        // rebase
        if (! (n < Iterations && Zz2 < ER2 && iters_ptb < PerturbIterations))
        {
          break;
        }
        if (! (m < Zp[phase].size()))
        {
          break;
        }
        complex<real> Z = Zp[phase][m];
        if constexpr (gather_statistics)
        {
          iters_ref = iters_ref > m ? iters_ref : m;
        }
        Zz = Z + z;
        Zz2 = normx(Zz);
        if (Zz2 < z2 || (ReferencePeriod == 0 && m + 1 == Zp[phase].size()))
        {
          z = Zz;
          phase = (phase + m) % Zp.size();
          m = 0;
          if constexpr (gather_statistics)
          {
            if (Zz2 < z2)
            {
              rebases_small++;
            }
            else
            {
              rebases_noref++;
            }
          }
        }
      }
    }

    // compute output
    complex<float> Z1 = complex<float>(float(Zz.x.x), float(Zz.y.x));
    mat2<float> J (float(Zz.x.dx[0]), float(Zz.x.dx[1]), float(Zz.y.dx[0]), float(Zz.y.dx[1]));
    complex<float> dC = Z1 * J * Kf;
    complex<float> de = norm(Z1) * log(abs(Z1)) / dC;
    float nf = std::min(std::max(1 - log(log(norm(Z1)) / log(float(ER2))) / log(degree), 0.f), 1.f);
    float t = arg(Z1) / (2.0f * 3.141592653f);
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
    const count_t count = pixels.fetch_add(1);
    progress[0] = count / progress_t(width * height);
    return stats(iters_ptb + iters_bla, iters_ptb, iters_bla, steps_ptb + steps_bla, steps_ptb, steps_bla, rebases_small + rebases_noref, rebases_small, rebases_noref, iters_ref);
  });
#undef normx
}

template <typename real>
void hybrid_render(map &out, stats &sta, const phybrid &H, const std::vector<blasR2<real>> &bla, const count_t subframe, const param &par, const real Zoom, const complex<real> offset, const std::vector<std::vector<complex<real>>> &Zp, volatile progress_t *progress, volatile bool *running)
{
  if (subframe == 0)
  {
    hybrid_render_stats<real, true>(out, sta, H, bla, subframe, par, Zoom, offset, Zp, progress, running);
  }
  else
  {
    stats dummy;
    hybrid_render_stats<real, false>(out, dummy, H, bla, subframe, par, Zoom, offset, Zp, progress, running);
  }
}
