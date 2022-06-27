// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <cmath>
#include <thread>

#include "bla.h"
#include "hybrid.h"
#include "map.h"
#include "parallel.h"
#include "softfloat.h"
#include "stats.h"

template <typename t>
void hybrid_blas(std::vector<blasR2<t>> &B, const std::vector<std::vector<complex<t>>> &Z, const phybrid &H, t h, t k, t L, volatile progress_t *progress, volatile bool *running)
{
  count_t count = H.per.size();
  for (count_t phase = 0; phase < count; ++phase)
  {
    B.push_back(blasR2(Z[phase], H, phase, h, k, L, &progress[phase], running));
  }
}

template void hybrid_blas(std::vector<blasR2<float>> &B, const std::vector<std::vector<complex<float>>> &Z, const phybrid &H, float h, float k, float L, volatile progress_t *progress, volatile bool *running);
template void hybrid_blas(std::vector<blasR2<double>> &B, const std::vector<std::vector<complex<double>>> &Z, const phybrid &H, double h, double k, double L, volatile progress_t *progress, volatile bool *running);
template void hybrid_blas(std::vector<blasR2<long double>> &B, const std::vector<std::vector<complex<long double>>> &Z, const phybrid &H, long double h, long double k, long double L, volatile progress_t *progress, volatile bool *running);
template void hybrid_blas(std::vector<blasR2<floatexp>> &B, const std::vector<std::vector<complex<floatexp>>> &Z, const phybrid &H, floatexp h, floatexp k, floatexp L, volatile progress_t *progress, volatile bool *running);
template void hybrid_blas(std::vector<blasR2<softfloat>> &B, const std::vector<std::vector<complex<softfloat>>> &Z, const phybrid &H, softfloat h, softfloat k, softfloat L, volatile progress_t *progress, volatile bool *running);
#ifdef HAVE_FLOAT128
template void hybrid_blas(std::vector<blasR2<float128>> &B, const std::vector<std::vector<complex<float128>>> &Z, const phybrid &H, float128 h, float128 k, float128 L, volatile progress_t *progress, volatile bool *running);
#endif

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

template count_t hybrid_reference(complex<float> *Zp, const struct phybrid &H, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template count_t hybrid_reference(complex<double> *Zp, const struct phybrid &H, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template count_t hybrid_reference(complex<long double> *Zp, const struct phybrid &H, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template count_t hybrid_reference(complex<floatexp> *Zp, const struct phybrid &H, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template count_t hybrid_reference(complex<softfloat> *Zp, const struct phybrid &H, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
#ifdef HAVE_FLOAT128
template count_t hybrid_reference(complex<float128> *Zp, const struct phybrid &H, const count_t &phase, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
#endif

template <typename t>
void hybrid_references(std::vector<std::vector<complex<t>>> &Zp, const struct phybrid &H, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running)
{
  parallel1d(std::thread::hardware_concurrency(), 0, H.per.size(), 1, running, [&](count_t phase)
  {
    count_t M = hybrid_reference(&Zp[phase][0], H, phase, MaxRefIters, C, &progress[phase], running);
    Zp[phase].resize(M);
  });
}

template void hybrid_references(std::vector<std::vector<complex<float>>> &Zp, const struct phybrid &H, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template void hybrid_references(std::vector<std::vector<complex<double>>> &Zp, const struct phybrid &H, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template void hybrid_references(std::vector<std::vector<complex<long double>>> &Zp, const struct phybrid &H, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template void hybrid_references(std::vector<std::vector<complex<floatexp>>> &Zp, const struct phybrid &H, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
template void hybrid_references(std::vector<std::vector<complex<softfloat>>> &Zp, const struct phybrid &H, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
#ifdef HAVE_FLOAT128
template void hybrid_references(std::vector<std::vector<complex<float128>>> &Zp, const struct phybrid &H, const count_t &MaxRefIters, const complex<mpreal> &C, volatile progress_t *progress, volatile bool *running);
#endif

template <typename real, bool gather_statistics>
void hybrid_render_stats(coord_t frame, map &out, stats &sta, const phybrid &H, const std::vector<blasR2<real>> &bla, const count_t subframe, const param &par, const real Zoom, const complex<real> offset, const std::vector<std::vector<complex<real>>> &Zp, volatile progress_t *progress, volatile bool *running)
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
  const real IR = par.p.bailout.inscape_radius;
  const real pixel_spacing = 4 / Zoom / height;
  const mat2<real> K (real(par.transform.x[0][0]), real(par.transform.x[0][1]), real(par.transform.x[1][0]), real(par.transform.x[1][1]));
  const mat2<float> Kf (float(par.transform.x[0][0]), float(par.transform.x[0][1]), float(par.transform.x[1][0]), float(par.transform.x[1][1]));
  const float degree (2); // FIXME
  std::atomic<count_t> pixels = 0;
  sta += parallel2dr<stats>(std::thread::hardware_concurrency(), 0, width, 32, 0, height, 32, running, [&](coord_t i, coord_t j) -> stats
  {
    // statistics
    count_t iters_ptb = 1;
    count_t iters_bla = 0;
    count_t steps_ptb = 1;
    count_t steps_bla = 0;
    count_t rebases_small = 0;
    count_t rebases_noref = 0;
    count_t iters_ref = 2;
    double di, dj;
    jitter(width, height, frame, i, j, subframe, di, dj);
	  dual<4, real> u0(real(i + di)); u0.dx[0] = real(1);
	  dual<4, real> v0(real(j + dj)); v0.dx[1] = real(1);
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
    dual<4, real> cx (u0 * pixel_spacing + offset.x);
    dual<4, real> cy (v0 * pixel_spacing + offset.y);
    const complex<real> C (Zp[0][1]); // FIXME
    if constexpr (gather_statistics)
    {
      iters_ref = 2;
    }
    complex<dual<4, real>> c (cx, cy);
    c = K * c;
    count_t phase = 0;
    count_t m = 1;
    count_t n = 1;
    complex<real> Z (Zp[phase][m]);
    complex<dual<4, real>> z (c);
    z.x.dx[2] = real(1);
    z.y.dx[3] = real(1);
    real z2 (normx(z));
    complex<dual<4, real>> Zz (Z + z);
    real Zz2 (normx(Zz));
    real dZ (sup(mat2<real>(Zz.x.dx[2], Zz.x.dx[3], Zz.y.dx[2], Zz.y.dx[3])));
    while
      ( n < Iterations &&
        Zz2 < ER2 &&
        IR  < dZ &&
        iters_ptb < PerturbIterations
      )
    {
      // bla steps
      const blaR2<real> *b = 0;
      while
        ( n < Iterations &&
          Zz2 < ER2 &&
          IR  < dZ &&
          (b = bla[phase].lookup(m, z2))
        )
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
        if (!
          ( n < Iterations &&
            Zz2 < ER2 &&
            IR  < dZ &&
            iters_ptb < PerturbIterations)
          )
        {
          break;
        }
        if (! (m < count_t(Zp[phase].size())))
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
        dZ = sup(mat2<real>(Zz.x.dx[2], Zz.x.dx[3], Zz.y.dx[2], Zz.y.dx[3]));
        if (Zz2 < z2 || (ReferencePeriod == 0 && m + 1 == count_t(Zp[phase].size())))
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
        if (! (m < count_t(Zp[phase].size())))
        {
          break;
        }
        // z = f(C, Z, c, z)
        z = hybrid_perturb(H.per[n % H.per.size()], C, Zp[phase][m], c, z);
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
        if (!
          ( n < Iterations &&
            Zz2 < ER2 &&
            IR < dZ &&
            iters_ptb < PerturbIterations)
          )
        {
          break;
        }
        if (! (m < count_t(Zp[phase].size())))
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
        dZ = sup(mat2<real>(Zz.x.dx[2], Zz.x.dx[3], Zz.y.dx[2], Zz.y.dx[3]));
        if (Zz2 < z2 || (ReferencePeriod == 0 && m + 1 == count_t(Zp[phase].size())))
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
    return stats(iters_ptb + iters_bla, iters_ptb, iters_bla, steps_ptb + steps_bla, steps_ptb, steps_bla, rebases_small + rebases_noref, rebases_small, rebases_noref, iters_ref, ! (Zz2 < ER2), ! (IR < dZ));
  });
#undef normx
}

template <typename real>
void hybrid_render(coord_t frame, map &out, stats &sta, const phybrid &H, const std::vector<blasR2<real>> &bla, const count_t subframe, const param &par, const real Zoom, const complex<real> offset, const std::vector<std::vector<complex<real>>> &Zp, volatile progress_t *progress, volatile bool *running)
{
  if (subframe == 0)
  {
    hybrid_render_stats<real, true>(frame, out, sta, H, bla, subframe, par, Zoom, offset, Zp, progress, running);
  }
  else
  {
    stats dummy;
    hybrid_render_stats<real, false>(frame, out, dummy, H, bla, subframe, par, Zoom, offset, Zp, progress, running);
  }
}

template void hybrid_render(coord_t frame, map &out, stats &sta, const phybrid &H, const std::vector<blasR2<float>> &bla, const count_t subframe, const param &par, const float Zoom, const complex<float> offset, const std::vector<std::vector<complex<float>>> &Zp, volatile progress_t *progress, volatile bool *running);
template void hybrid_render(coord_t frame, map &out, stats &sta, const phybrid &H, const std::vector<blasR2<double>> &bla, const count_t subframe, const param &par, const double Zoom, const complex<double> offset, const std::vector<std::vector<complex<double>>> &Zp, volatile progress_t *progress, volatile bool *running);
template void hybrid_render(coord_t frame, map &out, stats &sta, const phybrid &H, const std::vector<blasR2<long double>> &bla, const count_t subframe, const param &par, const long double Zoom, const complex<long double> offset, const std::vector<std::vector<complex<long double>>> &Zp, volatile progress_t *progress, volatile bool *running);
template void hybrid_render(coord_t frame, map &out, stats &sta, const phybrid &H, const std::vector<blasR2<floatexp>> &bla, const count_t subframe, const param &par, const floatexp Zoom, const complex<floatexp> offset, const std::vector<std::vector<complex<floatexp>>> &Zp, volatile progress_t *progress, volatile bool *running);
template void hybrid_render(coord_t frame, map &out, stats &sta, const phybrid &H, const std::vector<blasR2<softfloat>> &bla, const count_t subframe, const param &par, const softfloat Zoom, const complex<softfloat> offset, const std::vector<std::vector<complex<softfloat>>> &Zp, volatile progress_t *progress, volatile bool *running);
#ifdef HAVE_FLOAT128
template void hybrid_render(coord_t frame, map &out, stats &sta, const phybrid &H, const std::vector<blasR2<float128>> &bla, const count_t subframe, const param &par, const float128 Zoom, const complex<float128> offset, const std::vector<std::vector<complex<float128>>> &Zp, volatile progress_t *progress, volatile bool *running);
#endif

template <typename t>
count_t hybrid_period(const phybrid &H, const std::vector<std::vector<complex<t>>> &Zp, const complex<floatexp> &c0, const count_t &Iterations, const count_t &ReferencePeriod, const floatexp &s, const mat2<double> &K, volatile progress_t *progress, volatile bool *running)
{
  complex<floatexp> C (floatexp(Zp[0][1].x), floatexp(Zp[0][1].y)); // FIXME
  dual<2, floatexp> cx (c0.x); cx.dx[0] = 1; cx.dx[1] = 0;
  dual<2, floatexp> cy (c0.y); cy.dx[0] = 0; cy.dx[1] = 1;
  const complex<dual<2, floatexp>> c (cx, cy);
  complex<dual<2, floatexp>> z (0);
  const mat2<floatexp> K1 (inverse(mat2<floatexp>(floatexp(K.x[0][0]), floatexp(K.x[0][1]), floatexp(K.x[1][0]), floatexp(K.x[1][1]))));
  floatexp Zz2 = 0;
  const floatexp r2 = e10(1, 10000);
  bool p = true;
  count_t i = 0;
  count_t n = 0;
  count_t m = 0;
  count_t phase = 0;
  while (i < Iterations && Zz2 < r2 && p && *running)
  {
    // progress
    progress[0] = i / progress_t(Iterations);
    // formula
    if (! (m < count_t(Zp[phase].size())))
    {
      break;
    }
    {
      complex<floatexp> Z(floatexp(Zp[phase][m].x), floatexp(Zp[phase][m].y));
      z = hybrid_perturb(H.per[n % H.per.size()], C, Z, c, z);
    }
    m++;
    n++;
    // rebase
    if (m == ReferencePeriod)
    {
      phase = (phase + m) % Zp.size();
      m = 0;
    }
    if (! (m < count_t(Zp[phase].size())))
    {
      break;
    }
    complex<floatexp> Z(floatexp(Zp[phase][m].x), floatexp(Zp[phase][m].y));
    const complex<dual<2, floatexp>> Zz = Z + z;
    Zz2 = norm(complex<floatexp>(Zz.x.x, Zz.y.x));
    const floatexp z2 = norm(complex<floatexp>(z.x.x, z.y.x));
    if (Zz2 < z2)
    {
      z = Zz;
      phase = (phase + m) % Zp.size();
      m = 0;
    }
    if (ReferencePeriod == 0 && m + 1 == count_t(Zp[phase].size()))
    {
      z = Zz;
      phase = (phase + m) % Zp.size();
      m = 0;
    }
    // (u1 v1) = s^{-1} K^{-1} J^{-1} (u0 v0)
    const mat2<floatexp> J(z.x.dx[0], z.x.dx[1], z.y.dx[0], z.y.dx[1]);
    complex<floatexp> w = (K1 * (inverse(J) * complex<floatexp>(Zz.x.x, Zz.y.x)));
    floatexp q = floatexp(norm(w)) / (s * s);
    p = 1 <= q;
    ++i;
    {
      complex<floatexp> Z(floatexp(Zp[phase][m].x), floatexp(Zp[phase][m].y));
      Zz2 = norm(Z + complex<floatexp>(z.x.x, z.y.x));
    }
  }
  if (i == Iterations || ! (Zz2 < r2) || p || ! *running)
  {
    return 0;
  }
  return i;
}

template count_t hybrid_period(const phybrid &H, const std::vector<std::vector<complex<float>>> &Zp, const complex<floatexp> &c0, const count_t &N, const count_t &ReferencePeriod, const floatexp &s, const mat2<double> &K, volatile progress_t *progress, volatile bool *running);
template count_t hybrid_period(const phybrid &H, const std::vector<std::vector<complex<double>>> &Zp, const complex<floatexp> &c0, const count_t &N, const count_t &ReferencePeriod, const floatexp &s, const mat2<double> &K, volatile progress_t *progress, volatile bool *running);
template count_t hybrid_period(const phybrid &H, const std::vector<std::vector<complex<long double>>> &Zp, const complex<floatexp> &c0, const count_t &N, const count_t &ReferencePeriod, const floatexp &s, const mat2<double> &K, volatile progress_t *progress, volatile bool *running);
template count_t hybrid_period(const phybrid &H, const std::vector<std::vector<complex<floatexp>>> &Zp, const complex<floatexp> &c0, const count_t &N, const count_t &ReferencePeriod, const floatexp &s, const mat2<double> &K, volatile progress_t *progress, volatile bool *running);
template count_t hybrid_period(const phybrid &H, const std::vector<std::vector<complex<softfloat>>> &Zp, const complex<floatexp> &c0, const count_t &N, const count_t &ReferencePeriod, const floatexp &s, const mat2<double> &K, volatile progress_t *progress, volatile bool *running);
#ifdef HAVE_FLOAT128
template count_t hybrid_period(const phybrid &H, const std::vector<std::vector<complex<float128>>> &Zp, const complex<floatexp> &c0, const count_t &N, const count_t &ReferencePeriod, const floatexp &s, const mat2<double> &K, volatile progress_t *progress, volatile bool *running);
#endif

bool hybrid_center(const phybrid &h, complex<mpreal> &C0, const count_t period, volatile progress_t *progress, volatile bool *running)
{
  mpfr_prec_t prec = std::max(mpfr_get_prec(C0.x.mpfr_srcptr()), mpfr_get_prec(C0.y.mpfr_srcptr()));
  const floatexp epsilon2 = floatexp(1, 16 - 2 * prec);
  double lepsilon2 = double(log(epsilon2));
  double ldelta0 = 0;
  double ldelta1 = 0;
  progress_t eta = 0;
  bool converged = false;
  const count_t maxsteps = 64; // FIXME
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
      z = hybrid_plain(h.per[i % h.per.size()], c, z);
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

bool hybrid_size(floatexp &s, mat2<double> &K, const phybrid &h, const complex<mpreal> &C0, count_t period, volatile progress_t *progress, volatile bool *running)
{
  using std::abs, ::abs;
  using std::exp, ::exp;
  using std::log, ::log;
  using std::sqrt, ::sqrt;
  double log_degree = log(double(h.per[0].power));
  complex<dual<2, mpreal>> C(C0.x, C0.y);
  C.x.dx[0] = 0;
  C.y.dx[1] = 0;
  complex<dual<2, mpreal>> Z(C);
  Z.x.dx[0] = 1;
  Z.y.dx[1] = 1;
  mat2<floatexp> b (1);
  count_t j = 1;
  count_t m = 1;
  if (m == period)
  {
    m = 0;
  }
  while (j < period && *running)
  {
    progress[0] = j / progress_t(period);
    Z = hybrid_plain(h.per[j % h.per.size()], C, Z);
    log_degree += log(double(h.per[j % h.per.size()].power));
    const mat2<floatexp> l (floatexp(Z.x.dx[0]), floatexp(Z.x.dx[1]), floatexp(Z.y.dx[0]), floatexp(Z.y.dx[1]));
    b += inverse(l);
    ++j;
  }
  const double degree = exp(log_degree / period);
  // l^d b
  if (*running)
  {
    double d = degree / (degree - 1);
    const mat2<floatexp> l (floatexp(Z.x.dx[0]), floatexp(Z.x.dx[1]), floatexp(Z.y.dx[0]), floatexp(Z.y.dx[1]));
    const floatexp lambda = sqrt(abs(determinant(l)));
    const floatexp beta = sqrt(abs(determinant(b)));
    const floatexp llb = exp(log(lambda) * d) * beta;
    s = floatexp(1 / llb);
    b = inverse(transpose(b)) / beta;
    K = mat2<double>(double(b.x[0][0]), double(b.x[0][1]), double(b.x[1][0]), double(b.x[1][1]));
    return true;
  }
  return false;
}
