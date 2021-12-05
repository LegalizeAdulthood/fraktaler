// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "bla.h"
#include "complex.h"
#include "floatexp.h"
#include "formula.h"

template <typename real>
static void blas_init1(blas<real> *BLA, const formulaC *formula, const complex<real> *Z, const real h, const real k, real L, progress_t *progress, bool *running)
{
  using std::max;
  const count_t M = BLA->M;
  count_t total = 0;
  #pragma omp parallel for
  for (count_t m = 1; m < M; ++m) if (running)
  {
    BLA->b[0][m - 1] = formula->bla1(h, k, L, Z[m]);
    count_t done;
    #pragma omp atomic capture
    done = total++;
    progress[0] = done / progress_t(2 * M);
  }
}

template <typename real>
static void blas_merge(blas<real> *BLA, const real h, const real k, const real L, progress_t *progress, bool *running)
{
  (void) L;
  using std::abs;
  using std::max;
  using std::min;
  using std::sqrt;
  count_t M = BLA->M;
  count_t src = 0;
  count_t total = M;
  for (count_t msrc = M - 1; msrc > 1; msrc = (msrc + 1) >> 1) if (*running)
  {
    count_t dst = src + 1;
    count_t mdst = (msrc + 1) >> 1;
    #pragma omp parallel for
    for (count_t m = 0; m < mdst; ++m) if (*running)
    {
      const count_t mx = m * 2;
      const count_t my = m * 2 + 1;
      if (my < msrc)
      {
        const bla<real> x = BLA->b[src][mx];
        const bla<real> y = BLA->b[src][my];
        const count_t l = x.l + y.l;
        const complex<real> A = y.A * x.A;
        const complex<real> B = y.A * x.B + y.B;
        const real xA = abs(x.A);
        const real xB = abs(x.B);
        const real r = min(sqrt(x.r2), max(real(0), (sqrt(y.r2) - xB * h * k) / xA));
        const real r2 = r * r;
        bla<real> b = { A, B, r2, l };
        BLA->b[dst][m] = b;
      }
      else
      {
        BLA->b[dst][m] = BLA->b[src][mx];
      }
      count_t done;
      #pragma omp atomic capture
      done = total++;
      progress[0] = done / progress_t(2 * M);
    }
    src++;
  }
  progress[0] = 1;
}

template <typename real>
blas<real>::blas(const count_t M0, const complex<real> *Z, const formulaC *formula, const real h, const real k, const real stepcount, progress_t *progress, bool *running)
{
  M = M0;
  count_t total = 1;
  count_t count = 1;
  count_t m = M - 1;
  for ( ; m > 1; m = (m + 1) >> 1)
  {
    total += m;
    count++;
  }
  L = count;
  b = new bla<real> *[count];
  b[0] = new bla<real>[total];
  count_t ix = 1;
  m = M - 1;
  for ( ; m > 1; m = (m + 1) >> 1)
  {
    b[ix] = b[ix - 1] + m;
    ix++;
  }
  blas_init1(this, formula, Z, h, k, stepcount, progress, running);
  blas_merge(this, h, k, stepcount, progress, running);
}

template <typename real>
const bla<real> *blas<real>::lookup(const count_t m, const real z2) const
{
  if (m <= 0)
  {
    return 0;
  }
  if (! (m < M))
  {
    return 0;
  }
  const bla<real> *ret = 0;
  count_t ix = m - 1;
  for (count_t level = 0; level < L; ++level)
  {
    count_t ixm = (ix << level) + 1;
    if (m == ixm && z2 < b[level][ix].r2)
    {
      ret = &b[level][ix];
    }
    else
    {
      break;
    }
    ix = ix >> 1;
  }
  return ret;
}

template struct blas<float>;
template struct blas<double>;
template struct blas<long double>;
template struct blas<floatexp>;
