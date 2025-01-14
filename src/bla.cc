// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "bla.h"
#include "complex.h"
#include "float128.h"
#include "floatexp.h"
#include "hybrid.h"
#include "parallel.h"
#include "softfloat.h"

template <typename real>
void blas_init1(blasR2<real> &Bp, const std::vector<std::vector<opcode>> &opss, const std::vector<int> &degrees, const count_t phase, const std::vector<complex<real>> &Zp, const real h, const real k, const real L, volatile progress_t *progress, volatile bool *running) noexcept
{
  using std::max;
  const count_t M = Bp.M;
  std::atomic<count_t> total {0};
  parallel1d(std::thread::hardware_concurrency(), 1, M, 65536, running, [&](count_t m)
  {
    int w = (phase + m) % opss.size();
    Bp.b[0][m - 1] = hybrid_bla(opss[w], degrees[w], h, k, L, Zp[m]);
    const count_t done = total.fetch_add(1);
    progress[0] = done / progress_t(2 * M);
  });
}

template <typename real>
static void blas_merge(blasR2<real> &BLA, const real h, const real k, const real L, volatile progress_t *progress, volatile bool *running) noexcept
{
  (void) L;
  using std::abs;
  using ::abs;
  using std::max;
  using std::min;
  using std::sqrt;
  using ::sqrt;
  count_t M = BLA.M;
  count_t src = 0;
  std::atomic<count_t> total {M};
  real hk = h * k;
  for (count_t msrc = M - 1; msrc > 1; msrc = (msrc + 1) >> 1) if (*running)
  {
    count_t dst = src + 1;
    count_t mdst = (msrc + 1) >> 1;
    parallel1d(std::thread::hardware_concurrency(), 0, mdst, 65536, running, [&](coord_t m)
    {
      const count_t mx = m * 2;
      const count_t my = m * 2 + 1;
      if (my < msrc)
      {
        BLA.b[dst][m] = merge(BLA.b[src][my], BLA.b[src][mx], hk);
      }
      else
      {
        BLA.b[dst][m] = BLA.b[src][mx];
      }
      const count_t done = total.fetch_add(1);
      progress[0] = done / progress_t(2 * M);
      return 0;
    });
    src++;
  }
  progress[0] = 1;
}

template <typename real>
blasR2<real>::blasR2(const std::vector<complex<real>> &Z, const std::vector<std::vector<opcode>> &opss, const std::vector<int> &degrees, const count_t phase, const real h, const real k, const real epsL, volatile progress_t *progress, volatile bool *running)
{
  M = count_t(Z.size()) - 1;
  count_t count = M > 0;
  count_t m = M - 1;
  for ( ; m > 1; m = (m + 1) >> 1)
  {
    count++;
  }
  L = count;
  b.resize(count);
  m = M - 1;
  for (count_t ix = 0; ix < count; ++ix, m = (m + 1) >> 1)
  {
    b[ix].resize(m);
  }
  blas_init1(*this, opss, degrees, phase, Z, h, k, epsL, progress, running);
  blas_merge(*this, h, k, epsL, progress, running);
}

template <typename real>
const blaR2<real> *blasR2<real>::lookup(const count_t m, const real z2) const noexcept
{
  if (m <= 0)
  {
    return 0;
  }
  if (! (m < M))
  {
    return 0;
  }
  const blaR2<real> *ret = 0;
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

template struct blasR2<float>;
template struct blasR2<double>;
template struct blasR2<long double>;
template struct blasR2<floatexp>;
template struct blasR2<softfloat>;
#ifdef HAVE_FLOAT128
template struct blasR2<float128>;
#endif
