// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "floatexp.h"
#include "newton.h"

template <typename real>
count_t find_period(const complex<real> *Zp, const count_t M, const complex<real> c, const count_t N, const real r, progress_t *progress, bool *running)
{
  // perturbed version of knighty's Taylor ball with Zhuoran's rebasing
  complex<real> z(0), dz(0), Z(0), Zz(0);
  real rdz = 0, rz = 0, rzd = 0, rr = 0, rrdz = 0, Ei = 0, ER = 65536;
  count_t n = 0, m = 0;
  while (n < N && rz < ER)
  {
    progress[0] = n / progress_t(N);
    if (! *running)
    {
      break;
    }
    // step
    Ei = rdz * rdz + (2 * rz + r * (2 * rdz + r * Ei)) * Ei;
    dz = 2 * Zz * dz + 1;
    z = (2 * Z + z) * z + c;
    n++;
    m++;
    Z = Zp[m];
    Zz = Z + z;
    // calculate ball
    rdz = abs(dz);
    rz = abs(Zz);
    rzd = abs(z);
    rr = r * (rdz + r * Ei);
    rrdz = r * (rdz - r  * Ei);
    // rebase
    if (rz < rzd || m == M - 1)
    {
      z = Zz;
      m = 0;
      Z = Zp[m];
    }
    if (rz - rr > 2)
    {
      // escaped
      break;
    }
    if (rz <= rr)
    {
      count_t period = n;
      if (! (rz <= rrdz))
      {
        period = -period;
      }
      return period;
    }    
  }
  return 0;
}

template count_t find_period(const complex<float> *Zp, const count_t M, const complex<float> c, const count_t N, const float r, progress_t *progress, bool *running);
template count_t find_period(const complex<double> *Zp, const count_t M, const complex<double> c, const count_t N, const double r, progress_t *progress, bool *running);
template count_t find_period(const complex<long double> *Zp, const count_t M, const complex<long double> c, const count_t N, const long double r, progress_t *progress, bool *running);
template count_t find_period(const complex<floatexp> *Zp, const count_t M, const complex<floatexp> c, const count_t N, const floatexp r, progress_t *progress, bool *running);
