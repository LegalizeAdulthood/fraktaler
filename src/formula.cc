// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "formula.h"

#if 0
reference::reference(const mpfr_t Cx0, const mpfr_t Cy0)
{
  mpfr_prec_t prec = std::max(mpfr_get_prec(Cx0), mpfr_get_prec(Cy0));
  mpfr_init2(Cx, prec); mpfr_set(Cx, Cx0, MPFR_RNDN);
  mpfr_init2(Cy, prec); mpfr_set(Cy, Cy0, MPFR_RNDN);
  mpfr_init2(Zx, prec); mpfr_set_ui(Zx, 0, MPFR_RNDN);
  mpfr_init2(Zy, prec); mpfr_set_ui(Zy, 0, MPFR_RNDN);
}

reference::~reference()
{
  mpfr_clear(Cx);
  mpfr_clear(Cy);
  mpfr_clear(Zx);
  mpfr_clear(Zy);
}

complex<floatexp> reference::get() const
{
  long ex = 0, ey = 0;
  double x = mpfr_get_d_2exp(&ex, Zx, MPFR_RNDN);
  double y = mpfr_get_d_2exp(&ey, Zy, MPFR_RNDN);
  return complex<floatexp>(floatexp(x, ex), floatexp(y, ey));
}

bool reference::escaped() const
{
  return 4 < norm(get()); // FIXME hardcoded
}

#endif

#if 0
template <complex<dual<2, mpreal>> F(const complex<dual<2, mpreal>>&, const complex<dual<2, mpreal>>&)>
count_t period(const mpreal &A, const mpreal &B, const count_t N, const floatexp &S, const mat2<double> &K, progress_t *progress, bool *running)
{
  using R = dual<2, mpreal>;
  using C = complex<R>;
  R cx(A); cx.dx[0] = 1;
  R cy(B); cy.dx[1] = 1;
  C c(cx, cy);
  C z(0, 0);
  mat2<double> K1(inverse(K));
  double z2 = 0;
  double r2 = 1e50;
  bool p = true;
  count_t i = 0;
  while (i < N && z2 < r2 && p && *running)
  {
    // progress
    progress[0] = i / progress_t(N);
    // formula
    z = F(c, z);
    // unpack
    const mpreal &x = z.x.x;
    const mpreal &y = z.y.x;
    const mpreal &xa = z.x.dx[0];
    const mpreal &xb = z.x.dx[1];
    const mpreal &ya = z.y.dx[0];
    const mpreal &yb = z.y.dx[1];
    // J^{-1}
    const mpreal detJ = xa * yb - xb * ya;
    const mpreal xa1 =  yb / detJ;
    const mpreal xb1 = -ya / detJ;
    const mpreal ya1 = -xb / detJ;
    const mpreal yb1 =  xa / detJ;
    // (u0 v0) = J^{-1} * (x y)
    const mpreal u0 = xa1 * x + xb1 * y;
    const mpreal v0 = ya1 * x + yb1 * y;
    // (u1 v1) = s^{-1} K^{-1} * (u0 v0)
    const complex<mpreal> w0 (u0, v0);
    const complex<mpreal> w1 = (K1 * w0) / S;
    const double u2 = double(w1.x);
    const double v2 = double(w1.y);
    const double uv = u2 * u2 + v2 * v2;
    p = 1 <= uv;
    ++i;
    const double xf = double(x);
    const double yf = double(y);
    z2 = xf * xf + yf * yf;
  }
  if (i == N || r2 <= z2 || p || ! *running)
  {
    i = -1;
  }
  return i;
}
#endif

#if 0
template <complex<dual<2, mpreal>> F(const complex<dual<2, mpreal>>&, const complex<dual<2, mpreal>>&)>
bool center(mpreal &Cx, mpreal &Cy, const count_t period, progress_t *progress, bool *running) const
{
  using R = dual<2, mpreal>;
  using C = complex<R>;
  mpfr_prec_t prec = std::max(mpfr_get_prec(Cx.mpfr_srcptr()), mpfr_get_prec(Cy.mpfr_srcptr()));
  const floatexp epsilon2 = floatexp(1, 16 - 2 * prec);
  double lepsilon2 = double(log(epsilon2)); // FIXME slow? precision?
  double ldelta0 = 0;
  double ldelta1 = 0;
  progress_t eta = 0;
  bool converged = false;
  const count_t maxsteps = 64;
  for (count_t j = 0; j < maxsteps && *running && ! converged; ++j)
  {
    progress[0] = j / eta;
    progress[1] = 0;
    R cx(Cx); cx.dx[0] = 1;
    R cy(Cy); cy.dx[1] = 1;
    C c(cx, cy);
    C z(0, 0);
    // iteration
    for (count_t i = 0; i < period && *running; ++i)
    {
      progress[1] = i / progress_t(period);
      z = F(c, z);
    }
    if (*running)
    {
      const mpreal &x = z.x.x;
      const mpreal &y = z.y.x;
      const mpreal &dxa = z.x.dx[0];
      const mpreal &dxb = z.x.dx[1];
      const mpreal &dya = z.y.dx[0];
      const mpreal &dyb = z.Y.dx[1];
      // Newton step
      const mpreal det = dxa * dyb - dxb * dya;
      const mpreal u = -( dyb * x - dxb * y) / det;
      const mpreal v = -(-dya * x + dxa * y) / det;
      Cx += u;
      Cy += v;
      // check convergence
      floatexp uf = floatexp(u);
      floatexp vf = floatexp(uv;
      floatexp delta = sqr(uf) + sqr(vf);
      converged = delta < epsilon2;
      ldelta0 = ldelta1;
      ldelta1 = double(log(delta));
      eta = log2((lepsilon2 - ldelta0) / (ldelta1 - ldelta0));
    }
  }
  return converged;
}
#endif

#if 0
template <complex<dual<4, mpreal>> F(const complex<dual<4, mpreal>>&, const complex<dual<4, mpreal>>&)>
extern bool skew(const CDecNumber &cr, const CDecNumber &ci, bool useDZ, double *skew_matrix, volatile int *running, int *progress)
{
  Precision prec(std::max(cr.m_dec.precision(), ci.m_dec.precision()));
  using N = int;
  using R = dual<2, CDecNumber>;
  using C = complex<R>;
  int count = 0;
  int stanza = 0;
  // FIXME assumes seed is 0+0i and/or first iteration result is C
  R cx_dc(cr); cx_dc.dx[0] = 1;
  R cy_dc(ci); cy_dc.dx[1] = 1;
  C c_dc(cx_dc, cy_dc);
  C z_dc(c_dc);
  R cx_dz(cr);
  R cy_dz(ci);
  C c_dz(cx_dz, cy_dz);
  C z_dz(c_dc);
  for (N j = 0; j < maxiters && *running; ++j)
  {
    progress[0] = j;
    progress[1] = 0;
    progress[2] = 0;
    progress[3] = 0;
    if (++count >= h.stanzas[stanza].repeats)
    {
      count = 0;
      if (++stanza >= (ssize_t) h.stanzas.size())
      {
        stanza = h.loop_start;
      }
    }
    if (sqr(z_dc.m_r.x) + sqr(z_dc.m_i.x) > 65536.0)
    {
      break;
    }
    z_dc = hybrid_f(h.stanzas[stanza], z_dc, c_dc);
    if (useDZ)
    {
      z_dz = hybrid_f(h.stanzas[stanza], z_dz, c_dz);
    }
  }
  if (running)
  {
    const CDecNumber &dxa = z_dc.m_r.dx[0];
    const CDecNumber &dxb = z_dc.m_r.dx[1];
    const CDecNumber &dya = z_dc.m_i.dx[0];
    const CDecNumber &dyb = z_dc.m_i.dx[1];
    CDecNumber sii, sij, sji, sjj;
    if (useDZ)
    {
      const CDecNumber &dxx = z_dz.m_r.dx[0];
      const CDecNumber &dxy = z_dz.m_r.dx[1];
      const CDecNumber &dyx = z_dz.m_i.dx[0];
      const CDecNumber &dyy = z_dz.m_i.dx[1];
      sii = dyb * dxx - dxb * dyx;
      sij = dyb * dxy - dxb * dyy;
      sji = dxa * dyx - dya * dxx;
      sjj = dxa * dyy - dya * dxy;
    }
    else
    {
      sii =   dyb;
      sij = - dxb;
      sji = - dya;
      sjj =   dxa;
    }
    CDecNumber det = sqrt(abs(sii * sjj - sij * sji));
    sii = sii / det;
    sij = sij / det;
    sji = sji / det;
    sjj = sjj / det;
    skew_matrix[0] = double(sii);
    skew_matrix[1] = double(sij);
    skew_matrix[2] = double(sji);
    skew_matrix[3] = double(sjj);
    return true;
  }
  return false;
}
#endif

#if 0
template <complex<dual<2, mpreal>> F(const complex<dual<2, mpreal>>&, const complex<dual<2, mpreal>>&)>
bool size(const mpreal &A, const mpreal &B, const double degree, floatexp &S, mat2<double> &K, progress_t *progress, bool *running)
{
  using R = dual<2, mpreal>;
  using C = complex<R>;
  R x(A); x.dx[0] = 1;
  R y(B); y.dx[1] = 1;
  C z(x, y);
  C c(A, B);
  floatexp bxa = 1;
  floatexp bxb = 0;
  floatexp bya = 0;
  floatexp byb = 1;
  count_t j = 1;
  while (j < period && *running)
  {
    progress[0] = j / progress_t(period);
    z = F(c, z);
    const floatexp lxa = floatexp(z.x.dx[0]);
    const floatexp lxb = floatexp(z.x.dx[1]);
    const floatexp lya = floatexp(z.y.dx[0]);
    const floatexp lyb = floatexp(z.y.dx[1]);
    // step b
    const floatexp det = lxa * lyb - lxb * lya;
    bxa = bxa + lyb / det;
    bxb = bxb - lxb / det;
    bya = bya - lya / det;
    byb = byb + lxa / det;
    ++j;
  }
  double deg = degree / (degree - 1);
  if (isnan(deg) || isinf(deg)) deg = 0;
  // l^d b
  if (*running)
  {
    const floatexp lxa = floatexp(z.x.dx[0]);
    const floatexp lxb = floatexp(z.x.dx[1]);
    const floatexp lya = floatexp(z.y.dx[0]);
    const floatexp lyb = floatexp(z.y.dx[1]);
    const floatexp l = sqrt(abs(lxa * lyb - lxb * lya));
    const floatexp beta = sqrt(abs(bxa * byb - bxb * bya));
    const floatexp llb = exp(log(l) * deg) * beta;
    S = 1 / llb;
    byb =   byb / beta;
    bxb = - bxb / beta;
    bya = - bya / beta;
    bxa =   bxa / beta;
    K = mat2<double(byb, bxb, bya, bxa);
    return true;
  }
  return false;
}
#endif

#if 0
template <complex<dual<2, mpreal>> F(const complex<dual<2, mpreal>>&, const complex<dual<2, mpreal>>&)>
bool domain_size(const mpreal &A, const mpreal &B, const count_t period, floatexp &S, progress_t *progress, bool *running)
{
  using R = dual<2, mpreal>;
  using C = complex<R>;
  R x(A); x.dx[0] = 1;
  R y(B); y.dx[1] = 1;
  C c(x, y);
  C z(c);
  count_t j = 2;
  int count = 0;
  int stanza = 0;
  floatexp zq2 = sqr(floatexp(z.x.x)) + sqr(floatexp(z.y.x));
  while (j <= period && *running)
  {
    // progress
    progress[0] = j / progress_t(period);
    // formula
    z = F(c, z);
    // capture penultimate minimum |z|
    floatexp zp2 = sqr(floatexp(z.x.x)) + sqr(floatexp(z.y.x));
    if (j < period && zp2 < zq2)
    {
      zq2 = zp2;
    }
    ++j;
  }
  if (*running)
  {
    const floatexp lxa = floatexp(z.x.dx[0]);
    const floatexp lxb = floatexp(z.x.dx[1]);
    const floatexp lya = floatexp(z.y.dx[0]);
    const floatexp lyb = floatexp(z.y.dx[1]);
    const floatexp det = lxa * lyb - lxb * lya;
    S = sqrt(zq2) / sqrt(abs(det));
    return true;
  }
  return false;
}
#endif

#include "formula_mandelbrot.h"
#include "formula_burningship.h"

std::vector<formula *> formulas;

void formulas_init()
{
  #define C(stem) formulas.push_back(new formulaC< \
    stem ## _name, \
    stem ## _plain <complex<mpreal>>, \
    stem ## _plain <dual<1, complex<mpreal>>>, \
    stem ## _perturb <complex<float>, dual<1, complex<float>>>, \
    stem ## _perturb <complex<double>, dual<1, complex<double>>>, \
    stem ## _perturb <complex<long double>, dual<1, complex<long double>>>, \
    stem ## _perturb <complex<floatexp>, dual<1, complex<floatexp>>>, \
    stem ## _bla <float>, \
    stem ## _bla <double>, \
    stem ## _bla <long double>, \
    stem ## _bla <floatexp> \
  >());
  #define R2(stem) formulas.push_back(new formulaR2< \
    stem ## _name, \
    stem ## _plain <mpreal>, \
    stem ## _plain <dual<2, mpreal>>, \
    stem ## _perturb <float, dual<2, float>>, \
    stem ## _perturb <double, dual<2, double>>, \
    stem ## _perturb <long double, dual<2, long double>>, \
    stem ## _perturb <floatexp, dual<2, floatexp>>, \
    stem ## _bla <float>, \
    stem ## _bla <double>, \
    stem ## _bla <long double>, \
    stem ## _bla <floatexp> \
  >());
  C(mandelbrot)
  R2(burningship)
}
