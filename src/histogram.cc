// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

// just for square root
#include <complex>

#include "histogram.h"
#include "image_raw.h"
#include "matrix.h"

mat2<double> unskew_de(const image_raw &img)
{
  const float *DEX = img.DEX;
  const float *DEY = img.DEY;
  const coord_t width = img.width;
  const coord_t height = img.height;
  if (! (DEX && DEY))
  {
    return mat2<double>(1, 0, 0, 1);
  }
  // compute statistics pass one: mean and mean(C^2)
  const double squared_length_threshold = 4;
  glm::dvec2 mean(0.0, 0.0);
  glm::dvec2 mean_csqr(0.0, 0.0);
  glm::dmat2 covariance(0.0, 0.0, 0.0, 0.0);
  double count = 0;
  for (int j = 0; j < height; ++j)
  {
    for (int i = 0; i < width; ++i)
    {
      int ix = j * width + i;
      double x = DEX[ix];
      double y = DEY[ix];
      double squared_length = x * x + y * y;
      if (squared_length < squared_length_threshold)
      {
        if (0 < squared_length)
        {
          mean[0] += x;
          mean[1] += y;
          count += 1;
          // complex squaring of DE
          // makes limit distribution of hard-skewed locations unipolar
          double t = x * x - y * y;
          y = 2 * x * y;
          x = t;
          mean_csqr[0] += x;
          mean_csqr[1] += y;
        }
      }
    }
  }
  mean[0] /= count;
  mean[1] /= count;
  mean_csqr[0] /= count;
  mean_csqr[1] /= count;
  // compute statistics pass 2: covariance
  for (int j = 0; j < height; ++j)
  {
    for (int i = 0; i < width; ++i)
    {
      int ix = j * width + i;
      double x = DEX[ix];
      double y = DEY[ix];
      double squared_length = x * x + y * y;
      if (squared_length < squared_length_threshold)
      {
        if (0 < squared_length)
        {
          double dx = x - mean[0];
          double dy = y - mean[1];
          covariance[0][0] += dx * dx;
          covariance[0][1] += dx * dy;
          covariance[1][1] += dy * dy;
        }
      }
    }
  }
  covariance[0][0] /= count;
  covariance[0][1] /= count;
  covariance[1][1] /= count;
  covariance[1][0] = covariance[0][1];
  // get mean of complex de^2
  double x = mean_csqr.x;
  double y = mean_csqr.y;
  // the branch cut in square root doesn't matter
  // stretching up is the same as stretching down
  std::complex<double> dir = std::sqrt(std::complex<double>(x, y));
  // normalize
  x = std::real(dir);
  y = std::imag(dir);
  double d = 1 / std::sqrt(x * x + y * y);
  x *= d;
  y *= d;
  // formulate skew matrix
  double stretch = std::sqrt(eigenvalue_ratio(covariance));
  mat2<double> P(x, -y, y, x);
  mat2<double> P1(x, y, -y, x);
  mat2<double> S(1 / stretch, 0, 0, stretch);
  // no need for Jacobians, because it is a directional DE to a shape
  mat2<double> Q(P1 * S * P);
  return Q;
}

histogram histogram_uint(const uint32_t *data, const image_raw &img, int bins, count_t limit)
{
  histogram h = { 0, double(limit), false, 0, { }, false };
  h.data.resize(bins);
  std::fill(h.data.begin(), h.data.end(), 0.0f);
  if (! data)
  {
    return h;
  }
  const coord_t width = img.width;
  const coord_t height = img.height;
  double range = bins / (h.maximum - h.minimum);
  for (coord_t j = 0; j < height; ++j)
  {
    for (coord_t i = 0; i < width; ++i)
    {
      int ix = j * width + i;
      double x = data[ix];
      int bin = (x - h.minimum) * range;
      bin = std::min(std::max(bin, 0), bins - 1);
      h.data[bin] += 1.0f;
      h.total += 1.0f;
    }
  }
  return h;
}

histogram histogram_bla(const image_raw &img, int bins, count_t limit)
{
  return histogram_uint(img.BLA, img, bins, limit);
}

histogram histogram_ptb(const image_raw &img, int bins, count_t limit)
{
  return histogram_uint(img.PTB, img, bins, limit);
}

histogram histogram_n(const image_raw &img, int bins, count_t lower_limit, count_t upper_limit)
{
  histogram h = { double(lower_limit), double(upper_limit), false, 0, { }, false };
  h.data.resize(bins);
  std::fill(h.data.begin(), h.data.end(), 0.0f);
  if (! img.N0)
  {
    return h;
  }
  const coord_t width = img.width;
  const coord_t height = img.height;
  for (coord_t j = 0; j < height; ++j)
  {
    for (coord_t i = 0; i < width; ++i)
    {
      int ix = j * width + i;
      count_t n = img.N0[ix];
      if (img.N1)
      {
        n |= count_t(img.N1[ix]) << 32;
      }
      else if (n == count_t(~uint32_t(0)))
      {
        n = ~count_t(0);
      }
      if (n != ~count_t(0))
      {
        h.minimum = std::min(h.minimum, double(n));
      }
    }
  }
  double range = bins / (h.maximum - h.minimum);
  for (coord_t j = 0; j < height; ++j)
  {
    for (coord_t i = 0; i < width; ++i)
    {
      int ix = j * width + i;
      count_t n = img.N0[ix];
      if (img.N1)
      {
        n |= count_t(img.N1[ix]) << 32;
      }
      else if (n == count_t(~uint32_t(0)))
      {
        n = ~count_t(0);
      }
      if (n != ~count_t(0))
      {
        int bin = (n - h.minimum) * range;
        bin = std::min(std::max(bin, 0), bins - 1);
        h.data[bin] += 1.0f;
        h.total += 1.0f;
      }
    }
  }
  return h;
}

void histogram_log2(histogram &h)
{
  if (h.logdata)
  {
    return;
  }
  for (auto & b : h.data)
  {
    if (b >= 1)
    {
      b = 1 + std::log2(b);
    }
  }
  h.logdata = true;
}

void histogram_exp2(histogram &h)
{
  if (! h.logdata)
  {
    return;
  }
  for (auto & b : h.data)
  {
    if (b >= 1)
    {
      b = std::exp2(b - 1);
    }
  }
  h.logdata = false;
}

histogram2d histogram_logde(const image_raw &img)
{
  histogram2d h = { 32, 18, 0.0f, { }, 0.0f, false };
  h.data.resize(h.width * h.height);
  std::fill(h.data.begin(), h.data.end(), 0.0f);
  const float *DEX = img.DEX;
  const float *DEY = img.DEY;
  if (! (DEX && DEY)) return h;
  const float xmin = std::log(1.0e-24f);
  const float xmax = std::log(1.0e6f);
  const float xmul = h.width / (xmax - xmin);
  const float ymin = -3.141592653f;
  const float ymax = 3.141592653f;
  const float ymul = h.height / (ymax - ymin);
  const coord_t width = img.width;
  const coord_t height = img.height;
  for (coord_t j = 0; j < height; ++j)
  {
    for (coord_t i = 0; i < width; ++i)
    {
      int ix = j * width + i;
      const float x = DEX[ix];
      const float y = DEY[ix];
      const float de2 = x * x + y * y;
      if (de2 > 0)
      {
        int xbin = std::min(std::max(xmul * (std::log(de2) - xmin), 0.0f), h.width - 1.0f);
        int ybin = std::min(std::max(ymul * (std::atan2(y, x) - ymin), 0.0f), h.height - 1.0f);
        int bin = ybin * h.width + xbin;
        h.data[bin] += 1.0f;
        h.total += 1.0f;
      }
    }
  }
  h.peak = 0;
  for (const auto d : h.data)
  {
    h.peak = std::max(h.peak, d);
  }
  return h;
}

void histogram2d_log2(histogram2d &h)
{
  if (h.logdata)
  {
    return;
  }
  for (auto & b : h.data)
  {
    if (b >= 1)
    {
      b = 1 + std::log2(b);
    }
  }
  if (h.peak >= 1)
  {
    h.peak = 1 + std::log2(h.peak);
  }
  h.logdata = true;
}

void histogram2d_exp2(histogram2d &h)
{
  if (! h.logdata)
  {
    return;
  }
  for (auto & b : h.data)
  {
    if (b >= 1)
    {
      b = std::exp2(b - 1);
    }
  }
  if (h.peak >= 1)
  {
    h.peak = std::exp2(h.peak - 1);
  }
  h.logdata = false;
}
