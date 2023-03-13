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
