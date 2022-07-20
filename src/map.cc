// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#if 0

// just for square root
#include <complex>

#include "map.h"
#include "matrix.h"

#ifdef HAVE_EXR
#include <ImfNamespace.h>
#include <ImfOutputFile.h>
#include <ImfHeader.h>
#include <ImfChannelList.h>
#include <ImfIntAttribute.h>
#include <ImfStringAttribute.h>
#include <ImfArray.h>
#include <ImfFrameBuffer.h>

namespace IMF = OPENEXR_IMF_NAMESPACE;
using namespace IMF;
using namespace IMATH_NAMESPACE;

const char kf2plus[] = "KallesFraktaler2+";
const char fraktaler3[] = "Fraktaler3";
#endif

void map::saveEXR(const std::string &filename, const channel_mask_t channels, const int threads, const std::string &metadata, const std::string &kf2plus_metadata) const
{
#ifdef HAVE_EXR
  setGlobalThreadCount(threads);
  // prepare preview image
  Header header(width, height);
  // insert metadata
  if (metadata != "") header.insert(fraktaler3, StringAttribute(metadata));
  if (kf2plus_metadata != "") header.insert(kf2plus, StringAttribute(kf2plus_metadata));
  if (maxiters + Nbias < INT_MAX - 1)
  {
    header.insert("Iterations", IntAttribute(maxiters));
  }
  else
  {
    std::ostringstream s;
    s << maxiters;
    header.insert("Iterations", StringAttribute(s.str()));
  }
  if (Nbias < INT_MAX - 1)
  {
    header.insert("IterationsBias", IntAttribute(Nbias));
  }
  else
  {
    std::ostringstream s;
    s << Nbias;
    header.insert("IterationsBias", StringAttribute(s.str()));
  }
  // write image
  if (RGB && (channels & (1 << Channel_R))) header.channels().insert("R", Channel(IMF::HALF));
  if (RGB && (channels & (1 << Channel_G))) header.channels().insert("G", Channel(IMF::HALF));
  if (RGB && (channels & (1 << Channel_B))) header.channels().insert("B", Channel(IMF::HALF));
  bool twoN = maxiters + Nbias >= 0xFFffFFfeU;
  if (N0 && (channels & (1 << Channel_N0)) &&
      N1 && (channels & (1 << Channel_N1)) &&
      twoN)
  {
    header.channels().insert("N0",  Channel(IMF::UINT));
    header.channels().insert("N1",  Channel(IMF::UINT));
  }
  else if (N0 && (channels & (1 << Channel_N0)))
  {
    header.channels().insert("N",  Channel(IMF::UINT));
  }
  if (NF  && (channels & (1 << Channel_NF ))) header.channels().insert("NF",  Channel(IMF::FLOAT));
  if (T   && (channels & (1 << Channel_T  ))) header.channels().insert("T",   Channel(IMF::FLOAT));
  if (DEX && (channels & (1 << Channel_DEX))) header.channels().insert("DEX", Channel(IMF::FLOAT));
  if (DEY && (channels & (1 << Channel_DEY))) header.channels().insert("DEY", Channel(IMF::FLOAT));
  OutputFile of(filename.c_str(), header);
  FrameBuffer fb;
  if (RGB && (channels & (1 << Channel_R))) fb.insert("R", Slice(IMF::HALF, (char *)(RGB + 0), sizeof(*RGB) * 3, sizeof(*RGB) * 3 * width));
  if (RGB && (channels & (1 << Channel_G))) fb.insert("G", Slice(IMF::HALF, (char *)(RGB + 1), sizeof(*RGB) * 3, sizeof(*RGB) * 3 * width));
  if (RGB && (channels & (1 << Channel_B))) fb.insert("B", Slice(IMF::HALF, (char *)(RGB + 2), sizeof(*RGB) * 3, sizeof(*RGB) * 3 * width));
  if (N0 && (channels & (1 << Channel_N0)) && N1 && (channels & (1 << Channel_N1)) && twoN)
  {
    fb.insert("N0", Slice(IMF::UINT, (char *)(N0), sizeof(*N0), sizeof(*N0) * width));
    fb.insert("N1", Slice(IMF::UINT, (char *)(N1), sizeof(*N1), sizeof(*N1) * width));
  }
  else if (N0 && (channels & (1 << Channel_N0)))
  {
    fb.insert("N", Slice(IMF::UINT, (char *)(N0), sizeof(*N0), sizeof(*N0) * width));
  }

  if (NF  && (channels & (1 << Channel_NF ))) fb.insert("NF",  Slice(IMF::FLOAT, (char *)(NF ), sizeof(*NF ), sizeof(*NF ) * width));
  if (T   && (channels & (1 << Channel_T  ))) fb.insert("T",   Slice(IMF::FLOAT, (char *)(T  ), sizeof(*T  ), sizeof(*T  ) * width));
  if (DEX && (channels & (1 << Channel_DEY))) fb.insert("DEX", Slice(IMF::FLOAT, (char *)(DEX), sizeof(*DEX), sizeof(*DEX) * width));
  if (DEY && (channels & (1 << Channel_DEY))) fb.insert("DEY", Slice(IMF::FLOAT, (char *)(DEY), sizeof(*DEY), sizeof(*DEY) * width));
  of.setFrameBuffer(fb);
  of.writePixels(height);
#endif
}

mat2<double> map::unskewDE(void) const
{
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

#endif

