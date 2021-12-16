// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <filesystem>

#include <mpreal.h>

#include "complex.h"
#include "floatexp.h"
#include "matrix.h"
#include "types.h"

struct param
{
  complex<mpreal> C;
  floatexp Zoom;
  count_t Iterations;
  count_t MaxRefIters;
  count_t MaxPtbIters;
  double EscapeRadius;
  count_t ReferencePeriod;
  bool LockMaxRefItersToPeriod;
  bool ReuseReference;
  bool ReuseBLA;
  count_t MaxSubframes;
  bool ExponentialMap;
  bool ZoomOutSequence;
  channel_mask_t Channels;
  std::string Stem;
  coord_t Width;
  coord_t Height;
  mat2<double> K;
  std::string sRe;
  std::string sIm;
  std::string sZoom;
  std::string sIterations;
  std::string sMaxRefIters;
  std::string sMaxPtbIters;
  std::string sEscapeRadius;
  param();
};

void restring(param &par);
void home(param &par);
void zoom(param &par, double x, double y, double g, bool fixed_click = true);
void zoom(param &par, const mat3 &T, const mat3 &T0);
