// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <filesystem>

#include <mpfr.h>

#include "floatexp.h"
#include "types.h"

struct param
{
  mpfr_t Cx;
  mpfr_t Cy;
  floatexp Zoom;
  count_t Iterations;
  count_t ReferencePeriod;
  count_t MaximumReferenceIterations;
  count_t PerturbIterations;
  bool ExponentialMap;
  bool ZoomOutSequence;
  channel_mask_t Channels;
  std::filesystem::path Stem;
  coord_t Width;
  coord_t Height;
};

void home(param &par);
void zoom(param &par, double x, double y, double g, bool fixed_click = true);
