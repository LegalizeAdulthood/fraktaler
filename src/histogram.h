// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <vector>

#include "types.h"

struct image_raw;

mat2<double> unskew_de(const image_raw &img);

struct histogram
{
  double minimum;
  double maximum;
  bool logarithmic;
  float total;
  std::vector<float> data;
};

enum neighbourhood { none = 0, four = 4, eight = 8 };

histogram histogram_de_magnitude(const image_raw &img, int bins, neighbourhood next_to_interior = none);
