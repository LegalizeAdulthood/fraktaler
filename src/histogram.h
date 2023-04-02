// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <vector>

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
  bool logdata;
};

enum neighbourhood { none = 0, four = 4, eight = 8 };

histogram histogram_de_magnitude(const image_raw &img, int bins, neighbourhood next_to_interior = none);
histogram histogram_n(const image_raw &img, int bins, count_t lower_limit, count_t upper_limit);
histogram histogram_bla(const image_raw &img, int bins, count_t limit);
histogram histogram_ptb(const image_raw &img, int bins, count_t limit);

void histogram_log2(histogram &h);
void histogram_exp2(histogram &h);
