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

histogram histogram_n(const image_raw &img, int bins, count_t lower_limit, count_t upper_limit);
histogram histogram_bla(const image_raw &img, int bins, count_t limit);
histogram histogram_ptb(const image_raw &img, int bins, count_t limit);

void histogram_log2(histogram &h);
void histogram_exp2(histogram &h);

struct histogram2d
{
  int width;
  int height;
  float total;
  std::vector<float> data;
  float peak;
  bool logdata;
};

histogram2d histogram_logde(const image_raw &img);

void histogram2d_log2(histogram2d &h);
void histogram2d_exp2(histogram2d &h);
