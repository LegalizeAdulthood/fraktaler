// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <cstdint>

typedef int64_t coord_t;

typedef float smooth_t;
typedef int64_t count_t;

typedef float progress_t;

typedef float mantissa;
typedef int exponent;
struct floatexp;

template <typename real> struct complex;

struct map;
struct param;