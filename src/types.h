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

template <int D, typename T> struct dual;

typedef int channel_bit_t;
typedef int channel_mask_t;
struct map;

struct stats;

struct param;

template <typename real> struct blaC;
template <typename real> struct blasC;

template <typename real> struct blaR2;
template <typename real> struct blasR2;

struct reference;
struct formula;
struct formulaR2;
struct formulaC;
