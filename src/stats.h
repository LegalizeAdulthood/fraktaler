// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "types.h"

struct stats
{
  count_t pixels;
  count_t iterations;
  count_t perturb_iterations;
  count_t bla_iterations;
  count_t bla_steps;
  count_t rebases;
  count_t minimum_iterations;
  count_t maximum_iterations;
};

void reset(stats &sta);
