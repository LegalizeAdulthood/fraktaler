// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "stats.h"

void reset(stats &sta)
{
  sta.pixels = 0;
  sta.iterations = 0;
  sta.perturb_iterations = 0;
  sta.bla_iterations = 0;
  sta.bla_steps = 0;
  sta.rebases = 0;
  sta.minimum_iterations = 0x7fffFFFFffffFFFFLL;
  sta.maximum_iterations = 0;
}
