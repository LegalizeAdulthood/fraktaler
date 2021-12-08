// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "colour.h"
#include "colour_monochrome.h"
#include "colour_rainbow.h"

std::vector<colour *> colours;

void colours_init()
{
  colours.push_back(new colour_monochrome());
  colours.push_back(new colour_rainbow());
}
