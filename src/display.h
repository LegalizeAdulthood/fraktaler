// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "types.h"

struct display
{
  display() { }
  virtual ~display() { }
  virtual void resize(coord_t width, coord_t height) = 0;
  virtual void set_colour(const colour *clr) = 0;
  virtual void clear() = 0;
  virtual void accumulate(const map &out) = 0;
  virtual void get_rgb(map &out) const = 0;
};
