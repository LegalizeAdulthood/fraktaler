// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <vector>

#include "display.h"

struct display_cpu : public display
{
  coord_t width;
  coord_t height;
  bool do_clear;
  count_t subframes;
  std::vector<vec3> RGB;
  const colour *clr;

  display_cpu(const colour *clr);
  virtual ~display_cpu();

  virtual void resize(coord_t width, coord_t height);
  virtual void set_colour(const colour *clr);
  virtual void clear();
  virtual void accumulate(const map &out);
  virtual void get_rgb(map &out) const;
};
