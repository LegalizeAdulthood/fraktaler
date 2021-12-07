// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <string>
#include <vector>

#include "types.h"

struct colour
{
  colour() { }
  virtual ~colour() { }
  virtual std::string name() const = 0;
  virtual vec3 rgb(const count_t &n, const vec2 &coord, const vec2 &de) const noexcept = 0;
  virtual std::string frag() const = 0;
};

extern std::vector<colour *> colours;
void colours_init();
