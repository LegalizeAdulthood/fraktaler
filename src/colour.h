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

inline float linear_to_srgb(float c) noexcept
{
  c = std::min(std::max(c, 0.0f), 1.0f);
  if (c <= 0.0031308f)
  {
    return 12.92f * c;
  }
  else
  {
    return 1.055f * std::pow(c, 1.0f / 2.4f) - 0.055f;
  }
}

extern std::vector<colour *> colours;
void colours_init();
