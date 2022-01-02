// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "colour.h"

struct colour_monochrome : public colour
{
  colour_monochrome() { }
  virtual ~colour_monochrome() { }
  virtual const char *name() const { return "Monochrome"; }
  virtual vec3 rgb(const count_t &n, const vec2 &coord, const vec2 &de) const noexcept
  {
    (void) n;
    (void) coord;
    using std::log;
    using glm::clamp;
    using glm::length;
    const float v = clamp(0.75f + 0.125f * log(4.0f * length(de)), 0.0f, 1.0f);
    return vec3(v);
  }
  virtual std::string frag() const
  {
    return
      "vec3 colour(uvec2 n, vec2 coord, vec2 de)\n"
      "{\n"
      "  float v = clamp(0.75 + 0.125 * log(4.0 * length(de)), 0.0, 1.0);\n"
      "  return vec3(v);\n"
      "}\n"
      ;
  }
};
