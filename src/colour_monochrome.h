// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "colour.h"

struct colour_monochrome : public colour
{
  colour_monochrome() { }
  virtual ~colour_monochrome() { }
  virtual std::string name() const { return "Monochrome"; }
  virtual vec3 rgb(const count_t &n, const vec2 &coord, const vec2 &de) const noexcept
  {
    (void) n;
    using std::log;
    using glm::clamp;
    using glm::length;
    const float v = clamp(0.75f + 0.125f * log(length(16.0f * de)), 0.0f, 1.0f);
    return vec3(v);
  }
  virtual std::string frag() const
  {
    return
      "vec3 colour(uint n, vec2 coord, vec2 de)\n"
      "{\n"
      "  float v = clamp(0.75 + 0.125 * log(16.0 * length(de)), 0.0, 1.0);\n"
      "  return vec3(v);\n"
      "}\n"
      ;
  }
};
