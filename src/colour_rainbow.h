// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "colour.h"

struct colour_rainbow : public colour
{
  colour_rainbow() { }
  virtual ~colour_rainbow() { }
  virtual std::string name() const { return "Rainbow"; }
  virtual vec3 rgb(const count_t &n, const vec2 &coord, const vec2 &de) const noexcept
  {
    (void) n;
    using std::atan2;
    using std::floor;
    using std::log;
    using std::pow;
    using glm::clamp;
    using glm::length;
    using glm::mix;
    const float pi = M_PI;
    const float k = pow(0.5f, 0.5f + coord.y);
    const float w = 0.05f;
    const float wk = w * k;
    const bool g =
      w < coord.y && coord.y < 1.0f - w &&
      wk < coord.x && coord.x < 1.0f - wk;
    float h = atan2(de.y, de.x) / (2.0f * pi);
    h -= floor(h);
    const float s = clamp(2.0f / (1.0f + length(de)) + (g ? 0.0f : 0.5f), 0.0f, 1.0f);
    const float v = clamp(0.75f + 0.125f * log(length(de)), 0.0f, 1.0f);
    const vec3 c = mix(vec3(1.0f), cos(2.0f * pi * (h + vec3(0.0f, 1.0f, 2.0f) / 3.0f)), 0.5f);
    return mix(vec3(1.0f), c, s) * v;
  }
  virtual std::string frag() const
  {
    return
      "vec3 colour(uint n, vec2 coord, vec2 de)\n"
      "{\n"
      "  float k = pow(0.5, 0.5 + coord.y);\n"
      "  float w = 0.05;\n"
      "  float wk = w * k;\n"
      "  bool g =\n"
      "    w < coord.y && coord.y < 1.0 - w &&\n"
      "    wk < coord.x  && coord.x < 1.0 - wk;\n"
      "  float h = atan(de.y, de.x) / (2.0 * pi);\n"
      "  h -= floor(h);\n"
      "  float s = clamp(2.0 / (1.0 + length(de)) + (g ? 0.0 : 0.5), 0.0, 1.0);\n"
      "  float v = clamp(0.75 + 0.125 * log(length(de)), 0.0, 1.0);\n"
      "  vec3 c = mix(vec3(1.0), cos(2.0 * pi * (h + vec3(0.0, 1.0, 2.0) / 3.0)), 0.5);\n"
      "  return mix(vec3(1.0), c, s) * v;\n"
      "}\n"
      ;
  }
};
