// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2024 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

const float brightness = 0.0;
const float contrast = 0.0;
const float exposure = 0.0;
const float gamma = 2.0;

vec3 colour(void)
{
  vec2 de = getDE();
  float v = 0.75 + 0.125 * 0.25 * log(4.0 * 4.0 * dot(de, de));
  v = exp2(exposure) * pow(max(0.0, (v + brightness - 0.5) * exp2(contrast) + 0.5), 1.0 / gamma);
  return vec3(v);
}
