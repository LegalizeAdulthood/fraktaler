// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2024 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "types.h"

struct image_raw;
struct image_rgb;
struct ppostprocessing;

struct display
{
  coord_t width;
  coord_t height;
  display()
  : width(0)
  , height(0)
  {
  }
  virtual ~display()
  {
  }
  virtual void resize(coord_t widthx, coord_t heightx)
  {
    width = widthx;
    height = heightx;
  }
  virtual void plot(const image_rgb &img, const ppostprocessing &post) = 0;
  virtual void plot(const image_raw &img, const ppostprocessing &post) = 0;
  virtual void draw(coord_t win_width, coord_t win_height, const mat3 &T, const int srgb_conversion = 0, bool capture = false) = 0;
  virtual void draw_rectangle(coord_t win_width, coord_t win_height, float x0, float y0, float x1, float y1, const int srgb_conversion = 0) = 0;
  virtual void draw_circles(coord_t win_width, coord_t win_height, const std::vector<glm::vec4> &circles, const int srgb_conversion = 0) = 0;

};
