// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <vector>

#include "display_cpu.h"

struct display_web : public display_cpu
{
  std::vector<unsigned char> pixels;
  GLuint texture;
#ifdef HAVE_VAO
  GLuint vao;
#endif
  GLuint vbo;
  GLuint p_display;
  GLint u_display_transform;
  GLint u_display_rgb;
  GLint u_display_rect;
  GLint u_display_subframes;
  GLint u_display_srgb;

  display_web(const colour *clr);
  virtual ~display_web();
  virtual void resize(coord_t width, coord_t height);
  virtual void accumulate(const map &out);
  virtual void draw(coord_t win_width, coord_t win_height, float x0, float y0, float x1, float y1, const mat3 &T, const int srgb_conversion);
};

bool is_webgl_1(const char *version);
