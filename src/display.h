// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <GL/glew.h>

#include "types.h"

struct display
{
  static constexpr int TEXTURE_RGB = 0;
  static constexpr int TEXTURE_N0  = 1;
  static constexpr int TEXTURE_N1  = 2;
  static constexpr int TEXTURE_NF  = 3;
  static constexpr int TEXTURE_T   = 4;
  static constexpr int TEXTURE_DEX = 5;
  static constexpr int TEXTURE_DEY = 6;
  static constexpr int TEXTURES    = 7;

  coord_t tex_width;
  coord_t tex_height;

  GLuint texture[TEXTURES];

  GLuint fbo;
  GLuint vao;
  GLuint vbo;
  GLuint p_colourize;
  GLuint p_display;
  GLint u_display_rect;

  display();
  ~display();

  void resize(const map &out);
  void upload_raw(const map &out);
  void colourize();
  void draw(coord_t win_width, coord_t win_height, float x0, float y0, float x1, float y1);
  void download_rgb(map &out);
};
