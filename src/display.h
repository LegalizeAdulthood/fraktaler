// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <GL/glew.h>

#include "types.h"

struct display
{
  static constexpr int TEXTURE_RGB0 = 0;
  static constexpr int TEXTURE_RGB1 = 1;
  static constexpr int TEXTURE_N0   = 2;
  static constexpr int TEXTURE_N1   = 3;
  static constexpr int TEXTURE_NF   = 4;
  static constexpr int TEXTURE_T    = 5;
  static constexpr int TEXTURE_DEX  = 6;
  static constexpr int TEXTURE_DEY  = 7;
  static constexpr int TEXTURES     = 8;

  coord_t tex_width;
  coord_t tex_height;

  GLuint texture[TEXTURES];
  int pingpong;
  bool clear;

  GLuint fbo[2];
  GLuint vao;
  GLuint vbo;
  GLuint p_colourize;
  GLint u_colourize_backbuffer;
  GLint u_colourize_clear;
  GLuint p_display;
  GLint u_display_rgb;
  GLint u_display_rect;

  display();
  ~display();

  void resize(const map &out);
  void upload_raw(const map &out);
  void colourize();
  void draw(coord_t win_width, coord_t win_height, float x0, float y0, float x1, float y1);
  void download_rgb(map &out);
};
