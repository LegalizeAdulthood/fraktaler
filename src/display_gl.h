// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <GL/glew.h>

#include "display.h"

struct display_gl : public display
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
  bool do_clear;
  int subframes;

  GLuint fbo[2];
  GLuint vao;
  GLuint vbo;
  GLuint p_colourize;
  GLint u_colourize_backbuffer;
  GLint u_colourize_clear;
  GLuint p_display;
  GLint u_display_rgb;
  GLint u_display_rect;
  GLint u_display_subframes;

  display_gl(const colour *clr);
  virtual ~display_gl();

  virtual void resize(coord_t width, coord_t height);
  virtual void set_colour(const colour *clr);
  virtual void clear();
  virtual void accumulate(const map &out);
  virtual void get_rgb(map &out) const;

  virtual void upload_raw(const map &out);
  virtual void colourize();
  virtual void draw(coord_t win_width, coord_t win_height, float x0, float y0, float x1, float y1);
};
