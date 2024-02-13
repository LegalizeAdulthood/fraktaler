// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2024 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "colour.h"

const char *version = "#version 300 es\n";
#include "colour.vert.h"
#include "colour.frag.h"
#include "colour_default.frag.h"

#include "render.h"
#include "glutil.h"

struct colour
{
  // the bigger picture
  int image_width, image_height;
  float zoom_log_2, time;
  // current size
  int width, height;
  // textures
  GLuint t_N0, t_N1, t_NF, t_T, t_DEX, t_DEY, t_BLA, t_PTB, t_RGBA;
  // shader
  GLuint program;
  GLint u_N0, u_N1, u_NF, u_T, u_DEX, u_DEY, u_BLA, u_PTB;
  GLint u_have_N0, u_have_N1, u_have_NF, u_have_T, u_have_DEX, u_have_DEY, u_have_BLA, u_have_PTB;
  GLint u_tile, u_tile_size, u_image_size, u_zoom_log_2, u_time;
  // drawing
  GLuint fbo, vbo, vao;
  // buffer
  float *RGBA;
};

colour *colour_new()
{
  colour *c = new colour();
  c->width = -1;
  c->height = -1;
  GLuint t[9];
  glGenTextures(9, &t[0]);
  c->t_N0 = t[0];
  c->t_N1 = t[1];
  c->t_NF = t[2];
  c->t_T  = t[3];
  c->t_DEX = t[4];
  c->t_DEY = t[5];
  c->t_BLA = t[6];
  c->t_PTB = t[7];
  c->t_RGBA = t[8];
  c->program = vertex_fragment_shader(version, src_colour_vert_glsl, src_colour_frag_glsl, src_colour_default_frag_glsl);
  glUseProgram(c->program);
  c->u_N0 = glGetUniformLocation(c->program, "Internal_N0");
  c->u_N1 = glGetUniformLocation(c->program, "Internal_N1");
  c->u_NF = glGetUniformLocation(c->program, "Internal_NF");
  c->u_T  = glGetUniformLocation(c->program, "Internal_T");
  c->u_DEX = glGetUniformLocation(c->program, "Internal_DEX");
  c->u_DEY = glGetUniformLocation(c->program, "Internal_DEY");
  c->u_BLA = glGetUniformLocation(c->program, "Internal_BLA");
  c->u_PTB = glGetUniformLocation(c->program, "Internal_PTB");
  c->u_have_N0 = glGetUniformLocation(c->program, "Internal_have_N0");
  c->u_have_N1 = glGetUniformLocation(c->program, "Internal_have_N1");
  c->u_have_NF = glGetUniformLocation(c->program, "Internal_have_NF");
  c->u_have_T  = glGetUniformLocation(c->program, "Internal_have_T");
  c->u_have_DEX = glGetUniformLocation(c->program, "Internal_have_DEX");
  c->u_have_DEY = glGetUniformLocation(c->program, "Internal_have_DEY");
  c->u_have_BLA = glGetUniformLocation(c->program, "Internal_have_BLA");
  c->u_have_PTB = glGetUniformLocation(c->program, "Internal_have_PTB");
  c->u_tile = glGetUniformLocation(c->program, "Internal_tile");
  c->u_tile_size = glGetUniformLocation(c->program, "Internal_tile_size");
  c->u_image_size = glGetUniformLocation(c->program, "Internal_image_size");
  c->u_zoom_log_2 = glGetUniformLocation(c->program, "Internal_zoom_log_2");
  c->u_time = glGetUniformLocation(c->program, "Internal_time");
  glUseProgram(0);
  glGenFramebuffers(1, &c->fbo);
  glGenBuffers(1, &c->vbo);
  glBindBuffer(GL_ARRAY_BUFFER, c->vbo);
  float data[] = { 0, 0, 0, 1, 1, 0, 1, 1 };
  glBufferData(GL_ARRAY_BUFFER, sizeof(data), data, GL_STATIC_DRAW);
  glGenVertexArrays(1, &c->vao);
  glBindVertexArray(c->vao);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(0);
  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  return c;
}

void colour_set_image_size(colour *c, int width, int height)
{
  c->image_width = width;
  c->image_height = height;
}

void colour_set_zoom_log_2(colour *c, float zoom_log_2)
{
  c->zoom_log_2 = zoom_log_2;
}

void colour_set_time(colour *c, float time)
{
  c->time = time;
}

void colour_resize(colour *c, int width, int height)
{
  c->width = width;
  c->height = height;
  glActiveTexture(GL_TEXTURE + 0);
  glBindTexture(GL_TEXTURE_2D, c->t_N0);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, width, height, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, nullptr);
#define TP \
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); \
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); \
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE); \
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  TP
  glBindTexture(GL_TEXTURE_2D, c->t_N1);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, width, height, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, nullptr);
  TP
  glBindTexture(GL_TEXTURE_2D, c->t_NF);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, width, height, 0, GL_RED, GL_FLOAT, nullptr);
  TP
  glBindTexture(GL_TEXTURE_2D, c->t_T);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, width, height, 0, GL_RED, GL_FLOAT, nullptr);
  TP
  glBindTexture(GL_TEXTURE_2D, c->t_DEX);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, width, height, 0, GL_RED, GL_FLOAT, nullptr);
  TP
  glBindTexture(GL_TEXTURE_2D, c->t_DEY);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, width, height, 0, GL_RED, GL_FLOAT, nullptr);
  TP
  glBindTexture(GL_TEXTURE_2D, c->t_BLA);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, width, height, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, nullptr);
  TP
  glBindTexture(GL_TEXTURE_2D, c->t_PTB);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, width, height, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, nullptr);
  TP
  glBindTexture(GL_TEXTURE_2D, c->t_RGBA);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, nullptr);
  TP
#undef TP
  glBindTexture(GL_TEXTURE_2D, 0);
  delete[] c->RGBA;
  c->RGBA = new float[width * height * 4];
}

void colour_tile(colour *c, int x, int y, int subframe, tile *data)
{
  if (! data->RGB)
  {
    return;
  }
  if (data->width != c->width || data->height != c->height)
  {
    colour_resize(c, data->width, data->height);
  }

  glUseProgram(c->program);

  glUniform1i(c->u_N0, 0);
  glUniform1i(c->u_have_N0, !!data->N0);
  if (data->N0)
  {
    glActiveTexture(GL_TEXTURE0 + 0);
    glBindTexture(GL_TEXTURE_2D, c->t_N0);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, data->width, data->height, GL_RED_INTEGER, GL_UNSIGNED_INT, data->N0);
  }

  glUniform1i(c->u_N1, 1);
  glUniform1i(c->u_have_N1, !!data->N1);
  if (data->N1)
  {
    glActiveTexture(GL_TEXTURE0 + 1);
    glBindTexture(GL_TEXTURE_2D, c->t_N1);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, data->width, data->height, GL_RED_INTEGER, GL_UNSIGNED_INT, data->N1);
  }

  glUniform1i(c->u_NF, 2);
  glUniform1i(c->u_have_NF, !!data->NF);
  if (data->NF)
  {
    glActiveTexture(GL_TEXTURE0 + 2);
    glBindTexture(GL_TEXTURE_2D, c->t_NF);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, data->width, data->height, GL_RED, GL_FLOAT, data->NF);
  }

  glUniform1i(c->u_T, 3);
  glUniform1i(c->u_have_T, !!data->T);
  if (data->T)
  {
    glActiveTexture(GL_TEXTURE0 + 3);
    glBindTexture(GL_TEXTURE_2D, c->t_T);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, data->width, data->height, GL_RED, GL_FLOAT, data->T);
  }

  glUniform1i(c->u_DEX, 4);
  glUniform1i(c->u_have_DEX, !!data->DEX);
  if (data->DEX)
  {
    glActiveTexture(GL_TEXTURE0 + 4);
    glBindTexture(GL_TEXTURE_2D, c->t_DEX);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, data->width, data->height, GL_RED, GL_FLOAT, data->DEX);
  }

  glUniform1i(c->u_DEY, 5);
  glUniform1i(c->u_have_DEY, !!data->DEY);
  if (data->DEY)
  {
    glActiveTexture(GL_TEXTURE0 + 5);
    glBindTexture(GL_TEXTURE_2D, c->t_DEY);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, data->width, data->height, GL_RED, GL_FLOAT, data->DEY);
  }

  glUniform1i(c->u_BLA, 6);
  glUniform1i(c->u_have_BLA, !!data->BLA);
  if (data->BLA)
  {
    glActiveTexture(GL_TEXTURE0 + 6);
    glBindTexture(GL_TEXTURE_2D, c->t_BLA);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, data->width, data->height, GL_RED_INTEGER, GL_UNSIGNED_INT, data->BLA);
  }

  glUniform1i(c->u_PTB, 7);
  glUniform1i(c->u_have_PTB, !!data->PTB);
  if (data->PTB)
  {
    glActiveTexture(GL_TEXTURE0 + 7);
    glBindTexture(GL_TEXTURE_2D, c->t_PTB);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, data->width, data->height, GL_RED_INTEGER, GL_UNSIGNED_INT, data->PTB);
  }

  glUniform2i(c->u_image_size, c->image_width, c->image_height);
  glUniform1f(c->u_zoom_log_2, c->zoom_log_2);
  glUniform1f(c->time, c->time);
  glUniform2i(c->u_tile_size, data->width, data->height);
  glUniform3i(c->u_tile, x, y, subframe);

  glBindFramebuffer(GL_FRAMEBUFFER, c->fbo);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, c->t_RGBA, 0);
  GLenum buffer = GL_COLOR_ATTACHMENT0;
  glDrawBuffers(1, &buffer);
  glViewport(0, 0, data->width, data->height);

  glBindVertexArray(c->vao);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glBindVertexArray(0);

  glReadBuffer(buffer);
  glReadPixels(0, 0, data->width, data->height, GL_RGBA, GL_FLOAT, c->RGBA);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  glUseProgram(0);

  for (int y = 0; y < data->height; ++y)
  for (int x = 0; x < data->width; ++x)
  for (int z = 0; z < 3; ++z)
  {
    data->RGB[(y * data->width + x) * 3 + z] =
      c->RGBA[(y * data->width + x) * 4 + z];
  }
}
