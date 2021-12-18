// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <glm/gtx/matrix_transform_2d.hpp>

#include "colour.h"
#include "display_gl.h"
#include "map.h"

static const char *version =
  "#version 330 core\n"
  ;

static const char *vert =
  "layout (location = 0) in vec2 v_position;\n"
  "layout (location = 1) in vec2 v_texcoord;\n"
  "out vec2 Internal_texcoord;\n"
  "uniform mat3 Internal_transform;\n"
  "void main(void)\n"
  "{\n"
  "  vec3 p = Internal_transform * vec3(v_position, 1.0);\n"
  "  gl_Position = vec4(p.xy / p.z, 0.0, 1.0);\n"
  "  Internal_texcoord = v_texcoord;\n"
  "}\n"
  ;

static const char *frag_display =
  "uniform sampler2D Internal_RGB;\n"
  "uniform vec4 Internal_rectangle;\n"
  "uniform int Internal_subframes;\n"
  "in vec2 Internal_texcoord;\n"
  "out vec4 Internal_colour;\n"
  "bool in_rectangle(vec2 p, vec4 r)\n"
  "{\n"
  "  return r.x < p.x && p.x < r.z && r.y < p.y && p.y < r.w;\n"
  "}\n"
  "void main(void)\n"
  "{\n"
  "  vec2 t = Internal_texcoord;\n"
  "  vec4 c = texture(Internal_RGB, vec2(t.x, 1.0 - t.y));\n"
  "  if (Internal_subframes == 0)\n"
  "  {\n"
  "    c = vec4(vec3(0.5), 1.0);\n"
  "  }\n"
  "  else\n"
  "  {\n"
  "    c /= float(Internal_subframes);\n"
  "  }\n"
  "  vec4 d = vec4(-dFdx(t.x), -dFdy(t.y), dFdx(t.x), dFdy(t.y));\n"
  "  if (in_rectangle(t, Internal_rectangle + d))\n"
  "  {\n"
  "    if (in_rectangle(t, Internal_rectangle - d))\n"
  "    {\n"
  "      c.rgb = mix(c.rgb, vec3(1.0, 0.8, 0.5), 0.5);\n"
  "    }\n"
  "    else\n"
  "    {\n"
  "      c.rgb = mix(c.rgb, vec3(1.0, 0.8, 0.5), 0.75);\n"
  "    }\n"
  "  }\n"
  "  Internal_colour = c;\n"
  "}\n"
  ;

static const char *frag_colourize =
  "uniform sampler2D Internal_BackBuffer;\n"
  "uniform sampler2D Internal_DEX;\n"
  "uniform sampler2D Internal_DEY;\n"
  "uniform sampler2D Internal_T;\n"
  "uniform sampler2D Internal_NF;\n"
  "uniform usampler2D Internal_N0;\n"
  "uniform usampler2D Internal_N1;\n"
  "uniform bool Internal_Clear;\n"
  "in vec2 Internal_texcoord;\n"
  "out vec4 Internal_colour;\n"
  "const float pi = 3.141592653;\n"
  "vec2 Internal_texcoord2 = vec2(Internal_texcoord.x, 1.0 - Internal_texcoord.y);\n"
  "vec3 colour(uvec2 n, vec2 coord, vec2 de);\n"
  "vec2 getDE(void)\n"
  "{\n"
  "  return vec2(texture(Internal_DEX, Internal_texcoord2).x, texture(Internal_DEY, Internal_texcoord2).x);\n"
  "}\n"
  "float getT(void)\n"
  "{\n"
  "  return texture(Internal_T, Internal_texcoord2).x;\n"
  "}\n"
  "float getNF(void)\n"
  "{\n"
  "  return texture(Internal_NF, Internal_texcoord2).x;\n"
  "}\n"
  "uvec2 getN(void)\n"
  "{\n"
  "  return uvec2(texture(Internal_N0, Internal_texcoord2).x, texture(Internal_N1, Internal_texcoord2).x);\n"
  "}\n"
  "void main(void)\n"
  "{\n"
  "  vec4 back = texture(Internal_BackBuffer, Internal_texcoord);\n"
  "  vec4 front = vec4(colour(getN(), vec2(getT(), getNF()), getDE()), 1.0);\n"
  "  Internal_colour = (Internal_Clear ? vec4(0.0) : back) + front;\n"
  "}\n"
  ;

display_gl::display_gl(const colour *clr)
: tex_width(0)
, tex_height(0)
, pingpong(0)
, do_clear(false)
, subframes(0)
, p_colourize(0)
{
  glGenTextures(TEXTURES, &texture[0]);
  for (int t = 0; t < TEXTURES; ++t)
  {
    glActiveTexture(GL_TEXTURE0 + t);
    glBindTexture(GL_TEXTURE_2D, texture[t]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, t == TEXTURE_RGB0 || t == TEXTURE_RGB1 ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  }
  glGenFramebuffers(2, &fbo[0]);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo[0]);
  glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture[TEXTURE_RGB0], 0);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo[1]);
  glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture[TEXTURE_RGB1], 0);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  const GLfloat data[] = { -1, -1, 0, 0,  -1, 1, 0, 1,  1, -1, 1, 0,  1, 1, 1, 1 };
  glBufferData(GL_ARRAY_BUFFER, sizeof(data), data, GL_STATIC_DRAW);
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), ((char*)0) + 0 * sizeof(GLfloat)); // vertex
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), ((char*)0) + 2 * sizeof(GLfloat)); // texcoord
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  p_display = vertex_fragment_shader(version, vert, frag_display);
  glUseProgram(p_display);
  u_display_transform = glGetUniformLocation(p_display, "Internal_transform");
  u_display_rgb = glGetUniformLocation(p_display, "Internal_RGB");
  u_display_rect = glGetUniformLocation(p_display, "Internal_rectangle");
  u_display_subframes = glGetUniformLocation(p_display, "Internal_subframes");
  glUseProgram(0);
  set_colour(clr);
}

void display_gl::set_colour(const colour *clr)
{
  if (p_colourize)
  {
    glDeleteProgram(p_colourize);
    p_colourize = 0;
  }
  std::string frag_colourize_user = clr->frag();
  p_colourize = vertex_fragment_shader(version, vert, frag_colourize, frag_colourize_user.c_str());
  glUseProgram(p_colourize);
  mat3 id3(1.0f);
  glUniformMatrix3fv(glGetUniformLocation(p_colourize, "Internal_transform"), 1, false, &id3[0][0]);
  u_colourize_backbuffer = glGetUniformLocation(p_colourize, "Internal_BackBuffer");
  u_colourize_clear = glGetUniformLocation(p_colourize, "Internal_Clear");
  glUniform1i(glGetUniformLocation(p_colourize, "Internal_DEX"), TEXTURE_DEX);
  glUniform1i(glGetUniformLocation(p_colourize, "Internal_DEY"), TEXTURE_DEY);
  glUniform1i(glGetUniformLocation(p_colourize, "Internal_T"), TEXTURE_T);
  glUniform1i(glGetUniformLocation(p_colourize, "Internal_NF"), TEXTURE_NF);
  glUniform1i(glGetUniformLocation(p_colourize, "Internal_N0"), TEXTURE_N0);
  glUniform1i(glGetUniformLocation(p_colourize, "Internal_N1"), TEXTURE_N1);
  glUseProgram(0);
}

display_gl::~display_gl()
{
  glDeleteProgram(p_colourize);
  glDeleteProgram(p_display);
  glDeleteVertexArrays(1, &vao);
  glDeleteBuffers(1, &vbo);
  glDeleteFramebuffers(2, &fbo[0]);
  for (int t = 0; t < TEXTURES; ++t)
  {
    glActiveTexture(GL_TEXTURE0 + t);
    glBindTexture(GL_TEXTURE_2D, 0);
  }
  glDeleteTextures(TEXTURES, &texture[0]);
}

void display_gl::resize(coord_t out_width, coord_t out_height)
{
  if (out_width == tex_width && out_height == tex_height)
  {
    return;
  }
  coord_t w = out_width;
  coord_t h = out_height;
  tex_width = w;
  tex_height = h;
  do_clear = true;
  subframes = 0;
  // float RGB for output
  glActiveTexture(GL_TEXTURE0 + TEXTURE_RGB0);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, w, h, 0, GL_RGB, GL_FLOAT, nullptr);
  glGenerateMipmap(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0 + TEXTURE_RGB1);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, w, h, 0, GL_RGB, GL_FLOAT, nullptr);
  glGenerateMipmap(GL_TEXTURE_2D);
  // uint R for N0, N1
  GLuint zeroui = 0;
  glActiveTexture(GL_TEXTURE0 + TEXTURE_N0);
  if (true) // out.N0)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, w, h, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, 1, 1, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, &zeroui);
  }
  glActiveTexture(GL_TEXTURE0 + TEXTURE_N1);
  if (true) // out.N1)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, w, h, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, 1, 1, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, &zeroui);
  }
  // float R for NF, T, DEX, DEY
  GLfloat zerof = 0;
  glActiveTexture(GL_TEXTURE0 + TEXTURE_NF);
  if (true) // out.NF)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1, 1, 0, GL_RED, GL_FLOAT, &zerof);
  }
  glActiveTexture(GL_TEXTURE0 + TEXTURE_T);
  if (true) // out.T)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1, 1, 0, GL_RED, GL_FLOAT, &zerof);
  }
  glActiveTexture(GL_TEXTURE0 + TEXTURE_DEX);
  if (true) // out.DEX)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1, 1, 0, GL_RED, GL_FLOAT, &zerof);
  }
  glActiveTexture(GL_TEXTURE0 + TEXTURE_DEY);
  if (true) // out.DEY)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1, 1, 0, GL_RED, GL_FLOAT, &zerof);
  }
}

void display_gl::upload_raw(const map &out)
{
  assert(tex_width == out.width);
  assert(tex_height == out.height);
  const coord_t w = tex_width;
  const coord_t h = tex_height;
  if (out.N0)
  {
    glActiveTexture(GL_TEXTURE0 + TEXTURE_N0);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RED_INTEGER, GL_UNSIGNED_INT, out.N0);
  }
  if (out.N1)
  {
    glActiveTexture(GL_TEXTURE0 + TEXTURE_N1);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RED_INTEGER, GL_UNSIGNED_INT, out.N1);
  }
  if (out.NF)
  {
    glActiveTexture(GL_TEXTURE0 + TEXTURE_NF);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RED, GL_FLOAT, out.NF);
  }
  if (out.T)
  {
    glActiveTexture(GL_TEXTURE0 + TEXTURE_T);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RED, GL_FLOAT, out.T);
  }
  if (out.DEX)
  {
    glActiveTexture(GL_TEXTURE0 + TEXTURE_DEX);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RED, GL_FLOAT, out.DEX);
  }
  if (out.DEY)
  {
    glActiveTexture(GL_TEXTURE0 + TEXTURE_DEY);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RED, GL_FLOAT, out.DEY);
  }
}

void display_gl::colourize()
{
  if (do_clear)
  {
    subframes = 0;
  }
  glViewport(0, 0, tex_width, tex_height);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo[pingpong]);
  glBindVertexArray(vao);
  glUseProgram(p_colourize);
  glUniform1i(u_colourize_clear, do_clear);
  glUniform1i(u_colourize_backbuffer, TEXTURE_RGB0 + ! pingpong);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glUseProgram(0);
  glBindVertexArray(0);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
  glActiveTexture(GL_TEXTURE0 + TEXTURE_RGB0 + pingpong);
  glGenerateMipmap(GL_TEXTURE_2D);
  do_clear = false;
  pingpong = ! pingpong;
  subframes++;
}

void display_gl::clear()
{
  do_clear = true;
}

void display_gl::accumulate(const map &out)
{
  upload_raw(out);
  colourize();
}

void display_gl::get_rgb(map &out) const
{
#ifndef __EMSCRIPTEN__
  assert(tex_width == out.width);
  assert(tex_height == out.height);
  if (out.RGB)
  {
    glActiveTexture(GL_TEXTURE0 + TEXTURE_RGB0 + ! pingpong);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_HALF_FLOAT, out.RGB);
    #pragma omp parallel for
    for (coord_t j = 0; j < out.height; ++j)
    for (coord_t i = 0; i < out.width; ++i)
    for (coord_t c = 0; c < 3; ++c)
    {
      out.RGB[3 * (j * out.width + i) + c] /= subframes;
    }
  }
#endif
}

void display_gl::draw(coord_t win_width, coord_t win_height, float x0, float y0, float x1, float y1, const mat3 &T)
{
  glViewport(0, 0, win_width, win_height);
  glClearColor(0.5, 0.5, 0.5, 1);
  glClear(GL_COLOR_BUFFER_BIT);
  if (tex_width * win_height >= tex_height * win_width)
  {
    // texture is wider aspect than screen
    coord_t h = tex_height * win_width / tex_width;
    coord_t y = (win_height - h) / 2;
    glViewport(0, y, win_width, h);
  }
  else
  {
    // texture is taller aspect than screen
    coord_t w = tex_width * win_height / tex_height;
    coord_t x = (win_width - w) / 2;
    glViewport(x, 0, w, win_height);
  }
  glBindVertexArray(vao);
  glUseProgram(p_display);
  mat3 S = mat3(1.0f);
  // [0..w] x [0..h]
  S = glm::scale(S, vec2(float(win_width), float(win_height)));
  S = glm::scale(S, vec2(0.5f, 0.5f));
  S = glm::translate(S, vec2(1.0f));
  // [-1..1] x [-1..1]
  S = glm::inverse(S) * T * S;
  glUniformMatrix3fv(u_display_transform, 1, false, &S[0][0]);
  glUniform1i(u_display_rgb, TEXTURE_RGB0 + ! pingpong);
  glUniform4f(u_display_rect, (x0 + 1) / 2, 1 - (y1 + 1) / 2, (x1 + 1) / 2, 1 - (y0 + 1) / 2);
  glUniform1i(u_display_subframes, subframes);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glUseProgram(0);
  glBindVertexArray(0);
}
