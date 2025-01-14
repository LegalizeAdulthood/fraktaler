// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#ifdef HAVE_GUI

#define GLM_ENABLE_EXPERIMENTAL 1
#include <glm/gtx/matrix_transform_2d.hpp>

#include "display_gles.h"
#include "glutil.h"
#include "parallel.h"
#include "types.h"

inline float linear_to_srgb(float c) noexcept
{
  c = glm::clamp(c, 0.0f, 1.0f);
  if (c <= 0.0031308f)
  {
    return 12.92f * c;
  }
  else
  {
    return 1.055f * std::pow(c, 1.0f / 2.4f) - 0.055f;
  }
}

inline float srgb_to_linear(float c) noexcept
{
  c = glm::clamp(c, 0.0f, 1.0f);
  if (c <= 0.04045f)
  {
    return c / 12.92f;
  }
  else
  {
    return std::pow((c + 0.055f) / 1.055f, 2.4f);
  }
}

static const char *version = "#version 100\n";

static const char *vert =
  "attribute vec2 v_position;\n"
  "attribute vec2 v_texcoord;\n"
  "varying vec2 Internal_texcoord;\n"
  "uniform mat3 Internal_transform;\n"
  "void main(void)\n"
  "{\n"
  "  vec3 p = Internal_transform * vec3(v_position, 1.0);\n"
  "  gl_Position = vec4(p.xy / p.z, 0.0, 1.0);\n"
  "  Internal_texcoord = v_texcoord;\n"
  "}\n"
  ;

static const char *vert_simple =
  "attribute vec2 v_position;\n"
  "attribute vec2 v_texcoord;\n"
  "varying vec2 Internal_texcoord;\n"
  "void main(void)\n"
  "{\n"
  "  vec3 p = vec3(v_position, 1.0);\n"
  "  gl_Position = vec4(p.xy / p.z, 0.0, 1.0);\n"
  "  Internal_texcoord = v_texcoord;\n"
  "}\n"
  ;

static const char *frag_display =
  "precision highp float;\n"
  "uniform sampler2D Internal_RGB;\n"
  "uniform int Internal_subframes;\n"
  "uniform int Internal_srgb;\n"
  "varying vec2 Internal_texcoord;\n"
  "float srgb_to_linear(float c)\n"
  "{\n"
  "  c = clamp(c, 0.0, 1.0);\n"
  "  if (c <= 0.04045)\n"
  "  {\n"
  "    return c / 12.92;\n"
  "  }\n"
  "  else\n"
  "  {\n"
  "    return pow((c + 0.055) / 1.055, 2.4);\n"
  "  }\n"
  "}\n"
  "float linear_to_srgb(float c)\n"
  "{\n"
  "  c = clamp(c, 0.0, 1.0);\n"
  "  if (c <= 0.0031308)\n"
  "  {\n"
  "    return 12.92 * c;\n"
  "  }\n"
  "  else\n"
  "  {\n"
  "    return 1.055 * pow(c, 1.0 / 2.4) - 0.055;\n"
  "  }\n"
  "}\n"
  "vec3 srgb_to_linear(vec3 c)\n"
  "{\n"
  "  return vec3(srgb_to_linear(c.r), srgb_to_linear(c.g), srgb_to_linear(c.b));\n"
  "}\n"
  "vec3 linear_to_srgb(vec3 c)\n"
  "{\n"
  "  return vec3(linear_to_srgb(c.r), linear_to_srgb(c.g), linear_to_srgb(c.b));\n"
  "}\n"
  "void main(void)\n"
  "{\n"
  "  vec2 t = Internal_texcoord;\n"
  "  vec4 c = texture2D(Internal_RGB, vec2(t.x, 1.0 - t.y));\n"
  "  if (Internal_subframes == 0)\n"
  "  {\n"
  "    c = vec4(vec3(0.5), 1.0);\n"
  "  }\n"
  "  if (Internal_srgb > 0)\n"
  "  {\n"
  "    c.rgb = linear_to_srgb(c.rgb);\n"
  "  }\n"
  "  if (Internal_srgb < 0)\n"
  "  {\n"
  "    c.rgb = srgb_to_linear(c.rgb);\n"
  "  }\n"
  "  gl_FragColor = c;\n"
  "}\n"
  ;

static const char *frag_display_rectangle =
  "#extension GL_OES_standard_derivatives : require\n"
  "precision highp float;\n"
  "uniform vec4 Internal_rectangle;\n"
  "varying vec2 Internal_texcoord;\n"
  "bool in_rectangle(vec2 p, vec4 r)\n"
  "{\n"
  "  return r.x < p.x && p.x < r.z && r.y < p.y && p.y < r.w;\n"
  "}\n"
  "void main(void)\n"
  "{\n"
  "  vec2 t = Internal_texcoord;\n"
  "  vec4 d = vec4(-dFdx(t.x), -dFdy(t.y), dFdx(t.x), dFdy(t.y));\n"
  "  vec4 c = vec4(1.0, 0.8, 0.5, 0.0);\n"
  "  if (in_rectangle(t, Internal_rectangle + d))\n"
  "  {\n"
  "    if (in_rectangle(t, Internal_rectangle - d))\n"
  "    {\n"
  "      c.a = 0.25;\n"
  "    }\n"
  "    else\n"
  "    {\n"
  "      c.a = 0.75;\n"
  "    }\n"
  "  }\n"
  "  else\n"
  "  {\n"
  "    discard;\n"
  "  }\n"
  "  gl_FragColor = c;\n"
  "}\n"
  ;

static const char *frag_display_circles =
  "precision highp float;\n"
  "uniform vec4 Internal_circles[16];\n"
  "uniform int Internal_ncircles;\n"
  "varying vec2 Internal_texcoord;\n"
  "bool in_circle(vec2 p, vec4 c)\n"
  "{\n"
  "  float x = ((c.x + 1.0) / 2.0 - p.x) / c.z;\n"
  "  float y = ((c.y + 1.0) / 2.0 - p.y) / c.w;\n"
  "  return x * x + y * y < 1.0;\n"
  "}\n"
  "void main(void)\n"
  "{\n"
  "  vec2 t = Internal_texcoord;\n"
  "  vec4 c = vec4(1.0, 0.8, 0.5, 0.0);\n"
  "  for (int circle = 0; circle < 16; ++circle)\n"
  "  {\n"
  "    if (circle < Internal_ncircles)\n"
  "    {\n"
  "      if (in_circle(t, Internal_circles[circle]))\n"
  "      {\n"
  "        c.a = 0.5;\n"
  "      }\n"
  "    }\n"
  "  }\n"
  "  if (c.a == 0.0)\n"
  "  {\n"
  "    discard;\n"
  "  }\n"
  "  gl_FragColor = c;\n"
  "}\n"
  ;

display_gles::display_gles()
: display()
, pixels(0)
, have_data(false)
, texture(0)
#ifdef HAVE_VAO
, vao(0)
#endif
, vbo(0)
, p_display(0)
, u_display_rgb(0)
, u_display_rect(0)
, format(GL_RGBA)
{
  while (glGetError())
  {
  }
  glGenTextures(1, &texture);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  const GLfloat data[] = { -1, -1, 0, 0,  -1, 1, 0, 1,  1, -1, 1, 0,  1, 1, 1, 1 };
  glBufferData(GL_ARRAY_BUFFER, sizeof(data), data, GL_STATIC_DRAW);
#ifdef HAVE_VAO
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), ((char*)0) + 0 * sizeof(GLfloat)); // vertex
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), ((char*)0) + 2 * sizeof(GLfloat)); // texcoord
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glBindVertexArray(0);
#endif
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  p_display = vertex_fragment_shader(version, vert, frag_display);
  p_display = vertex_fragment_shader(version, vert, frag_display);
  p_display_rectangle = vertex_fragment_shader(version, vert_simple, frag_display_rectangle);
  p_display_circles = vertex_fragment_shader(version, vert_simple, frag_display_circles);
  glUseProgram(p_display);
  u_display_transform = glGetUniformLocation(p_display, "Internal_transform");
  u_display_rgb = glGetUniformLocation(p_display, "Internal_RGB");
  u_display_subframes = glGetUniformLocation(p_display, "Internal_subframes");
  u_display_srgb = glGetUniformLocation(p_display, "Internal_srgb");
  glUniform1i(u_display_rgb, 0);
  glUseProgram(p_display_rectangle);
  u_display_rect = glGetUniformLocation(p_display_rectangle, "Internal_rectangle");
  glUseProgram(p_display_circles);
  u_display_circles = glGetUniformLocation(p_display_circles, "Internal_circles");
  u_display_ncircles = glGetUniformLocation(p_display_circles, "Internal_ncircles");
  glUseProgram(0);
  glUseProgram(0);
  format = GL_RGBA;
#ifdef __EMSCRIPTEN__
  if (is_webgl_1((const char *) glGetString(GL_VERSION)))
  {
    format = GL_SRGB_ALPHA_EXT;
  }
#endif
  int e;
  while ((e = glGetError()))
  {
    std::fprintf(stderr, "GL ERROR %d display_gles\n", e);
  }
}

display_gles::~display_gles()
{
  while (glGetError())
  {
  }
  glDeleteProgram(p_display);
#ifdef HAVE_VAO
  glDeleteVertexArrays(1, &vao);
#endif
  glDeleteBuffers(1, &vbo);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
  glDeleteTextures(1, &texture);
  int e;
  while ((e = glGetError()))
  {
    std::fprintf(stderr, "GL ERROR %d ~display_gles\n", e);
  }
}

void display_gles::resize(coord_t width, coord_t height)
{
  while (glGetError())
  {
  }
  display::resize(width, height);
  pixels.resize(4 * width * height);
  have_data = false;
  glActiveTexture(GL_TEXTURE0);
  glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, nullptr);
  int e;
  while ((e = glGetError()))
  {
    std::fprintf(stderr, "GL ERROR %d resize\n", e);
  }
}

void display_gles::plot(image_rgb &out)
{
  while (glGetError())
  {
  }
#ifndef __EMSCRIPTEN__
  volatile bool running = true;
  parallel2d(std::thread::hardware_concurrency(), 0, width, 32, 0, height, 32, &running, [&](coord_t i, coord_t j) -> void
#else
  for (coord_t j = 0; j < height; ++j)
  for (coord_t i = 0; i < width; ++i)
#endif
  {
    coord_t k = 4 * (j * width + i);
    float A = out.RGBA[k + 3];
    if (A == 0)
    {
      for (coord_t c = 0; c < 3; ++c)
      {
        pixels[4 * ((height - 1 - j) * width + i) + c] = glm::clamp(255.0f * linear_to_srgb(0.5f), 0.0f, 255.0f);
      }
    }
    else
    {
      have_data = true;
      for (coord_t c = 0; c < 3; ++c)
      {
        pixels[4 * ((height - 1 - j) * width + i) + c] = glm::clamp(255.0f * linear_to_srgb(out.RGBA[k + c] / A), 0.0f, 255.0f);
      }
    }
    pixels[4 * ((height - 1 - j) * width + i) + 3] = 255;
#ifndef __EMSCRIPTEN__
  });
#else
  }
#endif
  glActiveTexture(GL_TEXTURE0);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, format, GL_UNSIGNED_BYTE, &pixels[0]);
  int e;
  while ((e = glGetError()))
  {
    std::fprintf(stderr, "GL ERROR %d plot rgb\n", e);
  }
}

void display_gles::plot(image_raw &out)
{
  while (glGetError())
  {
  }
#ifndef __EMSCRIPTEN__
  volatile bool running = true;
  parallel2d(std::thread::hardware_concurrency(), 0, width, 32, 0, height, 32, &running, [&](coord_t i, coord_t j) -> void
#else
  for (coord_t j = 0; j < height; ++j)
  for (coord_t i = 0; i < width; ++i)
#endif
  {
    coord_t k = j * width + i;
    coord_t w = 4 * ((height - 1 - j) * width + i);
    pixels[w + 0] = glm::clamp(255.0f * linear_to_srgb(out.R ? out.R[k] : 0.5f), 0.0f, 255.0f);
    pixels[w + 1] = glm::clamp(255.0f * linear_to_srgb(out.G ? out.G[k] : 0.5f), 0.0f, 255.0f);
    pixels[w + 2] = glm::clamp(255.0f * linear_to_srgb(out.B ? out.B[k] : 0.5f), 0.0f, 255.0f);
    pixels[w + 3] = 255;
#ifndef __EMSCRIPTEN__
  });
#else
  }
#endif
  glActiveTexture(GL_TEXTURE0);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, format, GL_UNSIGNED_BYTE, &pixels[0]);
  int e;
  while ((e = glGetError()))
  {
    std::fprintf(stderr, "GL ERROR %d plot raw\n", e);
  }
}

void display_gles::draw(coord_t win_width, coord_t win_height, const mat3 &T, const int srgb_conversion)
{
  while (glGetError())
  {
  }
  glViewport(0, 0, win_width, win_height);
  glClearColor(0.5, 0.5, 0.5, 1);
  glClear(GL_COLOR_BUFFER_BIT);
#ifdef HAVE_VAO
  glBindVertexArray(vao);
#else
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), ((char*)0) + 0 * sizeof(GLfloat)); // vertex
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), ((char*)0) + 2 * sizeof(GLfloat)); // texcoord
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
#endif
  glUseProgram(p_display);
  mat3 S = mat3(1.0f);
  // [0..w] x [0..h]
  S = glm::scale(S, vec2(float(win_width), float(win_height)));
  S = glm::scale(S, vec2(0.5f, 0.5f));
  S = glm::translate(S, vec2(1.0f));
  // [-1..1] x [-1..1]
  S = glm::inverse(S) * T * S;
  glUniformMatrix3fv(u_display_transform, 1, false, &S[0][0]);
  glUniform1i(u_display_subframes, have_data);
  glUniform1i(u_display_srgb, srgb_conversion);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glUseProgram(0);
#ifdef HAVE_VAO
  glBindVertexArray(0);
#else
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
#endif
  int e;
  while ((e = glGetError()))
  {
    std::fprintf(stderr, "GL ERROR %d draw\n", e);
  }
}

void display_gles::draw_rectangle(coord_t win_width, coord_t win_height, float x0, float y0, float x1, float y1, const int srgb_conversion)
{
  (void) srgb_conversion;
  while (glGetError())
  {
  }
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glViewport(0, 0, win_width, win_height);
#ifdef HAVE_VAO
  glBindVertexArray(vao);
#else
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), ((char*)0) + 0 * sizeof(GLfloat)); // vertex
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), ((char*)0) + 2 * sizeof(GLfloat)); // texcoord
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
#endif
  glUseProgram(p_display_rectangle);
  glUniform4f(u_display_rect, (x0 + 1) / 2, 1 - (y1 + 1) / 2, (x1 + 1) / 2, 1 - (y0 + 1) / 2);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glUseProgram(0);
#ifdef HAVE_VAO
  glBindVertexArray(0);
#else
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
#endif
  glDisable(GL_BLEND);
  int e;
  while ((e = glGetError()))
  {
    std::fprintf(stderr, "GL ERROR %d draw_rectangle\n", e);
  }
}

void display_gles::draw_circles(coord_t win_width, coord_t win_height, const std::vector<glm::vec4> &circles, const int srgb_conversion)
{
  (void) srgb_conversion;
  while (glGetError())
  {
  }
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glViewport(0, 0, win_width, win_height);
#ifdef HAVE_VAO
  glBindVertexArray(vao);
#else
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), ((char*)0) + 0 * sizeof(GLfloat)); // vertex
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), ((char*)0) + 2 * sizeof(GLfloat)); // texcoord
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
#endif
  glUseProgram(p_display_circles);
  glUniform4fv(u_display_circles, circles.size(), &circles[0][0]);
  glUniform1i(u_display_ncircles, circles.size());
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glUseProgram(0);
#ifdef HAVE_VAO
  glBindVertexArray(0);
#else
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
#endif
  glDisable(GL_BLEND);
  int e;
  while ((e = glGetError()))
  {
    std::fprintf(stderr, "GL ERROR %d draw_circles\n", e);
  }
}

bool is_webgl_1(const char *version)
{
  const char *webgl1 = "OpenGL ES 2.0 (WebGL 1.0";
  return 0 == std::strncmp(version, webgl1, std::strlen(webgl1));
}

#endif
