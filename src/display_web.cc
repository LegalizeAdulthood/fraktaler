// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <glm/gtx/matrix_transform_2d.hpp>

#include "colour.h"
#include "display_web.h"
#include "glutil.h"
#include "types.h"

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

static const char *frag_display =
  "#extension GL_OES_standard_derivatives : require\n"
  "precision highp float;\n"
  "uniform sampler2D Internal_RGB;\n"
  "uniform vec4 Internal_rectangle;\n"
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
  "bool in_rectangle(vec2 p, vec4 r)\n"
  "{\n"
  "  return r.x < p.x && p.x < r.z && r.y < p.y && p.y < r.w;\n"
  "}\n"
  "void main(void)\n"
  "{\n"
  "  vec2 t = Internal_texcoord;\n"
  "  vec4 c = texture2D(Internal_RGB, vec2(t.x, 1.0 - t.y));\n"
  "  if (Internal_subframes == 0)\n"
  "  {\n"
  "    c = vec4(vec3(0.5), 1.0);\n"
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

display_web::display_web(const colour *clr)
: display_cpu(clr)
, pixels(0)
, texture(0)
#ifdef HAVE_VAO
, vao(0)
#endif
, vbo(0)
, p_display(0)
, u_display_rgb(0)
, u_display_rect(0)
{
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
  glUseProgram(p_display);
  u_display_rgb = glGetUniformLocation(p_display, "Internal_RGB");
  u_display_transform = glGetUniformLocation(p_display, "Internal_transform");
  u_display_rect = glGetUniformLocation(p_display, "Internal_rectangle");
  u_display_subframes = glGetUniformLocation(p_display, "Internal_subframes");
  u_display_srgb = glGetUniformLocation(p_display, "Internal_srgb");
  glUniform1i(u_display_rgb, 0);
  glUseProgram(0);
}

display_web::~display_web()
{
  glDeleteProgram(p_display);
#ifdef HAVE_VAO
  glDeleteVertexArrays(1, &vao);
#endif
  glDeleteBuffers(1, &vbo);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
  glDeleteTextures(1, &texture);
}

void display_web::resize(coord_t width, coord_t height)
{
  display_cpu::resize(width, height);
  pixels.resize(4 * width * height);
  glActiveTexture(GL_TEXTURE0);
#ifdef __ANDROID__
  GLenum internal_format = GL_RGBA;
  GLenum format = GL_RGBA;
#else
  GLenum internal_format = GL_SRGB8_ALPHA8;
  GLenum format = GL_RGBA;
#endif
#ifdef __EMSCRIPTEN__
  if (is_webgl_1((const char *) glGetString(GL_VERSION)))
  {
    internal_format = GL_SRGB_ALPHA;
    format = GL_SRGB_ALPHA;
  }
#endif
  glTexImage2D(GL_TEXTURE_2D, 0, internal_format, width, height, 0, format, GL_UNSIGNED_BYTE, nullptr);
}

void display_web::accumulate(const map &out)
{
  display_cpu::accumulate(out);
  #pragma omp parallel for
  for (coord_t j = 0; j < height; ++j)
  for (coord_t i = 0; i < width; ++i)
  {
    vec3 rgb = RGB[j * width + i] / float(subframes);
    for (coord_t c = 0; c < 3; ++c)
    {
      pixels[4 * ((height - 1 - j) * width + i) + c] = glm::clamp(255.0f * linear_to_srgb(rgb[c]), 0.0f, 255.0f);
    }
    pixels[4 * ((height - 1 - j) * width + i) + 3] = 255;
  }
  glActiveTexture(GL_TEXTURE0);
  GLenum format = GL_RGBA;
#ifdef __EMSCRIPTEN__
  if (is_webgl_1((const char *) glGetString(GL_VERSION)))
  {
    format = GL_SRGB_ALPHA;
  }
#endif
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, format, GL_UNSIGNED_BYTE, &pixels[0]);
}

void display_web::draw(coord_t win_width, coord_t win_height, float x0, float y0, float x1, float y1, const mat3 &T, const int srgb_conversion)
{
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
  glUniform4f(u_display_rect, (x0 + 1) / 2, 1 - (y1 + 1) / 2, (x1 + 1) / 2, 1 - (y0 + 1) / 2);
  glUniform1i(u_display_subframes, subframes);
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
}

bool is_webgl_1(const char *version)
{
  const char *webgl1 = "OpenGL ES 2.0 (WebGL 1.0";
  return 0 == std::strncmp(version, webgl1, std::strlen(webgl1));
}
