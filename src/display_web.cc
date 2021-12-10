// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "display_web.h"
#include "glutil.h"
#include "types.h"

static const char *version = "#version 100\n";

static const char *vert =
  "attribute vec2 v_position;\n"
  "attribute vec2 v_texcoord;\n"
  "varying vec2 Internal_texcoord;\n"
  "void main(void)\n"
  "{\n"
  "  gl_Position = vec4(v_position, 0.0, 1.0);\n"
  "  Internal_texcoord = v_texcoord;\n"
  "}\n"
  ;

static const char *frag_display =
  "#extension GL_OES_standard_derivatives : require\n"
  "precision highp float;\n"
  "uniform sampler2D Internal_RGB;\n"
  "uniform vec4 Internal_rectangle;\n"
  "varying vec2 Internal_texcoord;\n"
  "bool in_rectangle(vec2 p, vec4 r)\n"
  "{\n"
  "  return r.x < p.x && p.x < r.z && r.y < p.y && p.y < r.w;\n"
  "}\n"
  "void main(void)\n"
  "{\n"
  "  vec4 c = texture2D(Internal_RGB, vec2(Internal_texcoord.x, 1.0 - Internal_texcoord.y));\n"
  "  vec2 t = Internal_texcoord;\n"
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
  "  gl_FragColor = c;\n"
  "}\n"
  ;

display_web::display_web(const colour *clr)
: display_cpu(clr)
, pixels(0)
, texture(0)
, vao(0)
, vbo(0)
, p_display(0)
, u_display_rgb(0)
, u_display_rect(0)
{
  glGenTextures(1, &texture);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
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
  u_display_rgb = glGetUniformLocation(p_display, "Internal_RGB");
  glUniform1i(u_display_rgb, 0);
  u_display_rect = glGetUniformLocation(p_display, "Internal_rectangle");
  glUseProgram(0);
}

display_web::~display_web()
{
  glDeleteProgram(p_display);
  glDeleteVertexArrays(1, &vao);
  glDeleteBuffers(1, &vbo);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
  glDeleteTextures(1, &texture);
}

void display_web::resize(coord_t width, coord_t height)
{
  display_cpu::resize(width, height);
  pixels.resize(3 * width * height);
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
      pixels[3 * ((height - 1 - j) * width + i) + c] = glm::clamp(255.0f * rgb[c], 0.0f, 255.0f);
    }
  }
  glActiveTexture(GL_TEXTURE0);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, &pixels[0]);
}

void display_web::draw(coord_t win_width, coord_t win_height, float x0, float y0, float x1, float y1)
{
  glViewport(0, 0, win_width, win_height);
  glUseProgram(p_display);
  glBindVertexArray(vao);
  glUseProgram(p_display);
  glUniform4f(u_display_rect, (x0 + 1) / 2, 1 - (y1 + 1) / 2, (x1 + 1) / 2, 1 - (y0 + 1) / 2);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glUseProgram(0);
  glBindVertexArray(0);
}
