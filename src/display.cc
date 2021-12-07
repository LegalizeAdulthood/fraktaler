// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <GL/glew.h>

#include "display.h"
#include "map.h"

//static ImGuiTextBuffer shader_log;

static bool debug_program(GLuint program) {
  GLint status = 0;
  glGetProgramiv(program, GL_LINK_STATUS, &status);
  GLint length = 0;
  glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);
  char *info = nullptr;
  if (length) {
    info = new char[length + 1];
    info[0] = 0;
    glGetProgramInfoLog(program, length, 0, info);
    info[length] = 0;
  }
  if ((info && info[0]) || ! status) {
//    shader_log.appendf("\nlink info:\n%s", info ? info : "(no info log)\n");
std::cerr << "\nlink info:\n" << (info ? info : "(no info log)\n");
  }
  delete[] info;
  return status;
}

static bool debug_shader(GLuint shader, GLenum type) {
  GLint status = 0;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
  GLint length = 0;
  glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
  char *info = nullptr;
  if (length) {
    info = new char[length + 1];
    info[0] = 0;
    glGetShaderInfoLog(shader, length, 0, info);
    info[length] = 0;
  }
  if ((info && info[0]) || ! status) {
    const char *type_str = "unknown";
    switch (type) {
      case GL_VERTEX_SHADER: type_str = "vertex"; break;
      case GL_FRAGMENT_SHADER: type_str = "fragment"; break;
    }
//    shader_log.appendf("\n%s info:\n%s", type_str, info ? info : "(no info log)\n");
std::cerr << "\n" << type_str << " info:\n" << (info ? info : "(no info log)\n");
  }
  delete[] info;
  return status;
}

static GLuint vertex_fragment_shader(const char *version, const char *vert, const char *frag, const char *frag2 = nullptr)
{
//  shader_log.clear();
  bool ok = true;
  GLuint program = glCreateProgram();
  {
    GLuint shader = glCreateShader(GL_VERTEX_SHADER);
    const char *sources[] = { version, vert };
    glShaderSource(shader, 2, sources, 0);
    glCompileShader(shader);
    ok &= debug_shader(shader, GL_VERTEX_SHADER);
    glAttachShader(program, shader);
    glDeleteShader(shader);
  }
  {
    GLuint shader = glCreateShader(GL_FRAGMENT_SHADER);
    const char *sources[] = { version, frag, frag2 };
    glShaderSource(shader, frag2 ? 3 : 2, sources, 0);
    glCompileShader(shader);
    ok &= debug_shader(shader, GL_FRAGMENT_SHADER);
    glAttachShader(program, shader);
    glDeleteShader(shader);
  }
  glLinkProgram(program);
  ok &= debug_program(program);
  if (! ok)
  {
    glDeleteProgram(program);
    program = 0;
  }
  return program;
}

display::display()
: tex_width(0)
, tex_height(0)
{
  glGenTextures(TEXTURES, &texture[0]);
  for (int t = 0; t < TEXTURES; ++t)
  {
    glActiveTexture(GL_TEXTURE0 + t);
    glBindTexture(GL_TEXTURE_2D, texture[t]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, t == TEXTURE_RGB ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, t == TEXTURE_RGB ? GL_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  }
  glGenFramebuffers(1, &fbo);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
  glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture[TEXTURE_RGB], 0);
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
  const char *version =
    "#version 330 core\n"
    ;
  const char *vert =
    "layout (location = 0) in vec2 v_position;\n"
    "layout (location = 1) in vec2 v_texcoord;\n"
    "out vec2 Internal_texcoord;\n"
    "void main(void)\n"
    "{\n"
    "  gl_Position = vec4(v_position, 0.0, 1.0);\n"
    "  Internal_texcoord = v_texcoord;\n"
    "}\n"
    ;
  const char *frag_display =
    "uniform sampler2D Internal_RGB;\n"
    "uniform vec4 Internal_rectangle;\n"
    "in vec2 Internal_texcoord;\n"
    "out vec4 Internal_colour;\n"
    "bool in_rectangle(vec2 p, vec4 r)\n"
    "{\n"
    "  return r.x < p.x && p.x < r.z && r.y < p.y && p.y < r.w;\n"
    "}\n"
    "void main(void)\n"
    "{\n"
    "  vec2 t = Internal_texcoord;\n"
    "  vec4 c = texture(Internal_RGB, t);\n"
    "  vec4 d = vec4(-dFdx(t.x), -dFdy(t.y), dFdx(t.x), dFdy(t.y));\n"
    "  if (in_rectangle(Internal_texcoord, Internal_rectangle + d))\n"
    "  {\n"
    "    if (in_rectangle(Internal_texcoord, Internal_rectangle - d))\n"
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
  const char *frag_colourize =
    "uniform sampler2D Internal_DEX;\n"
    "uniform sampler2D Internal_DEY;\n"
    "uniform sampler2D Internal_T;\n"
    "uniform sampler2D Internal_NF;\n"
    "in vec2 Internal_texcoord;\n"
    "out vec4 Internal_colour;\n"
    "const float pi = 3.141592653;\n"
    "vec3 colour(void);\n"
    "vec2 getDE(void)\n"
    "{\n"
    "  return vec2(texture(Internal_DEX, Internal_texcoord).x, texture(Internal_DEY, Internal_texcoord).x);\n"
    "}\n"
    "float getT(void)\n"
    "{\n"
    "  return texture(Internal_T, Internal_texcoord).x;\n"
    "}\n"
    "float getNF(void)\n"
    "{\n"
    "  return texture(Internal_NF, Internal_texcoord).x;\n"
    "}\n"
    "void main(void)\n"
    "{\n"
    "  Internal_colour = vec4(colour(), 1.0);\n"
    "}\n"
    ;
  const char *frag_colourize_user =
    "vec3 colour(void)\n"
    "{\n"
    "  vec2 de = getDE();\n"
    "  vec2 ex = vec2(getT(), 1.0 - getNF());\n"
    "  float k = pow(0.5, 0.5 - ex.y);\n"
    "  float w = 0.05;\n"
    "  bool g = w < ex.y && ex.y < 1.0 - w &&\n"
    "    w * k < ex.x  && ex.x < 1.0 - w * k;\n"
    "  float h = atan(de.y, de.x) / (2.0 * pi);\n"
    "  h -= floor(h);\n"
    "  float s = clamp(2.0 / (1.0 + length(de)) + (g ? 0.0 : 0.5), 0.0, 1.0);\n"
    "  float v = clamp(0.75 + 0.125 * log(length(de)), 0.0, 1.0);\n"
    "  vec3 c = mix(vec3(1.0), cos(2.0 * pi * (h + vec3(0.0, 1.0, 2.0) / 3.0)), 0.5);\n"
    "  return mix(vec3(1.0), c, s) * v;\n"
    "}\n"
    ;
  p_display = vertex_fragment_shader(version, vert, frag_display);
  glUseProgram(p_display);
  glUniform1i(glGetUniformLocation(p_display, "Internal_RGB"), TEXTURE_RGB);
  u_display_rect = glGetUniformLocation(p_display, "Internal_rectangle");
  glUseProgram(0);
  p_colourize = vertex_fragment_shader(version, vert, frag_colourize, frag_colourize_user);
  glUseProgram(p_colourize);
  // FIXME TODO
  glUniform1i(glGetUniformLocation(p_colourize, "Internal_DEX"), TEXTURE_DEX);
  glUniform1i(glGetUniformLocation(p_colourize, "Internal_DEY"), TEXTURE_DEY);
  glUniform1i(glGetUniformLocation(p_colourize, "Internal_T"), TEXTURE_T);
  glUniform1i(glGetUniformLocation(p_colourize, "Internal_NF"), TEXTURE_NF);
  glUseProgram(0);
}

display::~display()
{
  glDeleteProgram(p_colourize);
  glDeleteProgram(p_display);
  glDeleteVertexArrays(1, &vao);
  glDeleteBuffers(1, &vbo);
  glDeleteFramebuffers(1, &fbo);
  for (int t = 0; t < TEXTURES; ++t)
  {
    glActiveTexture(GL_TEXTURE0 + t);
    glBindTexture(GL_TEXTURE_2D, 0);
  }
  glDeleteTextures(TEXTURES, &texture[0]);
}

void display::resize(const map &out)
{
  if (out.width == tex_width && out.height == tex_height)
  {
    return;
  }
  coord_t w = out.width;
  coord_t h = out.height;
  tex_width = w;
  tex_height = h;
  // half RGB for output
  glActiveTexture(GL_TEXTURE0 + TEXTURE_RGB);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, w, h, 0, GL_RGB, GL_HALF_FLOAT, 0);
  // uint R for N0, N1
  GLuint zeroui = 0;
  glActiveTexture(GL_TEXTURE0 + TEXTURE_N0);
  if (out.N0)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, w, h, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, 1, 1, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, &zeroui);
  }
  glActiveTexture(GL_TEXTURE0 + TEXTURE_N1);
  if (out.N1)
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
  if (out.NF)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1, 1, 0, GL_RED, GL_FLOAT, &zerof);
  }
  glActiveTexture(GL_TEXTURE0 + TEXTURE_T);
  if (out.T)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1, 1, 0, GL_RED, GL_FLOAT, &zerof);
  }
  glActiveTexture(GL_TEXTURE0 + TEXTURE_DEX);
  if (out.DEX)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1, 1, 0, GL_RED, GL_FLOAT, &zerof);
  }
  glActiveTexture(GL_TEXTURE0 + TEXTURE_DEY);
  if (out.DEY)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, 0);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1, 1, 0, GL_RED, GL_FLOAT, &zerof);
  }
}

void display::upload_raw(const map &out)
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

void display::colourize()
{
  glViewport(0, 0, tex_width, tex_height);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
  glBindVertexArray(vao);
  glUseProgram(p_colourize);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glUseProgram(0);
  glBindVertexArray(0);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
  glActiveTexture(GL_TEXTURE0 + TEXTURE_RGB);
  glGenerateMipmap(GL_TEXTURE_2D);
}

void display::download_rgb(map &out)
{
  assert(tex_width == out.width);
  assert(tex_height == out.height);
  if (out.RGB)
  {
    glActiveTexture(GL_TEXTURE0 + TEXTURE_RGB);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_HALF_FLOAT, out.RGB);
  }
}

void display::draw(coord_t win_width, coord_t win_height, float x0, float y0, float x1, float y1)
{
  glViewport(0, 0, win_width, win_height);
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
  glUniform4f(u_display_rect, (x0 + 1) / 2, 1 - (y1 + 1) / 2, (x1 + 1) / 2, 1 - (y0 + 1) / 2);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glUseProgram(0);
  glBindVertexArray(0);
}
