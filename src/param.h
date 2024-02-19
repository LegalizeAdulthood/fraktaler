// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2024 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <iostream>
#include <map>
#include <vector>

#include <mpreal.h>

#include "complex.h"
#include "floatexp.h"
#include "matrix.h"
#include "types.h"

#include "colour_default.frag.h"

struct plocation
{
  std::string real = "0";
  std::string imag = "0";
  std::string zoom = "1";
};

struct preference
{
  std::string real = "0";
  std::string imag = "0";
  count_t period = 0;
};

struct pbailout
{
  count_t iterations = 1024;
  count_t maximum_reference_iterations = 1024;
  count_t maximum_perturb_iterations = 1024;
  count_t maximum_bla_steps = 1024;
  double escape_radius = 625.0;
  double inscape_radius = 1.0 / 1024.0;
};

struct palgorithm
{
  bool lock_maximum_reference_iterations_to_period = false;
  bool reuse_reference = false;
  bool reuse_bilinear_approximation = false;
  int bla_skip_levels = 3;
};

struct pnewton
{
  int action = 3;
  bool domain = false;
  bool absolute = false;
  float power = 0.5f;
  float factor = 4.0f;
};

struct pimage
{
  int width = 1024;
  int height = 576;
  int supersampling = 1;
  int subsampling = 1;
  int subframes = 1;
};

struct ptransform
{
  bool reflect = false;
  double rotate = 0.0;
  double stretch_angle = 0.0;
  double stretch_amount = 0.0;
  bool exponential_map = false;
};

struct ppostprocessing
{
  double brightness = 0.0;
  double contrast = 0.0;
  double gamma = 1.0;
  double exposure = 0.0;
};

struct prender
{
  std::string filename = "fraktaler-3";
  bool zoom_out_sequence = false;
  double zoom_out_factor = 2.0;
  int start_frame = 0;
  int frame_count = 0;
  int prng_seed = 0;
};

struct phybrid1
{
  bool abs_x = false;
  bool abs_y = false;
  bool neg_x = false;
  bool neg_y = false;
  int power = 2;
  std::vector<opcode> opcodes; // empty means use above, otherwise above ignored
};

inline bool operator==(const phybrid1 &a, const phybrid1 &b)
{
  return
    a.abs_x == b.abs_x &&
    a.abs_y == b.abs_y &&
    a.neg_x == b.neg_x &&
    a.neg_y == b.neg_y &&
    a.power == b.power &&
    a.opcodes == b.opcodes ;
}

struct phybrid
{
  // std::vector<hybrid_form> pre;
  std::vector<phybrid1> per = { phybrid1() };
};

inline bool operator==(const phybrid &a, const phybrid &b)
{
  return a.per == b.per;
}

struct popencl
{
  int platform = -1;
  int device = 0;
  int tile_width = 128;
  int tile_height = 128;
};

#include <unordered_map>
namespace toml
{
  template<typename C, template<typename ...> class T, template<typename ...> class A> class basic_value;
  struct discard_comments;
  using value = class toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>;
}

enum atomt { t_string, t_double, t_int64 };
struct patom
{
  atomt tag = t_string;
  std::string string_ = "";
  double double_ = 0.0;
  int64_t int64_ = 0;
  patom() { }
  patom(const std::string &x) : tag(t_string), string_(x) { }
  patom(const char *x) : tag(t_string), string_(x) { }
  patom(const double &x) : tag(t_double), double_(x) { }
  patom(const int64_t &x) : tag(t_int64), int64_(x) { }
  patom(const int &x) : patom(int64_t(x)) { }
  patom(const toml::value &x);
};

inline bool operator==(const patom &a, const patom &b)
{
  if (a.tag == t_string && b.tag == t_string)
  {
    return a.string_ == b.string_;
  }
  if (a.tag == t_double && b.tag == t_double)
  {
    return a.double_ == b.double_;
  }
  if (a.tag == t_int64 && b.tag == t_int64)
  {
    return a.int64_ == b.int64_;
  }
  return false;
}

struct pcolour
{
  std::string shader = std::string(examples_default_glsl);
  std::vector<std::map<std::string, patom>> uniforms =
    { {{"name", patom("brightness")}, {"index", patom(0)}, {"component", patom(0)}, {"type", patom("float")}, {"value", patom(0.0)}}
    , {{"name", patom("contrast"  )}, {"index", patom(0)}, {"component", patom(0)}, {"type", patom("float")}, {"value", patom(0.0)}}
    , {{"name", patom("exposure"  )}, {"index", patom(0)}, {"component", patom(0)}, {"type", patom("float")}, {"value", patom(0.0)}}
    , {{"name", patom("gamma"     )}, {"index", patom(0)}, {"component", patom(0)}, {"type", patom("float")}, {"value", patom(1.0)}}
    };
};

inline bool operator==(const pcolour &a, const pcolour &b)
{
  return
    a.shader == b.shader &&
    a.uniforms == b.uniforms ;
}

struct pparam
{
  phybrid formula;
  plocation location;
  preference reference;
  pbailout bailout;
  palgorithm algorithm;
  pimage image;
  ptransform transform;
  ppostprocessing postprocessing;
  prender render;
  pnewton newton;
  popencl opencl;
  pcolour colour;
};

std::ostream &operator<<(std::ostream &o, const pparam &p);
std::istream &operator>>(std::istream &i, pparam &p);

struct param
{
  pparam p;
  std::vector<std::vector<opcode>> opss;
  std::vector<int> degrees; // same length as opss, calculated from it
  complex<mpreal> center;
  floatexp zoom;
  complex<mpreal> reference;
  mat2<double> transform;
  std::string s_iterations, s_maximum_reference_iterations, s_maximum_perturb_iterations, s_maximum_bla_steps, s_escape_radius, s_inscape_radius, s_period, s_opss;
  param();
  std::string to_string() const;
  void from_string(const std::string &s);
  void load_toml(const std::string &filename);
  void load_exr(const std::string &filename);
  void save_toml(const std::string &filename) const;
};

void restring_locs(param &par);
void restring_vals(param &par);
void unstring_locs(param &par);
void unstring_vals(param &par);
void post_edit_formula(param &par);

void home(param &par);
void zoom(param &par, double x, double y, double g, bool fixed_click = true);
void zoom(param &par, const mat3 &T, const mat3 &T0);
complex<floatexp> get_delta_c(const param &par, double x, double y);

std::vector<std::vector<opcode>> compile_formula(const phybrid &H);
int opcodes_degree(const std::vector<opcode> &ops);
std::vector<int> opcodes_degrees(const std::vector<std::vector<opcode>> &ops);
std::string print_opcodess(const std::vector<std::vector<opcode>> &opss);
std::vector<std::vector<opcode>> parse_opcodess(const std::string &s);
std::string print_opcodes(const std::vector<opcode> &opss);
std::vector<opcode> parse_opcodes(const std::string &s);
bool validate_opcodess(const std::vector<std::vector<opcode>> &opss);
bool validate_opcodes(const std::vector<opcode> &opss);
