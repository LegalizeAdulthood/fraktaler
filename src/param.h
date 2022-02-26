// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <filesystem>
#include <iostream>

#include <mpreal.h>

#include "complex.h"
#include "floatexp.h"
#include "matrix.h"
#include "types.h"

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
  double escape_radius = 625.0;
};

struct palgorithm
{
  bool lock_maximum_reference_iterations_to_period = false;
  bool reuse_reference = false;
  bool reuse_bilinear_approximation = false;
  std::vector<std::string> number_types =
    { "float", "double", "long double", "floatexp", "softfloat"
#ifdef HAVE_FLOAT128
    , "float128"
#endif
    };
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

struct prender
{
  std::string filename = "fraktaler-3";
  bool zoom_out_sequence = false;
  double zoom_out_factor = 2.0;
  int start_frame = 0;
  int frame_count = 0;
};

struct phybrid1
{
  bool abs_x = false;
  bool abs_y = false;
  bool neg_x = false;
  bool neg_y = false;
  int power = 2;
};

inline bool operator==(const phybrid1 &a, const phybrid1 &b)
{
  return
    a.abs_x == b.abs_x &&
    a.abs_y == b.abs_y &&
    a.neg_x == b.neg_x &&
    a.neg_y == b.neg_y &&
    a.power == b.power ;
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


struct pparam
{
  phybrid formula;
  int colour_id = 0;
  plocation location;
  preference reference;
  pbailout bailout;
  palgorithm algorithm;
  pimage image;
  ptransform transform;
  prender render;
  pnewton newton;
};

std::ostream &operator<<(std::ostream &o, const pparam &p);
std::istream &operator>>(std::istream &i, pparam &p);

struct param
{
  pparam p;
  complex<mpreal> center;
  floatexp zoom;
  complex<mpreal> reference;
  mat2<double> transform;
  std::string s_iterations, s_maximum_reference_iterations, s_maximum_perturb_iterations, s_escape_radius, s_period;
  param();
  std::string to_string() const;
  void from_string(const std::string &s);
  void load_toml(const std::string &filename);
  void save_toml(const std::string &filename) const;
};

void restring_locs(param &par);
void restring_vals(param &par);
void unstring_locs(param &par);
void unstring_vals(param &par);
void home(param &par);
void zoom(param &par, double x, double y, double g, bool fixed_click = true);
void zoom(param &par, const mat3 &T, const mat3 &T0);
complex<floatexp> get_delta_c(const param &par, double x, double y);
