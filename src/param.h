// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
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
  int start_frame = 0;
  int frame_count = 0;
};

struct pparam
{
  int formula_id = 0;
  int colour_id = 0;
  plocation location;
  preference reference;
  pbailout bailout;
  palgorithm algorithm;
  pimage image;
  ptransform transform;
  prender render;
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
