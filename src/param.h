// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <filesystem>
#include <iostream>

#include <mpreal.h>
#include <toml.hpp>

#include "complex.h"
#include "floatexp.h"
#include "matrix.h"
#include "types.h"

struct plocation
{
  std::string real;
  std::string imag;
  std::string zoom;
  count_t period;
};
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(plocation, real, imag, zoom, period);

struct pbailout
{
  count_t iterations;
  count_t maximum_reference_iterations;
  count_t maximum_perturb_iterations;
  double escape_radius;
};
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(pbailout, iterations, maximum_reference_iterations, maximum_perturb_iterations, escape_radius);

struct palgorithm
{
  bool lock_maximum_reference_iterations_to_period;
  bool reuse_reference;
  bool reuse_bilinear_approximation;
  std::vector<std::string> number_types;
};
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(palgorithm, lock_maximum_reference_iterations_to_period, reuse_reference, reuse_bilinear_approximation, number_types);

struct pimage
{
  int width;
  int height;
  int subsampling;
  int subframes;
};
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(pimage, width, height, subsampling, subframes);

struct ptransform
{
  bool reflect;
  double rotate;
  double stretch_angle;
  double stretch_amount;
  bool exponential_map;
};
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(ptransform, reflect, rotate, stretch_angle, stretch_amount);

struct prender
{
  std::string filename;
  bool zoom_out_sequence;
  int start_frame;
  int frame_count;
};
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(prender, filename, zoom_out_sequence, start_frame, frame_count);

struct pparam
{
  plocation location;
  pbailout bailout;
  palgorithm algorithm;
  pimage image;
  ptransform transform;
  prender render;
};
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(pparam, location, bailout, algorithm, image, transform, render);

struct param
{
  pparam p;
  complex<mpreal> center;
  floatexp zoom;
  complex<mpreal> reference;
  mat2<double> transform;
  std::string s_iterations, s_maximum_reference_iterations, s_maximum_perturb_iterations, s_escape_radius;
  param();
  void parse(std::istream &s, const std::string &filename = "-");
};

void restring(param &par);
void unstring(param &par);
void home(param &par);
void zoom(param &par, double x, double y, double g, bool fixed_click = true);
void zoom(param &par, const mat3 &T, const mat3 &T0);
