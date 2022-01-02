// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <cstdint>

#ifdef HAVE_GLEW
#include <GL/glew.h>
#else
#ifdef __ANDROID__
#include <SDL_opengles2.h>
#endif
#endif
#include <glm/glm.hpp>
#include <mpreal.h>

typedef int64_t coord_t;

typedef float smooth_t;
typedef int64_t count_t;

typedef float progress_t;

typedef float mantissa;
typedef int exponent;
struct floatexp;

struct softfloat;

template <typename real> struct complex;

template <int D, typename T> struct dual;

template <typename real> struct mat2;

typedef int channel_bit_t;
typedef int channel_mask_t;
struct map;

struct colour;
struct stats;
struct param;

template <typename real> struct blaC;
template <typename real> struct blasC;

template <typename real> struct blaR2;
template <typename real> struct blasR2;

struct reference;
struct formula;
struct formulaCbase;
struct formulaR2base;

enum number_type
{ nt_none = 0
, nt_float = 1
, nt_double = 2
, nt_longdouble = 3
, nt_floatexp = 4
, nt_softfloat = 5
#ifdef HAVE_FLOAT128
, nt_float128 = 6
#endif
};

using mpreal = mpfr::mpreal;

using vec2 = glm::vec2;
using vec3 = glm::vec3;
using mat3 = glm::mat3;

#ifdef __ANDROID__
#define CONSTEXPR
#define CONSTEXPR_STATIC const
#else
#ifndef CONSTEXPR
#define CONSTEXPR constexpr
#endif
#define CONSTEXPR_STATIC constexpr
#endif

template<typename T> T convert(const mpreal &m) noexcept
{
  return T(m);
}

#ifndef M_PI
#define M_PI 3.141592653589793
#endif
