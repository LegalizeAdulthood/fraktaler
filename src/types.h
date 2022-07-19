// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <cstdint>

#ifdef HAVE_GLEW
#include <GL/glew.h>
#endif
#include <SDL_opengles2.h>
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
struct wlookup;

template <typename real> struct blaR2;
template <typename real> struct blasR2;

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

extern const char *nt_string[
#ifdef HAVE_FLOAT128
  7
#else
  6
#endif
];

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

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

inline void syncfs(void)
{
#ifdef __EMSCRIPTEN__
  EM_ASM(
    FS.syncfs(false, function (err) {
      /* ignore error, don't wait for done */
    });
  );
#endif
}
