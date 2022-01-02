// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <cassert>
#include <filesystem>

#ifdef HAVE_EXR
#include <half.h>
#else
typedef float half;
#endif

#include "types.h"
#include "complex.h"

#ifdef RGB
#undef RGB
#endif

constexpr channel_bit_t Channel_N0  = 0;
constexpr channel_bit_t Channel_N1  = 1;
constexpr channel_bit_t Channel_NF  = 2;
constexpr channel_bit_t Channel_T   = 3;
constexpr channel_bit_t Channel_DEX = 4;
constexpr channel_bit_t Channel_DEY = 5;
constexpr channel_bit_t Channel_R   = 6;
constexpr channel_bit_t Channel_G   = 7;
constexpr channel_bit_t Channel_B   = 8;

constexpr channel_bit_t Channel_Count = 9;

constexpr channel_mask_t Channels_RGB = (1 << Channel_R) | (1 << Channel_G) | (1 << Channel_B);
constexpr channel_mask_t Channels_default = (1 << Channel_Count) - 1;

constexpr count_t Nbias_default = 1024;

struct map
{
  coord_t width;
  coord_t height;
  count_t maxiters;
  count_t Nbias;
  uint32_t *N0, *N1;
  float *NF, *T, *DEX, *DEY;
  half *RGB;

  inline map(const coord_t w, const coord_t h, const count_t maxiters, const channel_mask_t channels = Channels_default, const count_t bias = Nbias_default)
  : width(w)
  , height(h)
  , maxiters(maxiters)
  , Nbias(bias)
  , N0(nullptr)
  , N1(nullptr)
  , NF(nullptr)
  , T(nullptr)
  , DEX(nullptr)
  , DEY(nullptr)
  , RGB(nullptr)
  {
    const coord_t count = w * h;
    if (channels & (1 << Channel_N0 )) { N0  = new uint32_t [count]; }
    if (channels & (1 << Channel_N1 )) { N1  = new uint32_t [count]; }
    if (channels & (1 << Channel_NF )) { NF  = new float    [count]; }
    if (channels & (1 << Channel_T  )) { T   = new float    [count]; }
    if (channels & (1 << Channel_DEX)) { DEX = new float    [count]; }
    if (channels & (1 << Channel_DEY)) { DEY = new float    [count]; }
    if ((channels & (1 << Channel_R)) ||
        (channels & (1 << Channel_G)) ||
        (channels & (1 << Channel_B))) { RGB = new half     [count * 3]; }
  }

  inline ~map()
  {
    delete[] N0;
    delete[] N1;
    delete[] NF;
    delete[] T;
    delete[] DEX;
    delete[] DEY;
    delete[] RGB;
  }

  inline void setN(const coord_t x, const coord_t y, const count_t n)
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    uint64_t nn = n + Nbias;
    if (n >= maxiters)
    {
      nn = ~(uint64_t(0));
    }
    if (N0)
    {
      N0[y * width + x] = nn;
    }
    if (N1)
    {
      N1[y * width + x] = nn >> 32;
    }
  }

  inline void setNF(const coord_t x, const coord_t y, const float nf)
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    if (NF)
    {
      NF[y * width + x] = nf;
    }
  }

  inline void setT(const coord_t x, const coord_t y, const float t)
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    if (T)
    {
      T[y * width + x] = t;
    }
  }

  inline void setDEX(const coord_t x, const coord_t y, const float dex)
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    if (DEX)
    {
      DEX[y * width + x] = dex;
    }
  }

  inline void setDEY(const coord_t x, const coord_t y, const float dey)
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    if (DEY)
    {
      DEY[y * width + x] = dey;
    }
  }

  inline void setDE(const coord_t x, const coord_t y, const complex<float> de)
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    setDEX(x, y, de.x);
    setDEY(x, y, de.y);
  }

  inline void setR(const coord_t x, const coord_t y, const half r)
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    if (RGB)
    {
      RGB[3 * (y * width + x) + 0] = r;
    }
  }

  inline void setG(const coord_t x, const coord_t y, const half g)
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    if (RGB)
    {
      RGB[3 * (y * width + x) + 1] = g;
    }
  }

  inline void setB(const coord_t x, const coord_t y, const half b)
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    if (RGB)
    {
      RGB[3 * (y * width + x) + 2] = b;
    }
  }

#if 0
  inline void setRGB(const coord_t x, const coord_t y, const half3 rgb)
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    setR(x, y, rgb[0]);
    setG(x, y, rgb[1]);
    setB(x, y, rgb[2]);
  }
#endif

  inline count_t getN(const coord_t x, const coord_t y) const
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    uint32_t n1 = N1 ? N1[y * width + x] : 0;
    uint32_t n0 = N0 ? N0[y * width + x] : 0;
    if (N1 && N0 && n1 == ~(uint32_t(0)) && n0 == ~(uint32_t(0)))
    {
      return maxiters;
    }
    else if (N0 && n0 == ~(uint32_t(0)))
    {
      return maxiters;
    }
    return count_t((uint64_t(n1) << 32) | n0) - Nbias;
  }

  inline float getNF(const coord_t x, const coord_t y) const
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    return NF ? NF[y * width + x] : 0;
  }

  inline float getT(const coord_t x, const coord_t y) const
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    return T ? T[y * width + x] : 0;
  }

  inline float getDEX(const coord_t x, const coord_t y) const
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    return DEX ? DEX[y * width + x] : 0;
  }

  inline float getDEY(const coord_t x, const coord_t y) const
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    return DEY ? DEY[y * width + x] : 0;
  }

  inline complex<float> getDE(const coord_t x, const coord_t y) const
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    return complex<float>(getDEX(x, y), getDEY(x, y));
  }

  inline half getR(const coord_t x, const coord_t y) const
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    return RGB ? RGB[3 * (y * width + x) + 0] : half(0);
  }

  inline half getG(const coord_t x, const coord_t y) const
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    return RGB ? RGB[3 * (y * width + x) + 1] : half(0);
  }

  inline half getB(const coord_t x, const coord_t y) const
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    return RGB ? RGB[3 * (y * width + x) + 2] : half(0);
  }

#if 0
  inline half3 getRGB(const coord_t x, const coord_t y) const
  {
    assert(0 <= x);
    assert(x < width);
    assert(0 <= y);
    assert(y < height);
    return half3(getR(x, y), getG(x, y), getB(x, y));
  }
#endif

  void saveEXR(const std::string &filename, const channel_mask_t channels, const int threads = 1, const std::string &metadata = "", const std::string &kf2plus_metadata = "") const;
};
