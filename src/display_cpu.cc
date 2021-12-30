// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "colour.h"
#include "display_cpu.h"
#include "map.h"
#include "parallel.h"

display_cpu::display_cpu(const colour *clr)
: width(0)
, height(0)
, do_clear(false)
, subframes(0)
, RGB(0)
, clr(clr)
{
}

display_cpu::~display_cpu()
{
}

void display_cpu::resize(coord_t new_width, coord_t new_height)
{
  if (new_width == width && new_height == height)
  {
    return;
  }
  width = new_width;
  height = new_height;
  RGB.resize(width * height);
  do_clear = true;
  subframes = 0;
}

void display_cpu::set_colour(const colour *new_clr)
{
  clr = new_clr;
}

void display_cpu::clear()
{
  do_clear = true;
}

void display_cpu::accumulate(const map &out)
{
  assert(out.width == width);
  assert(out.height == height);
  if (do_clear)
  {
    subframes = 0;
  }
  volatile bool running = true;
  parallel2d(std::thread::hardware_concurrency(), 0, width, 32, 0, height, 32, &running, [&](coord_t i, coord_t j) -> void
  {
    count_t n = out.getN(i, j);
    float t = out.getT(i, j);
    float f = out.getNF(i, j);
    complex<float> de = out.getDE(i, j);
    vec3 rgb = clr->rgb(n, vec2(t, f), vec2(de.x, de.y));
    if (do_clear)
    {
      RGB[j * width + i] = rgb;
    }
    else
    {
      RGB[j * width + i] += rgb;
    }
  });
  do_clear = false;
  subframes++;
}

void display_cpu::get_rgb(map &out) const
{
  assert(out.RGB);
  assert(out.width == width);
  assert(out.height == height);
  volatile bool running = true;
  parallel2d(std::thread::hardware_concurrency(), 0, width, 32, 0, height, 32, &running, [&](coord_t i, coord_t j) -> void
  {
    vec3 rgb = RGB[j * width + i] / float(subframes);
    out.setR(i, j, rgb[0]);
    out.setG(i, j, rgb[1]);
    out.setB(i, j, rgb[2]);
  });
}
