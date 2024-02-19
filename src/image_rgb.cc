// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2024 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <cstdio>

#include "exr.h"

#ifdef HAVE_EXR

static const char fraktaler3[] = "Fraktaler3";

#endif

#include "image_rgb.h"
#include "png.h"
#include "render.h"

image_rgb::image_rgb(coord_t width, coord_t height)
: image(width, height)
{
  RGBA = new float[4 * width * height]; // FIXME tile based alpha
}

image_rgb::image_rgb(image_rgb &source, bool vflip)
: image_rgb(source.width, source.height)
{
  std::lock_guard<std::mutex> lock(source.mutex);
  for (coord_t y = 0; y < height; ++y)
  {
    for (coord_t x = 0; x < width; ++x)
    {
      coord_t z = (y * width + x) * 4;
      coord_t w = vflip ? ((height - 1 - y) * width + x) * 4 : z;
      float A = source.RGBA[w + 3];
      if (A == 0)
      {
        A = 1;
      }
      RGBA[z] = source.RGBA[w] / A; z++; w++;
      RGBA[z] = source.RGBA[w] / A; z++; w++;
      RGBA[z] = source.RGBA[w] / A; z++; w++;
      RGBA[z] = 1;
    }
  }
}

image_rgb::~image_rgb()
{
  delete[] RGBA;
}

void image_rgb::clear()
{
  std::memset(RGBA, 0, 4 * width * height * sizeof(*RGBA));
}

void image_rgb::blit(coord_t tx, coord_t ty, const struct tile *t)
{
  if (t && t->RGB)
  {
    std::lock_guard<std::mutex> lock(mutex); // FIXME tile based locking
    for (coord_t j = 0; j < t->height; ++j)
    {
      coord_t y = ty * t->height + j;
      if (0 <= y && y < height)
      {
        for (coord_t i = 0; i < t->width; ++i)
        {
          coord_t x = tx * t->width + i;
          if (0 <= x && x < width)
          {
            coord_t k = (j * t->width + i) * 3;
            coord_t z = (y * width + x) * 4;
            RGBA[z++] += t->RGB[k++];
            RGBA[z++] += t->RGB[k++];
            RGBA[z++] += t->RGB[k++];
            RGBA[z++] += 1;
          }
        }
      }
    }
  }
}

bool image_rgb::save_exr(const std::string &filename, int threads, const std::string &metadata)
{
  std::lock_guard<std::mutex> lock(mutex); // FIXME tile based locking
#ifndef HAVE_EXR
  (void) filename;
  (void) threads;
  (void) metadata;
  (void) kf2plus_metadata;
#else
  try
  {
    setGlobalThreadCount(threads);
    Header header(width, height);
    if (metadata != "") header.insert(fraktaler3, StringAttribute(metadata));
    header.channels().insert("R", Channel(IMF::FLOAT));
    header.channels().insert("G", Channel(IMF::FLOAT));
    header.channels().insert("B", Channel(IMF::FLOAT));
    OutputFile of(filename.c_str(), header);
    FrameBuffer fb;
    fb.insert("R", Slice(IMF::FLOAT, (char *)(RGBA + 0), sizeof(*RGBA) * 4, sizeof(*RGBA) * 4 * width));
    fb.insert("G", Slice(IMF::FLOAT, (char *)(RGBA + 1), sizeof(*RGBA) * 4, sizeof(*RGBA) * 4 * width));
    fb.insert("B", Slice(IMF::FLOAT, (char *)(RGBA + 2), sizeof(*RGBA) * 4, sizeof(*RGBA) * 4 * width));
    of.setFrameBuffer(fb);
    of.writePixels(height);
    return true;
  }
  catch (...)
#endif
  {
    return false;
  }
}

image_rgb8::image_rgb8(coord_t width, coord_t height)
: width(width)
, height(height)
{
  RGB = new unsigned char[3 * width * height];
}

image_rgb8::~image_rgb8()
{
  delete[] RGB;
}

inline float linear_to_srgb(float c) noexcept
{
  c = std::fmin(std::fmax(c, 0.0f), 1.0f);
  if (c <= 0.0031308f)
  {
    return 12.92f * c;
  }
  else
  {
    return 1.055f * std::pow(c, 1.0f / 2.4f) - 0.055f;
  }
}

image_rgb8::image_rgb8(image_rgb &source, bool vflip)
: image_rgb8(source.width, source.height)
{
  std::lock_guard<std::mutex> lock(source.mutex);
  for (coord_t y = 0; y < height; ++y)
  {
    for (coord_t x = 0; x < width; ++x)
    {
      coord_t z = (y * width + x) * 3;
      coord_t w = vflip ? ((height - 1 - y) * width + x) * 4 : (y * width + x) * 4;;
      float A = source.RGBA[w + 3];
      if (A == 0)
      {
        A = 1;
      }
      RGB[z] = std::fmin(std::fmax(0.0f, 255.0f * linear_to_srgb(source.RGBA[w] / A)), 255.0f); z++; w++;
      RGB[z] = std::fmin(std::fmax(0.0f, 255.0f * linear_to_srgb(source.RGBA[w] / A)), 255.0f); z++; w++;
      RGB[z] = std::fmin(std::fmax(0.0f, 255.0f * linear_to_srgb(source.RGBA[w] / A)), 255.0f); z++; w++;
    }
  }
}

bool image_rgb8::save_png(const std::string &filename, const std::string &metadata)
{
  return save_png_rgb8(filename, RGB, width, height, metadata);
}
