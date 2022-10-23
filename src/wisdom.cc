// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <algorithm>
#include <atomic>
#include <iostream>
#include <map>

#include <toml.hpp>

#ifdef HAVE_CL
#define CL_TARGET_OPENCL_VERSION 200
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#ifdef HAVE_CLEW
#include "clew.h"
#else
#include <CL/cl.h>
#endif
#endif

#include "engine.h"
#include "param.h"
#include "render.h"
#include "wisdom.h"

TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(whardware, name, platform, device)
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(wdevice, platform, device, speed)
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(wtype, mantissa, exponent, device)
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(wisdom, hardware, type)

std::istream &operator>>(std::istream &ifs, wisdom &w)
{
  auto t = toml::parse(ifs);
  w.hardware = toml::find<std::map<std::string, std::vector<whardware>>>(t, "hardware");
  w.type = toml::find<std::map<std::string, wtype>>(t, "type");
  return ifs;
}

std::ostream &operator<<(std::ostream &ofs, const wisdom &w)
{
  const toml::value t{{"hardware", w.hardware}, {"type", w.type}};
  return ofs << t;
}

wisdom wisdom_load(const std::string &filename, bool &success)
{
  wisdom w;
  std::ifstream ifs;
  ifs.exceptions(std::ifstream::badbit);
  ifs.open(filename, std::ios_base::binary);
  if (! ifs.is_open())
  {
    success = false;
    return w;
  }
  ifs >> w;
  ifs.close();
  success = true;
  return w;
}

void wisdom_save(const wisdom &w, const std::string &filename)
{
  std::ofstream ofs;
  ofs.exceptions(std::ofstream::badbit | std::ofstream::failbit);
  ofs.open(filename, std::ios_base::binary);
  ofs << w;
}

wisdom wisdom_enumerate(bool use_opencl)
{
  wisdom w;
  int cpu = 0;
  (void) use_opencl;
#ifdef HAVE_CL
  if (use_opencl)
  {
    int gpu = 0;
    int other = 0;
    cl_platform_id platform_id[64];
    cl_uint platform_ids;
    if (CL_SUCCESS == clGetPlatformIDs(64, &platform_id[0], &platform_ids))
    {
      for (int platform = 0; platform < (int) platform_ids; ++platform)
      {
        char buf[1024];
        buf[0] = 0;
        if (CL_SUCCESS == clGetPlatformInfo(platform_id[platform], CL_PLATFORM_NAME, 1024, &buf[0], 0))
        {
          std::string platform_name = buf;
          cl_device_id device_id[64];
          cl_uint device_ids;
          if (CL_SUCCESS == clGetDeviceIDs(platform_id[platform], CL_DEVICE_TYPE_ALL, 64, &device_id[0], &device_ids))
          {
            for (int device = 0; device < (int) device_ids; ++device)
            {
              buf[0] = 0;
              if (CL_SUCCESS == clGetDeviceInfo(device_id[device], CL_DEVICE_NAME, 1024, &buf[0], 0))
              {
                std::string device_name = buf;
                cl_device_type device_type;
                if (CL_SUCCESS == clGetDeviceInfo(device_id[device], CL_DEVICE_TYPE, sizeof(device_type), &device_type, 0))
                {
                  std::ostringstream s;
                  switch (device_type)
                  {
                    case CL_DEVICE_TYPE_CPU:
                      s << "cpu" << cpu++;
                      break;
                    case CL_DEVICE_TYPE_GPU:
                      s << "gpu" << gpu++;
                      break;
                    default:
                      s << "other" << other++;
                      break;
                  }
                  w.hardware[s.str()] = std::vector<whardware>{ whardware{ device_name + " (" + platform_name + ")", platform, device } };
                }
              }
            }
          }
        }
      }
    }
  }
#endif

  {
    std::ostringstream s;
    s << "cpu" << cpu++;
    w.hardware[s.str()] = std::vector<whardware>{ whardware{ s.str(), -1, 0 } };
  }
  int nt_mantissa[] = { 0, 24, 53, 0, 24, 32, 113 };
  int nt_exponent[] = { 0,  8, 11, 0, 24, 31,  15 };
  double nt_speed[] = { 0, 110, 100, 10, 30, 20, 5 };
  {
    using std::isinf;
    using std::ldexp;
    int mantissa_bits = 0;
    int exponent_bits = 0;
    const long double one = 1;
    for (int bits = 0; ; bits++)
    {
      if (one + ldexp(one, -bits) == one)
      {
        mantissa_bits = bits;
        break;
      }
    }
    for (int bits = 0; ; bits++)
    {
      if (isinf(ldexp(one, 1 << bits)))
      {
        exponent_bits = bits + 1; // signed
        break;
      }
    }
    nt_mantissa[nt_longdouble] = mantissa_bits;
    nt_exponent[nt_longdouble] = exponent_bits;
  }

  for (const auto nt : { nt_float, nt_double, nt_longdouble, nt_floatexp, nt_softfloat
#ifdef HAVE_FLOAT128
  , nt_float128
#endif
  })
  {
    w.type[nt_string[nt]] = { nt_mantissa[nt], nt_exponent[nt], { } };
    for (const auto & taghws : w.hardware)
    {
      const auto & hws = taghws.second;
      for (const auto & hw : hws)
      {
        if (nt == nt_longdouble
#ifdef HAVE_FLOAT128
         || nt == nt_float128
#endif
        )
        {
          if (hw.platform == -1)
          {
            w.type[nt_string[nt]].device.push_back(wdevice{ hw.platform, hw.device, nt_speed[nt] });
          }
        }
        else
        {
          w.type[nt_string[nt]].device.push_back(wdevice{ hw.platform, hw.device, nt_speed[nt] });
        }
      }
    }
  }
  return w;
}

void wisdom_default(wisdom &w)
{
  w = wisdom_enumerate(false);
}

bool comparing_speed(const wlookup &a, const wlookup &b) 
{
  return a.speed < b.speed;
}

number_type nt_from_string(const std::string &s)
{
  for (int i = 0; i <
#ifdef HAVE_FLOAT128
  7
#else
  6
#endif
  ; ++i)
  {
    if (s == std::string(nt_string[i]))
    {
      return number_type(i);
    }
  }
  return nt_none;
}

wlookup wisdom_lookup(const wisdom &w, const std::set<number_type> &available, count_t pixel_spacing_exponent, count_t pixel_spacing_precision)
{
  // find suitable number types
  std::vector<wlookup> candidates;
  for (const auto & ntstype : w.type)
  {
    const auto & nts = ntstype.first;
    const auto & type = ntstype.second;
    number_type nt = nt_from_string(nts);
    if (pixel_spacing_exponent + 16 < (count_t(1) << type.exponent) >> 1 &&
        pixel_spacing_precision < type.mantissa &&
        available.find(nt) != available.end())
    {
      // pick fastest device in each hardware group
      wlookup candidate = { nt_from_string(nts), type.mantissa, type.exponent, 0.0, { } };
      for (const auto & namehardwares : w.hardware)
      {
        const auto & hardwares = namehardwares.second;
        int index = -1;
        double speed = 0.0;
        for (const auto & hardware : hardwares)
        {
          for (int ix = 0; ix < (int) type.device.size(); ++ix)
          {
            const auto & device = type.device[ix];
            if (hardware.platform == device.platform &&
                hardware.device == device.device)
            {
              if (speed <= device.speed)
              {
                speed = device.speed;
                index = ix;
              }
            }
          }
          if (index >= 0 && speed > 0.0)
          {
            candidate.speed += speed;
            candidate.device.push_back(type.device[index]);
          }
        }
      }
      if (! candidate.device.empty() && candidate.speed > 0.0)
      {
        candidates.push_back(candidate);
      }
    }
  }
  // pick fastest number type
  if (candidates.empty())
  {
    return wlookup{ nt_none, 0, 0, 0.0, { } };
  }
  return *std::max_element(candidates.begin(), candidates.end(), comparing_speed);
}

#if 0  
  /* check if device supports double precision */
  cl_uint dvecsize = 0;
  cl_int status = clGetDeviceInfo(device_id[device_index], CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(dvecsize), &dvecsize, 0);
  if (status != CL_SUCCESS)
  {
    dvecsize = 0;
  }
#endif

struct wisdom_hooks : public hooks
{
  int platform;
  int device;
  coord_t image_width;
  coord_t image_height;
  std::atomic<count_t> pixels;
  double seconds;
  float min, max;
  std::chrono::high_resolution_clock::time_point start_time;
  wisdom_hooks(int platform, int device, coord_t image_width, coord_t image_height)
  : platform(platform)
  , device(device)
  , image_width(image_width)
  , image_height(image_height)
  , pixels(0)
  , seconds(0.0)
  , min(1.0f/0.0f)
  , max(-1.0f/0.0f)
  {
  }
  virtual ~wisdom_hooks() { }
  virtual void pre_tile(int platformx, int devicex, int x, int y, int subframe) override
  {
    (void) x;
    (void) y;
    (void) subframe;
    if (platformx == platform && devicex == device)
    {
      start_time = std::chrono::high_resolution_clock::now();
    }
  }
  virtual void post_download(int platformx, int devicex, int x, int y, int subframe) override
  {
    (void) x;
    (void) y;
    (void) subframe;
    if (platformx == platform && devicex == device)
    {
      const auto end_time = std::chrono::high_resolution_clock::now();
      seconds += std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    }
  }
  // read tile RGB data to force transfer from mapped device
  virtual void tile(int platformx, int devicex, int x, int y, int subframe, const struct tile *data) override
  {
    (void) subframe;
    if (platformx == platform && devicex == device)
    {
      coord_t pixelsx = 0;
      for (coord_t j = 0; j < data->height; ++j)
      {
        if (y * data->height + j < image_height)
        {
          for (coord_t i = 0; i < data->width; ++i)
          {
            if (x * data->width + i < image_width)
            {
              for (coord_t c = 0; c < 3; ++c)
              {
                coord_t k = (j * data->width + i) * 3 + c;
                float v = data->RGB[k];
                min = std::min(min, v);
                max = std::max(max, v);
              }
              pixelsx++;
            }
          }
        }
      }
      pixels += pixelsx;
    }
  }
};

double wisdom_benchmark_device(const wlookup &l, const param &par0, volatile bool *running)
{
  std::fprintf(stderr, "%d.%d %s %d:%d ", l.device[0].platform, l.device[0].device, nt_string[l.nt], l.mantissa, l.exponent);
  count_t width = 32;
  count_t height = 18;
  count_t tile_width = width;
  count_t tile_height = height;
  count_t subframes = 16;
  double seconds = 0.0;
  double target_seconds = 10.0;
  double speed = 0.0;
  int bad = 0;
  std::vector<progress_t> progress(par0.p.formula.per.size() * 2 + 2 * l.device.size(), 0);
  do
  {
    std::fprintf(stderr, ".");
    width *= 2;
    height *= 2;
    if (tile_width < 2048)
    {
      tile_width *= 2;
      tile_height *= 2;
    }
    param par = par0;
    par.p.image.width = width;
    par.p.image.height = height;
    par.p.image.subframes = subframes;
    par.p.opencl.tile_width = tile_width;
    par.p.opencl.tile_height = tile_height;
    try
    {
      wisdom_hooks h(l.device[0].platform, l.device[0].device, width, height);
      render(l, par, &h, true, &progress[0], running);
      if (! (h.min < h.max))
      {
        bad++;
        if (bad >= 3)
        {
          speed = -1;
          std::fprintf(stderr, "(%g!<!%g)", h.min, h.max);
          break;
        }
      }
      if (running && h.min < h.max)
      {
        seconds = h.seconds;
        speed = h.pixels / seconds;
      }
    }
    catch (...)
    {
      break;
    }
  }
  while (seconds < target_seconds && *running);
  if (*running)
  {
    std::fprintf(stderr, " %.2f\n", speed);
    return speed;
  }
  else
  {
    return 0;
  }
}

wisdom wisdom_benchmark(const wisdom &wi, volatile bool *running)
{
  param par;
  par.from_string(
    "program = \"fraktaler-3\"\n"
    "location.real = \"0.352465816656845732\"\n"
    "location.imag = \"0.392188990843255425\"\n"
    "location.zoom = \"4.1943021e6\"\n"
    "algorithm.reuse_reference = true\n"
    "\n"
    "[[formula]]\n"
    "power = 2\n"
    "neg_y = false\n"
    "neg_x = false\n"
    "abs_y = false\n"
    "abs_x = false\n"
    "[[formula]]\n"
    "power = 2\n"
    "neg_y = false\n"
    "neg_x = false\n"
    "abs_y = true\n"
    "abs_x = false\n"
    "[[formula]]\n"
    "power = 2\n"
    "neg_y = false\n"
    "neg_x = false\n"
    "abs_y = false\n"
    "abs_x = false\n"
    "[[formula]]\n"
    "power = 2\n"
    "neg_y = false\n"
    "neg_x = false\n"
    "abs_y = true\n"
    "abs_x = true\n"
    "\n"
    "\n"
  );
  set_reference_to_image_center(par);
  wisdom wo;
  wo.hardware = wi.hardware;
  for (const auto & ntstype : wi.type)
  {
    const auto & nts = ntstype.first;
    const auto & type = ntstype.second;
    if (! running)
    {
      break;
    }
    wo.type[nts] = { type.mantissa, type.exponent, { } };
    number_type nt = nt_from_string(nts);
    if (nt != nt_none)
    {
      for (const auto device : type.device)
      {
        if (! running)
        {
          break;
        }
        wlookup l = { nt, type.mantissa, type.exponent, 0.0, { device } };
        double speed = wisdom_benchmark_device(l, par, running);
        if (running)
        {
          wo.type[nts].device.push_back(wdevice{ device.platform, device.device, speed });
        }
      }
    }
  }
  return wo;
}
