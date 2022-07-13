// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <chrono>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <thread>

#define CL_TARGET_OPENCL_VERSION 200
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#ifdef HAVE_CLEW
#include "clew.h"
#else
#include <CL/cl.h>
#endif

#include "bla.h"
#include "colour.h"
#include "display_cpu.h"
#include "engine.h"
#include "floatexp.h"
#include "map.h"
#include "parallel.h"
#include "param.h"
#include "softfloat.h"
#include "stats.h"
#include "types.h"
#include "version.h"

#include "cl-pre.h"
#include "cl-post.h"

extern std::vector<std::vector<complex<softfloat>>> Zsf;
extern std::vector<std::vector<complex<floatexp>>> Zfe;
extern std::vector<std::vector<complex<double>>> Zd;
extern std::vector<std::vector<complex<float>>> Zf;

extern std::vector<blasR2<softfloat>> Bsf;
extern std::vector<blasR2<floatexp>> Bfe;
extern std::vector<blasR2<double>> Bd;
extern std::vector<blasR2<float>> Bf;

std::string hybrid_perturb(const std::vector<phybrid1> &per)
{
  std::ostringstream s;
  s << "{\n";
  s << "  struct complex Z = { ref[config->ref_start[phase] + 2 * m], ref[config->ref_start[phase] + 2 * m + 1] };\n";
  s << "  real X = Z.x;\n";
  s << "  real Y = Z.y;\n";
  s << "  struct dual x = z.x;\n";
  s << "  struct dual y = z.y;\n";
  s << "  struct complexdual W = complexdual_add_complex_complexdual(Z, z);\n";
  s << "  struct complex B = Z;\n";
  s << "  switch (n % " << per.size() << ")\n";
  s << "  {\n";
  for (size_t k = 0; k < per.size(); ++k)
  {
    s << "  case " << k << ":\n";
    s << "    {\n";
    if (per[k].abs_x)
    {
      s << "      x = dual_diffabs_real_dual(X, x);\n";
      s << "      W.x = dual_abs_dual(W.x);\n";
      s << "      B.x = real_abs_real(B.x);\n";
    }
    if (per[k].abs_y)
    {
      s << "      y = dual_diffabs_real_dual(Y, y);\n";
      s << "      W.y = dual_abs_dual(W.y);\n";
      s << "      B.y = real_abs_real(B.y);\n";
    }
    if (per[k].neg_x)
    {
      s << "      x = dual_neg_dual(x);\n";
      s << "      W.x = dual_neg_dual(W.x);\n";
      s << "      B.x = real_neg_real(B.x);\n";
    }
    if (per[k].neg_y)
    {
      s << "      y = dual_neg_dual(y);\n";
      s << "      W.y = dual_neg_dual(W.y);\n";
      s << "      B.y = real_neg_real(B.y);\n";
    }
    s << "      struct complexdual P = { x, y };\n";
    s << "      struct complexdual S = { { real_from_int(0), { real_from_int(0), real_from_int(0) } }, { real_from_int(0), { real_from_int(0), real_from_int(0) } } };\n";
    s << "      struct complex Bn[" << per[k].power << "];\n";
    s << "      Bn[0].x = real_from_int(1); Bn[0].y = real_from_int(0);\n";
    for (int i = 1; i < per[k].power; ++i)
    {
      s << "      Bn["  << i << "] = complex_mul_complex_complex(Bn[" << (i - 1) << "], B);\n";
    }
    s << "      struct complexdual Wi = S; Wi.x.x = real_from_int(1);\n";
    for (int i = 0; i < per[k].power; ++i)
    {
      s << "      S = complexdual_add_complexdual_complexdual(S, complexdual_mul_complexdual_complex(Wi, Bn[" << (per[k].power - 1 - i) << "]));\n";
      if (i != per[k].power - 1)
      {
        s << "      Wi = complexdual_mul_complexdual_complexdual(Wi, W);\n";
      }
    }
    s << "      z = complexdual_add_complexdual_complexdual(complexdual_mul_complexdual_complexdual(P, S), c);\n";
    s << "    }\n";
    s << "    break;\n";
  }
  s << "  }\n";
  s << "}\n";
  return s.str();
}

struct opencl_error
{
  const char *file;
  int err, loc;
};

void e(const char *file, int loc, int err)
{
  if (err == CL_SUCCESS) { return; }
  throw opencl_error{ file, loc, err };
}

#define E(err) e("cl.cc", __LINE__, err)

#define MAX_PHASES 64
#define MAX_LEVELS 64

// used to find maximum possible size of config structure
union ufloat
{
  float f;
  double d;
  floatexp fe;
  softfloat sf;
};

template <typename real>
struct mat2_cl
{
  real a, b, c, d;
};

template <typename real>
struct blaR2_cl
{
  struct mat2_cl<real> A, B;
  real r2;
  cl_long l;
};

template <typename real>
struct config_cl
{
  cl_long size;
  cl_long number_type;
  /* shape */
  cl_long height;
  cl_long width;
  cl_long subframes;
  cl_long frame;
  /* bailout */
  cl_long Iterations;
  real ER2;
  cl_long PerturbIterations;
  /* transform */
  cl_long transform_exponential_map;
  struct mat2_cl<real> transform_K;
  real pixel_spacing;
  real offset_x;
  real offset_y;
  /* ref layout */
  cl_long number_of_phases;
  cl_long ref_size[MAX_PHASES];
  cl_long ref_start[MAX_PHASES];
  /* bla layout */
  cl_long bla_size[MAX_PHASES];
  cl_long bla_levels[MAX_PHASES];
  cl_long bla_start[MAX_PHASES][MAX_LEVELS];
};

void opencl_thread(map &out, param &par, progress_t *progress, bool *running, bool *ended)
{

  /* get platform */
  int platform_index = par.p.opencl.platform;
  int device_index = par.p.opencl.device;
  cl_platform_id platform_id[64];
  cl_uint platform_ids;
  E(clGetPlatformIDs(64, &platform_id[0], &platform_ids));
  if (! (0 <= platform_index && platform_index < (int) platform_ids))
  {
    fprintf(stderr, "error: platform out of range\n");
    abort();
  }
  char buf[1024];
  buf[0] = 0;
  E(clGetPlatformInfo(platform_id[platform_index], CL_PLATFORM_VENDOR, 1024, &buf[0], 0));
  fprintf(stderr, "platform vendor: %s\n", buf);
  buf[0] = 0;
  E(clGetPlatformInfo(platform_id[platform_index], CL_PLATFORM_VERSION, 1024, &buf[0], 0));
  fprintf(stderr, "platform version: %s\n", buf);
  /* get device */
  cl_device_id device_id[64];
  cl_uint device_ids;
  E(clGetDeviceIDs(platform_id[platform_index], CL_DEVICE_TYPE_ALL, 64, &device_id[0], &device_ids));
  if (! (0 <= device_index && device_index < (int) device_ids))
  {
    fprintf(stderr, "error: device out of range\n");
    abort();
  }
  buf[0] = 0;
  E(clGetDeviceInfo(device_id[device_index], CL_DEVICE_NAME, 1024, &buf[0], 0));
  fprintf(stderr, "device name: %s\n", buf);
  /* check if device supports double precision */
  cl_uint dvecsize = 0;
  cl_int status = clGetDeviceInfo(device_id[device_index], CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(dvecsize), &dvecsize, 0);
  if (status != CL_SUCCESS)
  {
    dvecsize = 0;
  }
  /* create context */
  cl_context_properties properties[] =
    {
      CL_CONTEXT_PLATFORM, (cl_context_properties) platform_id[platform_index]
    , 0
    };
  cl_int err;
  cl_context context = clCreateContext(properties, 1, &device_id[device_index], NULL, NULL, &err);
  if (! context) { E(err); }
  cl_command_queue commands = clCreateCommandQueue(context, device_id[device_index], 0, &err);
  if (! commands) { E(err); }

  const std::string body = hybrid_perturb(par.p.formula.per);
  unsigned int src_cl_body_cl_len = strlen(body.c_str());
  unsigned int source_len = src_cl_pre_cl_len + src_cl_body_cl_len + src_cl_post_cl_len + 1;
  char *source = (char *) malloc(source_len);
  if (! source)
  {
    fprintf(stderr, "error: out of memory\n");
    abort();
  }
  memcpy(source, src_cl_pre_cl, src_cl_pre_cl_len);
  memcpy(source + src_cl_pre_cl_len, body.c_str(), src_cl_body_cl_len);
  memcpy(source + src_cl_pre_cl_len + src_cl_body_cl_len, src_cl_post_cl, src_cl_post_cl_len);
  source[source_len - 1] = 0;
  cl_program program = clCreateProgramWithSource(context, 1, const_cast<const char **>(&source), 0, &err);
  if (! program) { E(err); }

  coord_t width = par.p.image.width / par.p.image.subsampling;
  coord_t height = par.p.image.height / par.p.image.subsampling;
  size_t rgb_bytes = sizeof(float) * width * height * 3;
  float *rgb_host = (float *) malloc(rgb_bytes);
  if (! rgb_host)
  {
    fprintf(stderr, "error: out of memory\n");
    abort();
  }

  cl_mem config_device = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(struct config_cl<ufloat>), 0, &err);
  if (! config_device) { E(err); }

  cl_mem rgb_device = clCreateBuffer(context, CL_MEM_READ_WRITE, rgb_bytes, 0, &err);
  if (! rgb_device) { E(err); }

  size_t raw_bytes = sizeof(float) * width * height;
  cl_mem N0_device = 0, N1_device = 0, NF_device = 0, T_device = 0, DEX_device = 0, DEY_device = 0;
  if (par.p.image.subframes == 1)
  {
    N0_device = clCreateBuffer(context, CL_MEM_WRITE_ONLY, raw_bytes, 0, &err);
    if (! N0_device) { E(err); }
    N1_device = clCreateBuffer(context, CL_MEM_WRITE_ONLY, raw_bytes, 0, &err);
    if (! N1_device) { E(err); }
    NF_device = clCreateBuffer(context, CL_MEM_WRITE_ONLY, raw_bytes, 0, &err);
    if (! NF_device) { E(err); }
    T_device = clCreateBuffer(context, CL_MEM_WRITE_ONLY, raw_bytes, 0, &err);
    if (! T_device) { E(err); }
    DEX_device = clCreateBuffer(context, CL_MEM_WRITE_ONLY, raw_bytes, 0, &err);
    if (! DEX_device) { E(err); }
    DEY_device = clCreateBuffer(context, CL_MEM_WRITE_ONLY, raw_bytes, 0, &err);
    if (! DEY_device) { E(err); }
  }

  using std::ceil;
  int threads = std::thread::hardware_concurrency();
  floatexp Zoom = par.zoom;
  floatexp ZoomedOut = 1 / 65536.0;
  count_t nframes = ceil(double(log(Zoom / ZoomedOut) / log(par.p.render.zoom_out_factor)));
  count_t start_frame = par.p.render.start_frame;
  count_t end_frame = par.p.render.frame_count == 0 ? nframes : par.p.render.start_frame + par.p.render.frame_count;
  nframes = end_frame - start_frame;
  if (nframes == 0)
  {
    nframes = 1;
  }
  const count_t count = par.p.formula.per.size();

  number_type nt_kernel = nt_none;
  cl_kernel kernel = 0;

  if (! par.p.render.zoom_out_sequence)
  {
    start_frame = 0;
    end_frame = 1;
  }
  for (count_t frame = start_frame; frame < end_frame; ++frame)
  {
    cl_mem ref_device = 0; // FIXME reuse
    cl_mem bla_device = 0; // FIXME reuse

    par.zoom = Zoom / pow(floatexp(par.p.render.zoom_out_factor), frame);
    progress[0] = (frame - start_frame) / progress_t(nframes);

    bool ref_ended = false;
    {
      stats sta;
      reference_thread(sta, par, false, &progress[1], running, &ref_ended);
    }

    if (nt_kernel != nt_current)
    {        
      std::ostringstream optionss;
      optionss << "-DNUMBER_TYPE=" << int(nt_current);
      err = clBuildProgram(program, 1, &device_id[device_index], optionss.str().c_str(), 0, 0);
      if (err != CL_SUCCESS)
      {
        char *error_log = (char *) malloc(1000000);
        error_log[0] = 0;
        E(clGetProgramBuildInfo(program, device_id[device_index], CL_PROGRAM_BUILD_LOG, 1000000, &error_log[0], 0));
        fprintf(stderr, "error: compile failed:\n%s\n", error_log);
        fprintf(stderr, "error: corresponding source:\n%s\n", source);
        E(err);
      }
      kernel = clCreateKernel(program, "fraktaler3", &err);
      if (! kernel) { E(err); }
      nt_kernel = nt_current;
    }

    complex<mpreal> offset;
    offset.x.set_prec(par.center.x.get_prec());
    offset.y.set_prec(par.center.y.get_prec());
    offset = par.center - par.reference;

    cl_event ready;
    switch (nt_current)
    {
      case nt_float:
        {
          struct config_cl<float> config_host =
            { sizeof(config_host)
            , nt_float
            , height
            , width
            , par.p.image.subframes
            , frame
            , par.p.bailout.iterations
            , float(par.p.bailout.escape_radius * par.p.bailout.escape_radius)
            , par.p.bailout.maximum_perturb_iterations
            , par.p.transform.exponential_map
            , { float(par.transform.x[0][0]), float(par.transform.x[0][1]), float(par.transform.x[1][0]), float(par.transform.x[1][1]) }
            , float(4 / par.zoom / height)
            , float(offset.x)
            , float(offset.y)
            , cl_long(par.p.formula.per.size())
            // ...
            };
    
          /* reference layout */
          cl_long ref_size = 0;
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            ref_size += config_host.ref_size[phase] = Zf[phase].size();
          }
          cl_long ref_bytes = 2 * ref_size * sizeof(float);
          config_host.ref_start[0] = 0;
          for (int phase = 1; phase < config_host.number_of_phases; ++phase)
          {
            config_host.ref_start[phase] = config_host.ref_start[phase - 1] + 2 * config_host.ref_size[phase - 1];
          }
    
          /* bla layout */
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            config_host.bla_size[phase] = Bf[phase].M;
            config_host.bla_levels[phase] = Bf[phase].L;
          }
          cl_long bla_start = 0;
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            for (int level = 0; level < config_host.bla_levels[phase]; ++level)
            {
              config_host.bla_start[phase][level] = bla_start;
              bla_start += Bf[phase].b[level].size();
            }
          }
          cl_long bla_bytes = bla_start * sizeof(struct blaR2_cl<float>);
        
          /* upload config */
          E(clEnqueueWriteBuffer(commands, config_device, CL_FALSE, 0, sizeof(config_host), &config_host, 0, 0, &ready));
    
          /* upload reference */
          ref_device = clCreateBuffer(context, CL_MEM_READ_ONLY, ref_bytes, 0, &err);
          if (! ref_device) { E(err); }
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            const cl_long start_bytes = config_host.ref_start[phase] * sizeof(float);
            const cl_long size_bytes = config_host.ref_size[phase] * 2 * sizeof(float);
            const void *ptr = &Zf[phase][0];
            cl_event done;
            E(clEnqueueWriteBuffer(commands, ref_device, CL_FALSE, start_bytes, size_bytes, ptr, 1, &ready, &done));
            E(clReleaseEvent(ready));
            ready = done;
          }
    
          /* upload bla */
          bla_device = clCreateBuffer(context, CL_MEM_READ_ONLY, bla_bytes, 0, &err);
          if (! bla_device) { E(err); }
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            for (int level = 0; level < config_host.bla_levels[phase]; ++level)
            {
              const cl_long start_bytes = config_host.bla_start[phase][level] * sizeof(blaR2_cl<float>);
              const cl_long size_bytes = Bf[phase].b[level].size() * sizeof(blaR2_cl<float>);
              const void *ptr = &Bf[phase].b[level][0];
              cl_event done;
              E(clEnqueueWriteBuffer(commands, bla_device, CL_FALSE, start_bytes, size_bytes, ptr, 1, &ready, &done));
              E(clReleaseEvent(ready));
              ready = done;
            }
          }
        }
        break;

      case nt_double:
        {
          if (dvecsize == 0)
          {
            fprintf(stderr, "error: device does not support double precision\n");
            abort();
          }
          struct config_cl<double> config_host =
            { sizeof(config_host)
            , nt_double
            , height
            , width
            , par.p.image.subframes
            , frame
            , par.p.bailout.iterations
            , double(par.p.bailout.escape_radius * par.p.bailout.escape_radius)
            , par.p.bailout.maximum_perturb_iterations
            , par.p.transform.exponential_map
            , { double(par.transform.x[0][0]), double(par.transform.x[0][1]), double(par.transform.x[1][0]), double(par.transform.x[1][1]) }
            , double(4 / par.zoom / height)
            , double(offset.x)
            , double(offset.y)
            , cl_long(par.p.formula.per.size())
            // ...
            };
    
          /* reference layout */
          cl_long ref_size = 0;
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            ref_size += config_host.ref_size[phase] = Zd[phase].size();
          }
          cl_long ref_bytes = 2 * ref_size * sizeof(double);
          config_host.ref_start[0] = 0;
          for (int phase = 1; phase < config_host.number_of_phases; ++phase)
          {
            config_host.ref_start[phase] = config_host.ref_start[phase - 1] + 2 * config_host.ref_size[phase - 1];
          }
    
          /* bla layout */
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            config_host.bla_size[phase] = Bd[phase].M;
            config_host.bla_levels[phase] = Bd[phase].L;
          }
          cl_long bla_start = 0;
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            for (int level = 0; level < config_host.bla_levels[phase]; ++level)
            {
              config_host.bla_start[phase][level] = bla_start;
              bla_start += Bd[phase].b[level].size();
            }
          }
          cl_long bla_bytes = bla_start * sizeof(struct blaR2_cl<double>);
    
          /* upload config */
          E(clEnqueueWriteBuffer(commands, config_device, CL_FALSE, 0, sizeof(config_host), &config_host, 0, 0, &ready));
    
          /* upload reference */
          ref_device = clCreateBuffer(context, CL_MEM_READ_ONLY, ref_bytes, 0, &err);
          if (! ref_device) { E(err); }
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            const cl_long start_bytes = config_host.ref_start[phase] * sizeof(double);
            const cl_long size_bytes = config_host.ref_size[phase] * 2 * sizeof(double);
            const void *ptr = &Zd[phase][0];
            cl_event done;
            E(clEnqueueWriteBuffer(commands, ref_device, CL_FALSE, start_bytes, size_bytes, ptr, 1, &ready, &done));
            E(clReleaseEvent(ready));
            ready = done;
          }
    
          /* upload bla */
          bla_device = clCreateBuffer(context, CL_MEM_READ_ONLY, bla_bytes, 0, &err);
          if (! bla_device) { E(err); }
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            for (int level = 0; level < config_host.bla_levels[phase]; ++level)
            {
              const cl_long start_bytes = config_host.bla_start[phase][level] * sizeof(blaR2_cl<double>);
              const cl_long size_bytes = Bd[phase].b[level].size() * sizeof(blaR2_cl<double>);
              const void *ptr = &Bd[phase].b[level][0];
              cl_event done;
              E(clEnqueueWriteBuffer(commands, bla_device, CL_FALSE, start_bytes, size_bytes, ptr, 1, &ready, &done));
              E(clReleaseEvent(ready));
              ready = done;
            }
          }
        }
        break;
      
      case nt_floatexp:
        {
          struct config_cl<floatexp> config_host =
            { sizeof(config_host)
            , nt_floatexp
            , height
            , width
            , par.p.image.subframes
            , frame
            , par.p.bailout.iterations
            , floatexp(par.p.bailout.escape_radius * par.p.bailout.escape_radius)
            , par.p.bailout.maximum_perturb_iterations
            , par.p.transform.exponential_map
            , { floatexp(par.transform.x[0][0]), floatexp(par.transform.x[0][1]), floatexp(par.transform.x[1][0]), floatexp(par.transform.x[1][1]) }
            , floatexp(4 / par.zoom / height)
            , floatexp(offset.x)
            , floatexp(offset.y)
            , cl_long(par.p.formula.per.size())
            // ...
            };
    
          /* reference layout */
          cl_long ref_size = 0;
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            ref_size += config_host.ref_size[phase] = Zfe[phase].size();
          }
          cl_long ref_bytes = 2 * ref_size * sizeof(floatexp);
          config_host.ref_start[0] = 0;
          for (int phase = 1; phase < config_host.number_of_phases; ++phase)
          {
            config_host.ref_start[phase] = config_host.ref_start[phase - 1] + 2 * config_host.ref_size[phase - 1];
          }
    
          /* bla layout */
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            config_host.bla_size[phase] = Bfe[phase].M;
            config_host.bla_levels[phase] = Bfe[phase].L;
          }
          cl_long bla_start = 0;
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            for (int level = 0; level < config_host.bla_levels[phase]; ++level)
            {
              config_host.bla_start[phase][level] = bla_start;
              bla_start += Bfe[phase].b[level].size();
            }
          }
          cl_long bla_bytes = bla_start * sizeof(struct blaR2_cl<floatexp>);
    
          /* upload config */
          E(clEnqueueWriteBuffer(commands, config_device, CL_FALSE, 0, sizeof(config_host), &config_host, 0, 0, &ready));
    
          /* upload reference */
          ref_device = clCreateBuffer(context, CL_MEM_READ_ONLY, ref_bytes, 0, &err);
          if (! ref_device) { E(err); }
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            const cl_long start_bytes = config_host.ref_start[phase] * sizeof(floatexp);
            const cl_long size_bytes = config_host.ref_size[phase] * 2 * sizeof(floatexp);
            const void *ptr = &Zfe[phase][0];
            cl_event done;
            E(clEnqueueWriteBuffer(commands, ref_device, CL_FALSE, start_bytes, size_bytes, ptr, 1, &ready, &done));
            E(clReleaseEvent(ready));
            ready = done;
          }
    
          /* upload bla */
          bla_device = clCreateBuffer(context, CL_MEM_READ_ONLY, bla_bytes, 0, &err);
          if (! bla_device) { E(err); }
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            for (int level = 0; level < config_host.bla_levels[phase]; ++level)
            {
              const cl_long start_bytes = config_host.bla_start[phase][level] * sizeof(blaR2_cl<floatexp>);
              const cl_long size_bytes = Bfe[phase].b[level].size() * sizeof(blaR2_cl<floatexp>);
              const void *ptr = &Bfe[phase].b[level][0];
              cl_event done;
              E(clEnqueueWriteBuffer(commands, bla_device, CL_FALSE, start_bytes, size_bytes, ptr, 1, &ready, &done));
              E(clReleaseEvent(ready));
              ready = done;
            }
          }
        }
        break;
      
      case nt_softfloat:
        {
          struct config_cl<softfloat> config_host =
            { sizeof(config_host)
            , nt_softfloat
            , height
            , width
            , par.p.image.subframes
            , frame
            , par.p.bailout.iterations
            , softfloat(par.p.bailout.escape_radius * par.p.bailout.escape_radius)
            , par.p.bailout.maximum_perturb_iterations
            , par.p.transform.exponential_map
            , { softfloat(par.transform.x[0][0]), softfloat(par.transform.x[0][1]), softfloat(par.transform.x[1][0]), softfloat(par.transform.x[1][1]) }
            , softfloat(4 / par.zoom / height)
            , softfloat(offset.x)
            , softfloat(offset.y)
            , cl_long(par.p.formula.per.size())
            // ...
            };
    
          /* reference layout */
          cl_long ref_size = 0;
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            ref_size += config_host.ref_size[phase] = Zsf[phase].size();
          }
          cl_long ref_bytes = 2 * ref_size * sizeof(double);
          config_host.ref_start[0] = 0;
          for (int phase = 1; phase < config_host.number_of_phases; ++phase)
          {
            config_host.ref_start[phase] = config_host.ref_start[phase - 1] + 2 * config_host.ref_size[phase - 1];
          }
    
          /* bla layout */
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            config_host.bla_size[phase] = Bsf[phase].M;
            config_host.bla_levels[phase] = Bsf[phase].L;
          }
          cl_long bla_start = 0;
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            for (int level = 0; level < config_host.bla_levels[phase]; ++level)
            {
              config_host.bla_start[phase][level] = bla_start;
              bla_start += Bsf[phase].b[level].size();
            }
          }
          cl_long bla_bytes = bla_start * sizeof(struct blaR2_cl<softfloat>);
    
          /* upload config */
          E(clEnqueueWriteBuffer(commands, config_device, CL_FALSE, 0, sizeof(config_host), &config_host, 0, 0, &ready));
    
          /* upload reference */
          ref_device = clCreateBuffer(context, CL_MEM_READ_ONLY, ref_bytes, 0, &err);
          if (! ref_device) { E(err); }
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            const cl_long start_bytes = config_host.ref_start[phase] * sizeof(softfloat);
            const cl_long size_bytes = config_host.ref_size[phase] * 2 * sizeof(softfloat);
            const void *ptr = &Zsf[phase][0];
            cl_event done;
            E(clEnqueueWriteBuffer(commands, ref_device, CL_FALSE, start_bytes, size_bytes, ptr, 1, &ready, &done));
            E(clReleaseEvent(ready));
            ready = done;
          }
    
          /* upload bla */
          bla_device = clCreateBuffer(context, CL_MEM_READ_ONLY, bla_bytes, 0, &err);
          if (! bla_device) { E(err); }
          for (int phase = 0; phase < config_host.number_of_phases; ++phase)
          {
            for (int level = 0; level < config_host.bla_levels[phase]; ++level)
            {
              const cl_long start_bytes = config_host.bla_start[phase][level] * sizeof(blaR2_cl<softfloat>);
              const cl_long size_bytes = Bsf[phase].b[level].size() * sizeof(blaR2_cl<softfloat>);
              const void *ptr = &Bsf[phase].b[level][0];
              cl_event done;
              E(clEnqueueWriteBuffer(commands, bla_device, CL_FALSE, start_bytes, size_bytes, ptr, 1, &ready, &done));
              E(clReleaseEvent(ready));
              ready = done;
            }
          }
        }
        break;

      default:
        {
          fprintf(stderr, "error: number type '%s' not supported\n", nt_string[nt_current]);
          abort();
        }
        break;
    }
    /* set up kernel arguments */
    E(clSetKernelArg(kernel, 0, sizeof(cl_mem), &config_device));
    E(clSetKernelArg(kernel, 1, sizeof(cl_mem), &ref_device));
    E(clSetKernelArg(kernel, 2, sizeof(cl_mem), &bla_device));
    E(clSetKernelArg(kernel, 3, sizeof(cl_mem), &rgb_device));
    E(clSetKernelArg(kernel, 4, sizeof(cl_mem), par.p.image.subframes == 1 ? &N0_device : nullptr));
    E(clSetKernelArg(kernel, 5, sizeof(cl_mem), par.p.image.subframes == 1 ? &N1_device : nullptr));
    E(clSetKernelArg(kernel, 6, sizeof(cl_mem), par.p.image.subframes == 1 ? &NF_device : nullptr));
    E(clSetKernelArg(kernel, 7, sizeof(cl_mem), par.p.image.subframes == 1 ? &T_device : nullptr));
    E(clSetKernelArg(kernel, 8, sizeof(cl_mem), par.p.image.subframes == 1 ? &DEX_device : nullptr));
    E(clSetKernelArg(kernel, 9, sizeof(cl_mem), par.p.image.subframes == 1 ? &DEY_device : nullptr));
    const cl_long tile_count =
      ((width + par.p.opencl.tile_width - 1) / par.p.opencl.tile_width) *
      ((height + par.p.opencl.tile_height - 1) / par.p.opencl.tile_height);
    cl_long tile = 0;
    for (cl_long y0 = 0; y0 < height; y0 += par.p.opencl.tile_height)
    {
      E(clSetKernelArg(kernel, 10, sizeof(cl_long), &y0));
      for (cl_long x0 = 0; x0 < width; x0 += par.p.opencl.tile_width)
      {
        E(clSetKernelArg(kernel, 11, sizeof(cl_long), &x0));
        progress[2 * count + 1] = tile++ / progress_t(tile_count);
        for (cl_long subframe = 0; subframe < par.p.image.subframes; ++subframe)
        {
          E(clSetKernelArg(kernel, 12, sizeof(cl_long), &subframe));
          progress[2 * count + 2] = subframe / progress_t(par.p.image.subframes);
          cl_event done;
          size_t global[2] = { (size_t) par.p.opencl.tile_height, (size_t) par.p.opencl.tile_width };
          E(clEnqueueNDRangeKernel(commands, kernel, 2, 0, global, 0, 1, &ready, &done));
          E(clReleaseEvent(ready));
          ready = done;
          if (! running)
          {
            break;
          }
          clFinish(commands);
        }
      }
    }
    if (running)
    {
      progress[2 * count + 1] = 1;
      /* synchronous transfer from device to host */
      E(clEnqueueReadBuffer(commands, rgb_device, CL_TRUE, 0, rgb_bytes, rgb_host, 1, &ready, 0));
      if (N0_device && out.N0)
      {
        E(clEnqueueReadBuffer(commands, N0_device, CL_TRUE, 0, raw_bytes, out.N0, 1, &ready, 0));
      }
      if (N1_device && out.N1)
      {
        E(clEnqueueReadBuffer(commands, N1_device, CL_TRUE, 0, raw_bytes, out.N1, 1, &ready, 0));
      }
      if (NF_device && out.NF)
      {
        E(clEnqueueReadBuffer(commands, NF_device, CL_TRUE, 0, raw_bytes, out.NF, 1, &ready, 0));
      }
      if (T_device && out.T)
      {
        E(clEnqueueReadBuffer(commands, T_device, CL_TRUE, 0, raw_bytes, out.T, 1, &ready, 0));
      }
      if (DEX_device && out.DEX)
      {
        E(clEnqueueReadBuffer(commands, DEX_device, CL_TRUE, 0, raw_bytes, out.DEX, 1, &ready, 0));
      }
      if (DEY_device && out.DEY)
      {
        E(clEnqueueReadBuffer(commands, DEY_device, CL_TRUE, 0, raw_bytes, out.DEY, 1, &ready, 0));
      }
    }
    E(clReleaseEvent(ready));
    if (running)
    {
      parallel2d(threads, 0, width, 32, 0, height, 32, running, [&](coord_t i, coord_t j) -> void
      {
        const size_t k = 3 * (j * width + i);
        out.setR(i, j, rgb_host[k + 0]);
        out.setG(i, j, rgb_host[k + 1]);
        out.setB(i, j, rgb_host[k + 2]);
      });
    }
    if (running)
    {
      if (par.p.render.zoom_out_sequence)
      {
        std::ostringstream s;
        s << par.p.render.filename << "." << std::setfill('0') << std::setw(8) << frame << ".exr";
        out.saveEXR(s.str(), par.p.image.subframes == 1 ? Channels_all : Channels_RGB, threads);
      }
      else
      {
        out.saveEXR(par.p.render.filename + ".exr", par.p.image.subframes == 1 ? Channels_all : Channels_RGB, threads);
      }
    }
    if (bla_device) { clReleaseMemObject(bla_device); bla_device = 0; }
    if (ref_device) { clReleaseMemObject(ref_device); ref_device = 0; }
  }

  if (rgb_host) { free(rgb_host); rgb_host = 0; }
  if (N0_device) { clReleaseMemObject(N0_device); N0_device = 0; }
  if (N1_device) { clReleaseMemObject(N1_device); N1_device = 0; }
  if (NF_device) { clReleaseMemObject(NF_device); NF_device = 0; }
  if (T_device) { clReleaseMemObject(T_device); T_device = 0; }
  if (DEX_device) { clReleaseMemObject(DEX_device); DEX_device = 0; }
  if (DEY_device) { clReleaseMemObject(DEY_device); DEY_device = 0; }
  if (rgb_device) { clReleaseMemObject(rgb_device); rgb_device = 0; }
  if (config_device) { clReleaseMemObject(config_device); config_device = 0; }

  /* cleanup */
  if (kernel) { /* FIXME */ }
  if (commands) { clReleaseCommandQueue(commands); commands = 0; }
  if (context) { clReleaseContext(context); context = 0; }
  if (source) { free(source); source = 0; }

  *ended = true;
}

int main(int argc, char **argv)
{
  if (! (argc == 2 || argc == 3))
  {
    std::cerr << version() << std::endl;
    std::cerr << "usage: " << argv[0] << " in.toml [outstem]" << std::endl;
    return 1;
  }

#ifdef HAVE_CLEW
  if (clewInit())
  {
    std::cerr << "error: clewInit() failed" << std::endl;
    return 1;
  }
#endif

  populate_number_type_wisdom();
  colours_init();

  param par;
  std::string default_filename = par.p.render.filename;
  if (argc > 1)
  {
    par.load_toml(argv[1]);
  }
  if (argc > 2)
  {
    par.p.render.filename = argv[2];
  }
  else if (par.p.render.filename == default_filename)
  {
    // FIXME remove extension?
  }

  map out((par.p.image.width + par.p.image.subsampling - 1) / par.p.image.subsampling, (par.p.image.height + par.p.image.subsampling - 1) / par.p.image.subsampling, par.p.bailout.iterations);

  const count_t count = par.p.formula.per.size();
  std::vector<progress_t> progress;
  progress.resize(2 * count + 3);
  for (count_t i = 0; i < count; ++i)
  {
    progress[i] = 0;
  }
  bool running = true;
  bool ended = false;

  std::thread bg(opencl_thread, std::ref(out), std::ref(par), &progress[0], &running, &ended);
  while (! ended)
  {
    for (count_t ms = 0; ms < 500 && ! ended; ++ms)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    std::ostringstream s;
    s << "Frame["         << std::setw(3) << int(progress[0] * 100) << "%] ";
    progress_t r = 0;
    for (count_t i = 0; i < count; ++i)
    {
      r += progress[1 + i];
    }
    s << "Reference[" << std::setw(3) << int(r * 100 / count) << "%] ";
    progress_t a = 0;
    for (count_t i = 0; i < count; ++i)
    {
      a += progress[1 + count + i];
    }
    s << "Approximation[" << std::setw(3) << int(a * 100 / count) << "%] ";
    s << "Tile["          << std::setw(3) << int(progress[2 * count + 1] * 100) << "%] ";
    s << "Subframe["      << std::setw(3) << int(progress[2 * count + 2] * 100) << "%] ";
    s << "\r";
    std::cerr << s.str();
  }
  std::cerr << "\n";
  bg.join();

  return 0;
}
