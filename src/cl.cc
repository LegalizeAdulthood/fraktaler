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
#include <CL/cl.h>

#include "bla.h"
#include "colour.h"
#include "display_cpu.h"
#include "engine.h"
#include "floatexp.h"
#include "map.h"
#include "parallel.h"
#include "param.h"
#include "stats.h"
#include "types.h"
#include "version.h"

#include "cl-pre.h"
#include "cl-post.h"

extern std::vector<std::vector<complex<double>>> Zd;
extern std::vector<blasR2<double>> Bd;

std::string hybrid_perturb(const std::vector<phybrid1> &per)
{
  std::ostringstream s;
  s << "{\n";
  s << "  struct complex Z = { ref[config->ref_start[phase] + 2 * m], ref[config->ref_start[phase] + 2 * m + 1] };\n";
  s << "  double X = Z.x;\n";
  s << "  double Y = Z.y;\n";
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
      s << "      x = dual_diffabs_double_dual(X, x);\n";
      s << "      W.x = dual_abs_dual(W.x);\n";
      s << "      B.x = fabs(B.x);\n";
    }
    if (per[k].abs_y)
    {
      s << "      y = dual_diffabs_double_dual(Y, y);\n";
      s << "      W.y = dual_abs_dual(W.y);\n";
      s << "      B.y = fabs(B.y);\n";
    }
    if (per[k].neg_x)
    {
      s << "      x = dual_neg_dual(x);\n";
      s << "      W.x = dual_neg_dual(W.x);\n";
      s << "      B.x = -B.x;\n";
    }
    if (per[k].neg_y)
    {
      s << "      y = dual_neg_dual(y);\n";
      s << "      W.y = dual_neg_dual(W.y);\n";
      s << "      B.y = -B.y;\n";
    }
    s << "      struct complexdual P = { x, y };\n";
    s << "      struct complexdual S = { { 0, { 0, 0 } }, { 0, { 0, 0 } } };\n";
    s << "      struct complex Bn[" << per[k].power << "];\n";
    s << "      Bn[0].x = 1; Bn[0].y = 0;\n";
    for (int i = 1; i < per[k].power; ++i)
    {
      s << "      Bn["  << i << "] = complex_mul_complex_complex(Bn[" << (i - 1) << "], B);\n";
    }
    s << "      struct complexdual Wi = S; Wi.x.x = 1;\n";
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

struct mat2_cl
{
  cl_double a, b, c, d;
};

struct blaR2_cl
{
  struct mat2_cl A, B;
  cl_double r2;
  cl_long l;
};

struct config_cl
{
  /* shape */
  cl_long height;
  cl_long width;
  cl_long subframes;
  cl_long frame;
  /* bailout */
  cl_long Iterations;
  cl_double ER2;
  cl_long PerturbIterations;
  /* transform */
  cl_long transform_exponential_map;
  struct mat2_cl transform_K;
  cl_double pixel_spacing;
  cl_double offset_x;
  cl_double offset_y;
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
  if (dvecsize == 0)
  {
    fprintf(stderr, "error: device does not support double precision\n");
    abort();
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
  const char *options = "";
  err = clBuildProgram(program, 1, &device_id[device_index], options, 0, 0);
  if (err != CL_SUCCESS)
  {
    char *error_log = (char *) malloc(1000000);
    error_log[0] = 0;
    E(clGetProgramBuildInfo(program, device_id[device_index], CL_PROGRAM_BUILD_LOG, 1000000, &error_log[0], 0));
    fprintf(stderr, "error: compile failed:\n%s\n", error_log);
    fprintf(stderr, "error: correspinding source:\n%s\n", source);
    E(err);
  }
  cl_kernel kernel = clCreateKernel(program, "fraktaler3", &err);
  if (! kernel) { E(err); }

  coord_t width = par.p.image.width / par.p.image.subsampling;
  coord_t height = par.p.image.height / par.p.image.subsampling;
  size_t grey_bytes = sizeof(float) * width * height;
  float *grey_host = (float *) malloc(grey_bytes);
  if (! grey_host)
  {
    fprintf(stderr, "error: out of memory\n");
    abort();
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
  if (par.p.render.zoom_out_sequence)
  {
#if 0
    for (count_t frame = start_frame; frame < end_frame; ++frame)
    {
      par.zoom = Zoom / pow(floatexp(par.p.render.zoom_out_factor), frame);
      progress[0] = (frame - start_frame) / progress_t(nframes);
      bool ref_ended = false;
      reference_thread(sta, par, false, &progress[1], running, &ref_ended);

      const float zero = 0;
      cl_event cleared;
      E(clEnqueueFillBuffer(commands, rgb, &zero, sizeof(zero), 0, rgb_bytes, 0, 0, &cleared));
      cl_event ready = cleared;
      for (int subframe = 0; subframe < par.p.image.subframes; ++subframe)
      {
        progress[2 * count + 1] = subframe / progress_t(par.p.image.subframes);
        /* prepare configuration data */
        struct config host_config =
          { height
          , width
          , subframe
          , subframes
          , ...
          };
        /* synchronous transfer from host to device */
        E(clEnqueueWriteBuffer(commands, config,  CL_TRUE, 0, sizeof(struct config), &host_config, 1, &ready, 0));
        cl_event done;
        size_t global[2] = { height, width };
        E(clEnqueueNDRangeKernel(commands, kernel, 2, 0, global, 0, 0, 0, &done));
        ready = done;
        if (! *running)
        {
          break;
        }
      }
      if (running)
      {
        progress[2 * count + 1] = 1;
        /* synchronous transfer from device to host */
        E(clEnqueueReadBuffer(commands, rgb, CL_TRUE, 0, sizeof(float) * par.p.image.width * par.p.image.height * 3, host_rgb, 1, &ready, 0));
      }
      if (running)
      {
        parallel2d(threads, 0, 32, width, 0, 32, height, running, [&](coord_t i, coord_t j) -> void
        {
          const float v = grey_host[j * width + i];
          out.setR(i, j, v);
          out.setG(i, j, v);
          out.setB(i, j, v);
        });
      }
      if (running)
      {
        std::ostringstream s;
        s << par.p.render.filename << "." << std::setfill('0') << std::setw(8) << frame << ".exr";
        out.saveEXR(s.str(), /*par.p.image.subframes == 1 ? Channels_all : */Channels_RGB, threads);
      }
    }
#endif
  }
  else
  {
    progress[0] = 0;
    bool ref_ended = false;
    {
      stats sta;
      reference_thread(sta, par, false, &progress[1], running, &ref_ended);
    }
    complex<mpreal> offset;
    offset.x.set_prec(par.center.x.get_prec());
    offset.y.set_prec(par.center.y.get_prec());
    offset = par.center - par.reference;

    struct config_cl config_host =
      { height
      , width
      , par.p.image.subframes
      , 0
      , par.p.bailout.iterations
      , par.p.bailout.escape_radius * par.p.bailout.escape_radius
      , par.p.bailout.maximum_perturb_iterations
      , par.p.transform.exponential_map
      , { par.transform.x[0][0], par.transform.x[0][1], par.transform.x[1][0], par.transform.x[1][1] }
      , double(4 / Zoom / height)
      , double(offset.x)
      , double(offset.y)
      , cl_long(par.p.formula.per.size())
      // ...
      };

std::cerr << (config_host.width) << std::endl;
std::cerr << (config_host.height) << std::endl;
std::cerr << (config_host.subframes) << std::endl;
std::cerr << (config_host.frame) << std::endl;
std::cerr << (config_host.Iterations) << std::endl;
std::cerr << (config_host.ER2) << std::endl;
std::cerr << (config_host.PerturbIterations) << std::endl;
std::cerr << (config_host.transform_exponential_map) << std::endl;
std::cerr << (config_host.transform_K.a) << std::endl;
std::cerr << (config_host.transform_K.b) << std::endl;
std::cerr << (config_host.transform_K.c) << std::endl;
std::cerr << (config_host.transform_K.d) << std::endl;
std::cerr << (config_host.pixel_spacing) << std::endl;
std::cerr << (config_host.offset_x) << std::endl;
std::cerr << (config_host.offset_y) << std::endl;
std::cerr << (config_host.number_of_phases) << std::endl;

    /* reference layout */
    cl_long ref_size = 0;
    for (int phase = 0; phase < config_host.number_of_phases; ++phase)
    {
std::cerr << Zd[phase].size() << std::endl;
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
std::cerr << Bd[phase].M << " " << Bd[phase].L << std::endl;
      config_host.bla_size[phase] = Bd[phase].M;
      config_host.bla_levels[phase] = Bd[phase].L;
    }
    cl_long bla_start = 0;
    for (int phase = 0; phase < config_host.number_of_phases; ++phase)
    {
      for (int level = 0; level < config_host.bla_levels[phase]; ++level)
      {
std::cerr << bla_start << " ";
        config_host.bla_start[phase][level] = bla_start;
        bla_start += Bd[phase].b[level].size();
      }
std::cerr << std::endl;
    }
    cl_long bla_bytes = bla_start * sizeof(struct blaR2_cl);

    cl_event ready;

    /* upload config */
    cl_mem config_device = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(struct config_cl), 0, &err);
    if (! config_device) { E(err); }
    E(clEnqueueWriteBuffer(commands, config_device, CL_FALSE, 0, sizeof(struct config_cl), &config_host, 0, 0, &ready));

    /* upload reference */
    cl_mem ref_device = clCreateBuffer(context, CL_MEM_READ_ONLY, ref_bytes, 0, &err);
    if (! ref_device) { E(err); }
    for (int phase = 0; phase < config_host.number_of_phases; ++phase)
    {
      fprintf(stderr, "ref upload phase %d\n", phase);
      const cl_long start_bytes = config_host.ref_start[phase] * sizeof(double);
      const cl_long size_bytes = config_host.ref_size[phase] * 2 * sizeof(double);
      const void *ptr = &Zd[phase][0];
      cl_event done;
      E(clEnqueueWriteBuffer(commands, ref_device, CL_FALSE, start_bytes, size_bytes, ptr, 1, &ready, &done));
      ready = done;
    }

    /* upload bla */
    cl_mem bla_device = clCreateBuffer(context, CL_MEM_READ_ONLY, bla_bytes, 0, &err);
    if (! bla_device) { E(err); }
    for (int phase = 0; phase < config_host.number_of_phases; ++phase)
    {
      for (int level = 0; level < config_host.bla_levels[phase]; ++level)
      {
        fprintf(stderr, "bla upload phase %d level %d\n", phase, level);
        const cl_long start_bytes = config_host.bla_start[phase][level] * sizeof(blaR2_cl);
        const cl_long size_bytes = Bd[phase].b[level].size() * sizeof(blaR2_cl);
        const void *ptr = &Bd[phase].b[level][0];
        cl_event done;
        E(clEnqueueWriteBuffer(commands, bla_device, CL_FALSE, start_bytes, size_bytes, ptr, 1, &ready, &done));
        ready = done;
      }
    }

    cl_mem grey_device = clCreateBuffer(context, CL_MEM_READ_WRITE, grey_bytes, 0, &err);
    if (! grey_device) { E(err); }

    /* set up kernel arguments */
    E(clSetKernelArg(kernel, 0, sizeof(cl_mem), &config_device));
    E(clSetKernelArg(kernel, 1, sizeof(cl_mem), &ref_device));
    E(clSetKernelArg(kernel, 2, sizeof(cl_mem), &bla_device));
    E(clSetKernelArg(kernel, 3, sizeof(cl_mem), &grey_device));
    for (cl_long subframe = 0; subframe < par.p.image.subframes; ++subframe)
    {
      fprintf(stderr, "subframe %d\n", (int) subframe);
      E(clSetKernelArg(kernel, 4, sizeof(cl_long), &subframe));
      progress[2 * count + 1] = subframe / progress_t(par.p.image.subframes);
      cl_event done;
      size_t global[2] = { (size_t) height, (size_t) width };
      E(clEnqueueNDRangeKernel(commands, kernel, 2, 0, global, 0, 1, &ready, &done));
      ready = done;
      if (! running)
      {
        break;
      }
      clFinish(commands);
    }
    if (running)
    {
      fprintf(stderr, "reading\n");
      progress[2 * count + 1] = 1;
      /* synchronous transfer from device to host */
      E(clEnqueueReadBuffer(commands, grey_device, CL_TRUE, 0, grey_bytes, grey_host, 1, &ready, 0));
    }
    if (running)
    {
      fprintf(stderr, "setting\n");
      parallel2d(threads, 0, width, 32, 0, height, 32, running, [&](coord_t i, coord_t j) -> void
      {
        const float v = grey_host[j * width + i];
        out.setR(i, j, v);
        out.setG(i, j, v);
        out.setB(i, j, v);
      });
    }
    if (running)
    {
      fprintf(stderr, "saving\n");
      out.saveEXR(par.p.render.filename + ".exr", /*par.p.image.subframes == 1 ? Channels_all : */Channels_RGB, threads);
    }
    if (grey_host) { free(grey_host); grey_host = 0; }
    if (grey_device) { clReleaseMemObject(grey_device); grey_device = 0; }
    if (bla_device) { clReleaseMemObject(bla_device); bla_device = 0; }
    if (ref_device) { clReleaseMemObject(ref_device); ref_device = 0; }
    if (config_device) { clReleaseMemObject(config_device); config_device = 0; }
  }

  /* cleanup */
  if (kernel) { /* ... */ }
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
    s << "Subframe["      << std::setw(3) << int(progress[2 * count + 1] * 100) << "%] ";
    s << "Pixels["        << std::setw(3) << int(progress[2 * count + 2] * 100) << "%] ";
    s << "\r";
    std::cerr << s.str();
  }
  std::cerr << "\n";
  bg.join();

  return 0;
}
