// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <iomanip>
#include <thread>

#include "image_raw.h"
#include "image_rgb.h"
#include "param.h"
#include "render.h"
#include "wisdom.h"

struct batch_hooks : public hooks
{
  image_rgb *img_rgb;
  image_raw *img_raw;
  const param &par;
  coord_t frame;
  int threads;
  batch_hooks(image_rgb *img_rgb, image_raw *img_raw, const param &par, coord_t frame, int threads)
  : img_rgb(img_rgb)
  , img_raw(img_raw)
  , par(par)
  , frame(frame)
  , threads(threads)
  {
  }
  virtual void start()
  {
    if (img_rgb)
    {
      img_rgb->clear();
    }
    if (img_raw)
    {
      img_raw->clear();
    }
  }
  virtual void stop()
  {
    std::ostringstream s;
    if (frame >= 0)
    {
      s << par.p.render.filename << "." << std::setfill('0') << std::setw(8) << (frame++) << ".exr";
    }
    else
    {
      s << par.p.render.filename << ".exr";
    }
    if (img_rgb)
    {
      // normalize
      image_rgb(*img_rgb).save_exr(s.str(), threads, par.to_string() /* , metakf */);
    }
    if (img_raw)
    {
      img_raw->save_exr(s.str(), Channels_all, par.p.bailout.iterations, threads, par.to_string() /* , metakf */);
    }
  }
  virtual void tile(int platform, int device, int x, int y, int subframe, const struct tile *data)
  {
    (void) platform;
    (void) device;
    (void) subframe;
    if (img_rgb)
    {
      img_rgb->blit(x, y, data);
    }
    if (img_raw)
    {
      img_raw->blit(x, y, data);
    }
  }
};

void batch_thread(const param &par0, progress_t *progress, volatile bool *running, volatile bool *ended)
{
  param par = par0;
  if (par.p.image.subframes <= 0)
  {
    *ended = true;
    return;
  }
  image_rgb *img_rgb = par.p.image.subframes != 1 ? new image_rgb(par.p.image.width / par.p.image.subsampling, par.p.image.height / par.p.image.subsampling) : nullptr;
  image_raw *img_raw = par.p.image.subframes == 1 ? new image_raw(par.p.image.width / par.p.image.subsampling, par.p.image.height / par.p.image.subsampling, Channels_all) : nullptr;
  if ((!!img_rgb) == (!!img_raw))
  {
    if (img_rgb)
    {
      delete img_rgb;
    }
    if (img_raw)
    {
      delete img_raw;
    }
    *ended = true;
    return;
  }
  int threads = std::thread::hardware_concurrency();
  floatexp Zoom = par.zoom;
  floatexp ZoomedOut = 1 / 65536.0;
  count_t nframes = std::ceil(double(log(Zoom / ZoomedOut) / log(par.p.render.zoom_out_factor)));
  count_t start_frame = par.p.render.start_frame;
  count_t end_frame = par.p.render.frame_count == 0 ? nframes : par.p.render.start_frame + par.p.render.frame_count;
  nframes = end_frame - start_frame;
  if (nframes == 0)
  {
    nframes = 1;
  }
  end_frame = start_frame + nframes;
  if (! par.p.render.zoom_out_sequence)
  {
    start_frame = 0;
    end_frame = 1;
    nframes = 1;
  }
  const std::set<number_type> available = { nt_float, nt_double, nt_longdouble, nt_floatexp, nt_softfloat
#ifdef HAVE_FLOAT128
  , nt_float128
#endif
  };
  for (count_t frame = start_frame; frame < end_frame; ++frame)
  {
    par.zoom = Zoom / pow(floatexp(par.p.render.zoom_out_factor), frame);
    batch_hooks h(img_rgb, img_raw, par, frame, threads);
    progress[0] = (frame - start_frame) / progress_t(nframes);
    count_t pixel_spacing_exp, pixel_precision_exp;
    get_required_precision(par, pixel_spacing_exp, pixel_precision_exp);
    auto l = wisdom_lookup(wdom, available, pixel_spacing_exp, pixel_precision_exp);
    render(l, par, &h, &progress[1], running); // FIXME frame
    if (! *running)
    {
      break;
    }
  }
  *ended = true;
}

int batch(int verbosity, const param &par)
{
  const count_t count = par.p.formula.per.size();
  std::vector<progress_t> progress;
  progress.resize(2 * count + 2);
  for (count_t i = 0; i < count; ++i)
  {
    progress[i] = 0;
  }
  volatile bool running = true;
  volatile bool ended = false;
  std::thread bg(batch_thread, std::cref(par), &progress[0], &running, &ended);
  while (! ended)
  {
    for (count_t ms = 0; ms < 500 && ! ended; ++ms)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (verbosity > 0)
    {
      std::ostringstream s;
      s << "Frame["         << std::setw(3) << int(progress[0] * 100) << "%] ";
      progress_t r = 0;
      for (count_t i = 0; i < count; ++i)
      {
        r += progress[1 + i];
      }
      s << "Ref[" << std::setw(3) << int(r * 100 / count) << "%] ";
      progress_t a = 0;
      for (count_t i = 0; i < count; ++i)
      {
        a += progress[1 + count + i];
      }
      s << "BLA[" << std::setw(3) << int(a * 100 / count) << "%] ";
      s << "Tile[" << std::setw(3) << int(progress[2 * count + 1] * 100) << "%] ";
      s << "\r";
      std::cerr << s.str();
    }
  }
  if (verbosity > 0)
  {
    std::cerr << "\n";
  }
  bg.join();
  return 0;
}
