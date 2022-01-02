// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <chrono>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <thread>

#ifdef HAVE_OMP
#include <omp.h>
#else
int omp_get_num_procs()
{
  return 1;
}
#endif

#include "colour.h"
#include "display_cpu.h"
#include "engine.h"
#include "floatexp.h"
#include "formula.h"
#include "map.h"
#include "param.h"
#include "stats.h"
#include "types.h"
#include "version.h"

void cli_thread(display_cpu &dsp, map &out, stats &sta, param &par, const formula *form, progress_t *progress, bool *running, bool *ended)
{
  int threads = omp_get_num_procs();
  floatexp Zoom = par.zoom;
  floatexp ZoomedOut = 1 / 65536.0;
  count_t nframes = (Zoom / ZoomedOut).exp + 1;
  count_t start_frame = par.p.render.start_frame;
  count_t end_frame = par.p.render.frame_count == 0 ? nframes : par.p.render.start_frame + par.p.render.frame_count;
  nframes = end_frame - start_frame;
  if (nframes == 0)
  {
    nframes = 1;
  }
  reset(sta);
  if (par.p.render.zoom_out_sequence)
  {
    for (count_t frame = start_frame; frame < end_frame; ++frame)
    {
      par.zoom.exp = Zoom.exp - frame;
      progress[0] = (frame - start_frame) / progress_t(nframes);
      bool ref_ended = false;
      reference_thread(sta, form, par, &progress[1], running, &ref_ended);
      dsp.clear();
      for (count_t subframe = 0; subframe < par.p.image.subframes; subframe++)
      {
        progress[3] = subframe / progress_t(par.p.image.subframes);
        bool sub_ended = false;
        subframe_thread(out, sta, form, par, subframe, &progress[4], running, &sub_ended);
        if (! *running)
        {
          break;
        }
        dsp.accumulate(out);
      }
      if (! *running)
      {
        break;
      }
      dsp.get_rgb(out);
      std::ostringstream s;
      s << par.p.render.filename << "." << std::setfill('0') << std::setw(8) << frame << ".exr";
      out.saveEXR(s.str(), Channels_RGB, threads);
    }
  }
  else
  {
    progress[0] = 0;
    bool ref_ended = false;
    reference_thread(sta, form, par, &progress[1], running, &ref_ended);
    dsp.clear();
    for (count_t subframe = 0; subframe < par.p.image.subframes; subframe++)
    {
      progress[3] = subframe / progress_t(par.p.image.subframes);
      bool sub_ended = false;
      subframe_thread(out, sta, form, par, subframe, &progress[4], running, &sub_ended);
      if (! *running)
      {
        break;
      }
      dsp.accumulate(out);
    }
    progress[3] = 1;
    if (running)
    {
      dsp.get_rgb(out);
      out.saveEXR(par.p.render.filename + ".exr", Channels_RGB, threads);
    }
  }
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

  formulas_init();
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

  formula *form = formulas[par.p.formula_id];
  colour *clr = colours[par.p.colour_id];

  display_cpu dsp(clr);
  dsp.resize(out.width, out.height);

  stats sta;

  progress_t progress[5] = { 0, 0, 0, 0, 0 };
  bool running = true;
  bool ended = false;

  std::thread bg(cli_thread, std::ref(dsp), std::ref(out), std::ref(sta), std::ref(par), form, &progress[0], &running, &ended);
  while (! ended)
  {
    for (count_t ms = 0; ms < 500 && ! ended; ++ms)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    std::cerr
      << "Frame["         << std::setw(3) << int(progress[0] * 100) << "%] "
      << "Reference["     << std::setw(3) << int(progress[1] * 100) << "%] "
      << "Approximation[" << std::setw(3) << int(progress[2] * 100) << "%] "
      << "Subframe["      << std::setw(3) << int(progress[3] * 100) << "%] "
      << "Pixels["        << std::setw(3) << int(progress[4] * 100) << "%] "
      << "\r";
  }
  std::cerr << "\n";
  bg.join();

  return 0;
}
