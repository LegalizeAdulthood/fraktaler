// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
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
  floatexp Zoom = par.Zoom;
  floatexp ZoomedOut = 1 / 65536.0;
  count_t nframes = (Zoom / ZoomedOut).exp + 1;
  reset(sta);
  if (par.ZoomOutSequence)
  {
    for (count_t frame = 0; frame < nframes; ++frame, --par.Zoom.exp)
    {
      progress[0] = frame / progress_t(nframes);
      bool ref_ended = false;
      reference_thread(sta, form, par, &progress[1], running, &ref_ended);
      dsp.clear();
      for (count_t subframe = 0; subframe < par.MaxSubframes; subframe++)
      {
        progress[3] = subframe / progress_t(par.MaxSubframes);
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
      s << par.Stem << "_" << std::setfill('0') << std::setw(8) << (frame++) << ".exr";
      out.saveEXR(s.str(), par.Channels, threads);
    }
  }
  else
  {
    progress[0] = 0;
    bool ref_ended = false;
    reference_thread(sta, form, par, &progress[1], running, &ref_ended);
    dsp.clear();
    for (count_t subframe = 0; subframe < par.MaxSubframes; subframe++)
    {
      progress[3] = subframe / progress_t(par.MaxSubframes);
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
      out.saveEXR(par.Stem, par.Channels, threads);
    }
  }
  *ended = true;
}

int main(int argc, char **argv)
{
  using std::isnan;
  using std::isinf;
  using std::log;
  using std::max;
  using std::min;

  const coord_t Width = 1920;
  const coord_t Height = 1080;
  const bool ExponentialMap = false;
  const bool ZoomOutSequence = false;

#if 1
  if (argc != 7)
  {
    std::cerr << version() << std::endl;
    std::cerr << "usage: " << argv[0] << " re im zoom angle iters out.exr" << std::endl;
    return 1;
  }
  const char *Re = argv[1];
  const char *Im = argv[2];
  floatexp Zoom = floatexp(mpreal(argv[3], 53));
  double Angle = atof(argv[4]);
  floatexp ZoomPrec = Zoom;
  count_t Iterations = atoll(argv[5]);
  count_t MaxRefIters = Iterations;
  count_t MaxPtbIters = 1024;
  count_t ReferencePeriod = 0;
  const std::string Stem = argv[6];
#else
#if 1
  // Dinkydau - Flake
  const char *Re = "-1.999966194450370304184346885063505796755312415407248515117619229448015842423426843813761297788689138122870464065609498643538105757447721664856724960928039200953321766548438990841";
  const char *Im = "3.001382436790938324072497303977592498734683119077333527017425728012047497561482358118564729928841407551922418650497818162547852894554815712214577268132087243118270209466408349072e-34";
  floatexp Zoom (e10(2.5620330788506154104770818136626, 157));
  floatexp ZoomPrec (e10(1, 171));
  double Angle = 0;
  count_t Iterations = 1100100;
  count_t MaxPtbIters = 1000;
  count_t ReferencePeriod = 7884;
  count_t MaxRefIters = ReferencePeriod + 1;
  const std::string Stem = argc > 1 ? argv[1] : "out.exr";
#else
#if 1
  // Dinkydau - Evolution of Trees
  const char *Re = "-0.74962449737876168207862620426836138684529527812364025754481571424558286479672801698750203976779135608148687747196595174858125388297577788573082753210469806739721377901297451172624762450940529814713048377873612297220313016947039287469999999999999999";
  const char *Im = "0.03427010874046016945172951545749474868051534672439111853174854738370595772935930171842999778713222794215453436628998200591914208207632674780780978496784843807690510401829301309069542198413574392425885166367231192087416338497065495509999999999999999";
  floatexp Zoom (7.58065474756E227);
  floatexp ZoomPrec = Zoom;
  double Angle = 0;
  count_t Iterations = 1200000;
  count_t MaxPtbIters = 8000;
  count_t ReferencePeriod = 567135;
  count_t MaxRefIters = ReferencePeriod + 1;
  const std::string Stem = argc > 1 ? argv[1] : "out.exr";
#else
  // Fractal Universe - Hard Location
  const char *Re = "-1.74992248092759992827133368754228945303043302447370334500650852139592486065065408129935547375121997659867849111435922542786389338654238247523763456e+00";
  const char *Im = "-9.59502198314327569948975707202650233401883670299418141500240641361234506320676962536124684582340235944852850785763764700482870741065313152008753180e-13";
  floatexp Zoom = 6.49E137;
  floatexp ZoomPrec = 6.367062888622018e138;
  double Angle = 0;
  count_t Iterations = 2100100100;
  count_t MaxPtbIters = 2000;
  count_t ReferencePeriod = 43292334;
  count_t MaxRefIters = ReferencePeriod + 1;
  const std::string Stem = argc > 1 ? argv[1] : "out.exr";
#endif
#endif
#endif

  formulas_init();
  colours_init();

  param par;
  home(par);
  par.Zoom = Zoom;
  par.Iterations = Iterations;
  par.ReferencePeriod = ReferencePeriod;
  par.MaxRefIters = MaxRefIters;
  par.MaxPtbIters = MaxPtbIters;
  par.ExponentialMap = ExponentialMap;
  par.ZoomOutSequence = ZoomOutSequence;
  par.Channels = Channels_default;
  par.Stem = Stem;
  par.Width = Width;
  par.Height = Height;
  par.MaxSubframes = 16;
  par.K = rotation(Angle);

  floatexp pixel_spacing = 4 / ZoomPrec / par.Height;
  mpfr_prec_t prec = 24 - pixel_spacing.exp;
  par.C.x.set_prec(prec);
  par.C.y.set_prec(prec);
  par.C.x = Re;
  par.C.y = Im;

  map out(par.Width, par.Height, par.Iterations);

  formula *form = formulas[0];
  colour *clr = colours[0];

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
