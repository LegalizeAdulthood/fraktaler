// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <chrono>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <thread>

#include <half.h>
#include <mpfr.h>
#include <omp.h>

#include "bla.h"
#include "complex.h"
#include "floatexp.h"
#include "main.h"
#include "map.h"
#include "param.h"
#include "stats.h"
#include "types.h"

#if 0
static inline double hypot2(double x, double y)
{
  return x * x + y * y;
}

void main_thread(map &out, stats &sta, const param &par, progress_t *progress, bool *running, bool *ended)
{
  int threads = omp_get_num_procs();
  floatexp Zoom = par.Zoom;
  floatexp ZoomedOut = 1 / 65536.0;
  progress_t nframes = (Zoom / ZoomedOut).exp + 1;
  reset(sta);
  if (par.ZoomOutSequence)
  {

    count_t frame = 0;
    count_t M = 0;
    complex<float> *Zf = nullptr;
    if (Zoom > e10(1, 30))
    {
      complex<double> *Zd = nullptr;
      if (Zoom > e10(1, 300))
      {
        complex<long double> *Zld = nullptr;
        if (Zoom > e10(1, 4900))
        {
          // calculate reference
          complex<floatexp> *Zfe = new complex<floatexp>[par.MaxRefIters];
          M = reference(Zfe, par.MaxRefIters, par.Cx, par.Cy, &progress[0], running);
          for (; Zoom > e10(1, 4900); Zoom >>= 1)
          {
            progress[1] = frame / nframes;
            progress[2] = 0;
            progress[3] = 0;
            render(out, sta, par, Zoom, M, Zfe, &progress[2], running);
            std::ostringstream s;
            s << par.Stem << "_" << std::setfill('0') << std::setw(8) << (frame++) << ".exr";
            out.saveEXR(s.str(), par.Channels, threads);
          }
          Zld = new complex<long double>[M];
          #pragma omp parallel for
          for (count_t m = 0; m < M; ++m)
          {
            complex<floatexp> Z = Zfe[m];
            Zld[m] = complex<long double>((long double)(Z.x), (long double)(Z.y));
          }
          delete[] Zfe;
        }
        else
        {
          // calculate reference
          complex<long double> *Zld = new complex<long double>[par.MaxRefIters];
          M = reference(Zld, par.MaxRefIters, par.Cx, par.Cy, &progress[0], running);
        }
        for (; Zoom > e10(1, 300); Zoom >>= 1)
        {
          progress[1] = frame / nframes;
          progress[2] = 0;
          progress[3] = 0;
          render(out, sta, par, (long double)(Zoom), M, Zld, &progress[2], running);
          std::ostringstream s;
          s << par.Stem << "_" << std::setfill('0') << std::setw(8) << (frame++) << ".exr";
          out.saveEXR(s.str(), par.Channels, threads);
        }
        Zd = new complex<double>[M];
        #pragma omp parallel for
        for (count_t m = 0; m < M; ++m)
        {
          complex<long double> Z = Zld[m];
          Zd[m] = complex<double>(double(Z.x), double(Z.y));
        }
        delete[] Zld;
      }
      else
      {
        // calculate reference
        complex<double> *Zd = new complex<double>[par.MaxRefIters];
        M = reference(Zd, par.MaxRefIters, par.Cx, par.Cy, &progress[0], running);
      }
      for (; Zoom > e10(1, 30); Zoom >>= 1)
      {
        progress[1] = frame / nframes;
        progress[2] = 0;
        progress[3] = 0;
        render(out, sta, par, double(Zoom), M, Zd, &progress[2], running);
        std::ostringstream s;
        s << par.Stem << "_" << std::setfill('0') << std::setw(8) << (frame++) << ".exr";
        out.saveEXR(s.str(), par.Channels, threads);
      }
      Zf = new complex<float>[M];
      #pragma omp parallel for
      for (count_t m = 0; m < M; ++m)
      {
        const complex<double> Z = Zd[m];
        Zf[m] = complex<float>(float(Z.x), float(Z.y));
      }
      delete[] Zd;
    }
    else
    {
      // calculate reference
      complex<float> *Zf = new complex<float>[par.MaxRefIters];
      M = reference(Zf, par.MaxRefIters, par.Cx, par.Cy, &progress[0], running);
    }
    for (; Zoom > ZoomedOut; Zoom >>= 1)
    {
      progress[1] = frame / nframes;
      progress[2] = 0;
      progress[3] = 0;
      render(out, sta, par, float(Zoom), M, Zf, &progress[2], running);
      std::ostringstream s;
      s << par.Stem << "_" << std::setfill('0') << std::setw(8) << (frame++) << ".exr";
      out.saveEXR(s.str(), par.Channels, threads);
    }
    delete[] Zf;

  }
  else
  {

    if (Zoom > e10(1, 4900))
    {
      complex<floatexp> *Zfe = new complex<floatexp>[par.MaxRefIters];
      const count_t M = reference(Zfe, par.MaxRefIters, par.Cx, par.Cy, &progress[0], running);
      progress[1] = 1;
      render(out, sta, par, Zoom, M, Zfe, &progress[2], running);
      delete[] Zfe;
    }
    else if (Zoom > e10(1, 300))
    {
      complex<long double> *Zld = new complex<long double>[par.MaxRefIters];
      const count_t M = reference(Zld, par.MaxRefIters, par.Cx, par.Cy, &progress[0], running);
      progress[1] = 1;
      render(out, sta, par, (long double)(Zoom), M, Zld, &progress[2], running);
      delete[] Zld;
    }
    else if (Zoom > e10(1, 30))
    {
      complex<double> *Zd = new complex<double>[par.MaxRefIters];
      const count_t M = reference(Zd, par.MaxRefIters, par.Cx, par.Cy, &progress[0], running);
      progress[1] = 1;
      render(out, sta, par, double(Zoom), M, Zd, &progress[2], running);
      delete[] Zd;
    }
    else
    {
      complex<float> *Zf = new complex<float>[par.MaxRefIters];
      const count_t M = reference(Zf, par.MaxRefIters, par.Cx, par.Cy, &progress[0], running);
      progress[1] = 1;
      render(out, sta, par, float(Zoom), M, Zf, &progress[2], running);
      delete[] Zf;
    }
    out.saveEXR(par.Stem, par.Channels, threads);

  }
  *ended = true;
}
#endif

int main(int argc, char **argv)
{
  return 0;
#if 0
  using std::isnan;
  using std::isinf;
  using std::log;
  using std::max;
  using std::min;

  const coord_t Width = 1920;
  const coord_t Height = 1080;
  const bool ExponentialMap = false;
  const bool ZoomOutSequence = false;

#if 0
  if (argc != 5)
  {
    std::cerr << "usage: " << argv[0] << "re im zoom iterations out.exr" << std::endl;
    return 1;
  }
  const char *Re = argv[1];
  const char *Im = argv[2];
  floatexp Zoom = atof(argv[3]); // FIXME TODO
  floatexp ZoomPrec = Zoom;
  count_t Iterations = atoi(argv[4]);
  count_t MaxRefIters = Iterations;
  const std::string Stem = argv[5];
#else
#if 1
  // Dinkydau - Flake
  const char *Re = "-1.999966194450370304184346885063505796755312415407248515117619229448015842423426843813761297788689138122870464065609498643538105757447721664856724960928039200953321766548438990841";
  const char *Im = "3.001382436790938324072497303977592498734683119077333527017425728012047497561482358118564729928841407551922418650497818162547852894554815712214577268132087243118270209466408349072e-34";
  floatexp Zoom (e10(2.5620330788506154104770818136626, 157));
  floatexp ZoomPrec (e10(1, 171));
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
  count_t Iterations = 2100100100;
  count_t MaxPtbIters = 2000;
  count_t ReferencePeriod = 43292334;
  count_t MaxRefIters = ReferencePeriod + 1;
  const std::string Stem = argc > 1 ? argv[1] : "out.exr";
#endif
#endif
#endif

  param par;
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

  floatexp pixel_spacing = 4 / ZoomPrec / par.Height;
  mpfr_prec_t prec = 24 - pixel_spacing.exp;
  mpfr_init2(par.Cx, prec);
  mpfr_init2(par.Cy, prec);
  mpfr_set_str(par.Cx, Re, 10, MPFR_RNDN);
  mpfr_set_str(par.Cy, Im, 10, MPFR_RNDN);

  map out(par.Width, par.Height, par.Iterations);

  stats sta;

  progress_t progress[4] = { 0, 0, 0, 0 };
  bool running = true;
  bool ended = false;

  std::thread bg(main_thread, std::ref(out), std::ref(sta), std::cref(par), &progress[0], &running, &ended);

  while (! ended)
  {
    for (count_t ms = 0; ms < 500 && ! ended; ++ms)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    std::cerr
      << "Reference["     << std::setw(3) << int(progress[0] * 100) << "%] "
      << "Frame["         << std::setw(3) << int(progress[1] * 100) << "%] "
      << "Approximation[" << std::setw(3) << int(progress[2] * 100) << "%] "
      << "Pixels["        << std::setw(3) << int(progress[3] * 100) << "%] "
      << "\r";
  }
  if (ended)
  {
    std::cerr << "\n";
  }

  bg.join();

  mpfr_clear(par.Cx);
  mpfr_clear(par.Cy);
  return 0;
#endif
}
