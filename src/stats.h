// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2024 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#if 0

#include "types.h"

struct stat1
{
  double s0, s1, s2, mi, ma;
  inline ~stat1() = default;
  inline stat1(const stat1 &o) = default;
  inline stat1 &operator=(const stat1 &o) = default;
  inline stat1() noexcept
  : s0(0)
  , s1(0)
  , s2(0)
  , mi(1.0/0.0)
  , ma(-1.0/0.0)
  {
  }
  inline stat1(double x) noexcept
  : s0(1)
  , s1(x)
  , s2(x * x)
  , mi(x)
  , ma(x)
  {
  }
  inline stat1& operator+=(const stat1& o) noexcept
  {
    s0 += o.s0;
    s1 += o.s1;
    s2 += o.s2;
    mi = mi < o.mi ? mi : o.mi;
    ma = ma > o.ma ? ma : o.ma;
    return *this;
  }
  inline double mean() noexcept
  {
    return s1 / s0;
  }
  inline double stddev() noexcept
  {
    return std::sqrt(s0 * s2 - s1 * s1) / s0;
  }
};

struct stats
{
  stat1 iters;
  stat1 iters_ptb;
  stat1 iters_bla;
  stat1 steps;
  stat1 steps_ptb;
  stat1 steps_bla;
  stat1 uiters;
  stat1 uiters_ptb;
  stat1 uiters_bla;
  stat1 usteps;
  stat1 usteps_ptb;
  stat1 usteps_bla;
  stat1 iiters;
  stat1 iiters_ptb;
  stat1 iiters_bla;
  stat1 isteps;
  stat1 isteps_ptb;
  stat1 isteps_bla;
  stat1 rebases;
  stat1 rebases_small;
  stat1 rebases_noref;
  stat1 iters_ref;
  inline stats() noexcept
  : iters()
  , iters_ptb()
  , iters_bla()
  , steps()
  , steps_ptb()
  , steps_bla()
  , uiters()
  , uiters_ptb()
  , uiters_bla()
  , usteps()
  , usteps_ptb()
  , usteps_bla()
  , iiters()
  , iiters_ptb()
  , iiters_bla()
  , isteps()
  , isteps_ptb()
  , isteps_bla()
  , rebases()
  , rebases_small()
  , rebases_noref()
  , iters_ref()
  {
  }
  inline stats
  ( double iters
  , double iters_ptb
  , double iters_bla
  , double steps
  , double steps_ptb
  , double steps_bla
  , double rebases
  , double rebases_small
  , double rebases_noref
  , double iters_ref
  , bool escape
  , bool inscape
  ) noexcept
  : iters(escape ? stat1(iters) : stat1())
  , iters_ptb(escape ? stat1(iters_ptb) : stat1())
  , iters_bla(escape ? stat1(iters_bla) : stat1())
  , steps(escape ? stat1(steps) : stat1())
  , steps_ptb(escape ? stat1(steps_ptb) : stat1())
  , steps_bla(escape ? stat1(steps_bla) : stat1())
  , uiters(! (escape || inscape) ? stat1(iters) : stat1())
  , uiters_ptb(! (escape || inscape) ? stat1(iters_ptb) : stat1())
  , uiters_bla(! (escape || inscape) ? stat1(iters_bla) : stat1())
  , usteps(! (escape || inscape) ? stat1(steps) : stat1())
  , usteps_ptb(! (escape || inscape) ? stat1(steps_ptb) : stat1())
  , usteps_bla(! (escape || inscape) ? stat1(steps_bla) : stat1())
  , iiters(inscape ? stat1(iters) : stat1())
  , iiters_ptb(inscape ? stat1(iters_ptb) : stat1())
  , iiters_bla(inscape ? stat1(iters_bla) : stat1())
  , isteps(inscape ? stat1(steps) : stat1())
  , isteps_ptb(inscape ? stat1(steps_ptb) : stat1())
  , isteps_bla(inscape ? stat1(steps_bla) : stat1())
  , rebases(rebases)
  , rebases_small(rebases_small)
  , rebases_noref(rebases_noref)
  , iters_ref(iters_ref)
  {
  }
  inline stats& operator+=(const stats& o) noexcept
  {
    iters += o.iters;
    iters_ptb += o.iters_ptb;
    iters_bla += o.iters_bla;
    steps += o.steps;
    steps_ptb += o.steps_ptb;
    steps_bla += o.steps_bla;
    uiters += o.uiters;
    uiters_ptb += o.uiters_ptb;
    uiters_bla += o.uiters_bla;
    usteps += o.usteps;
    usteps_ptb += o.usteps_ptb;
    usteps_bla += o.usteps_bla;
    iiters += o.iiters;
    iiters_ptb += o.iiters_ptb;
    iiters_bla += o.iiters_bla;
    isteps += o.isteps;
    isteps_ptb += o.isteps_ptb;
    isteps_bla += o.isteps_bla;
    rebases += o.rebases;
    rebases_small += o.rebases_small;
    rebases_noref += o.rebases_noref;
    iters_ref += o.iters_ref;
    return *this;
  }
};

inline void reset(stats &sta) noexcept
{
  sta = stats();
}

#endif

