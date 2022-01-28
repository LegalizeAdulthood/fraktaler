// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <chrono>
#include <cinttypes>
#include <thread>
#include <vector>

#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"
#include "imgui_stdlib.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <SDL2/SDL_opengles2.h>
#else
#include <SDL2/SDL_opengl.h>
#endif
#include <mpreal.h>
#include <omp.h>

#include "colour.h"
#include "display.h"
#include "engine.h"
#include "floatexp.h"
#include "main.h"
#include "map.h"
#include "param.h"
#include "stats.h"

bool newton_ok = true;
count_t newton_period = 0;

void newton_thread(param &out, const param &par, const formula *form, const complex<floatexp> c, const floatexp r, progress_t *progress, bool *running, bool *ended)
{
#if 0
  newton_period = 0;
  if (newton_action >= 0 && *running && newton_ok)
  {
    switch (nt)
    {
      case nt_float:
        newton_period = form->period(&Zf[0], Zf.size(), complex<float>(float(c.x), float(c.y)), par.Iterations, float(r / par.Zoom), mat2<float>(1), &progress[0], running);
        break;
      case nt_double:
        newton_period = form->period(&Zd[0], Zd.size(), complex<double>(double(c.x), double(c.y)), par.Iterations, double(r / par.Zoom), mat2<double>(1), &progress[0], running);
        break;
      case nt_longdouble:
        newton_period = form->period(&Zld[0], Zld.size(), complex<long double>((long double)(c.x), (long double)(c.y)), par.Iterations, (long double)(r / par.Zoom), mat2<long double>(1), &progress[0], running);
        break;
      case nt_floatexp:
        newton_period = form->period(&Zfe[0], Zfe.size(), complex<floatexp>(floatexp(c.x), floatexp(c.y)), par.Iterations, floatexp(r / par.Zoom), mat2<floatexp>(1), &progress[0], running);
        break;
    }
    newton_ok = newton_period > 0;
  }
  if (newton_action >= 1 && *running && newton_ok)
  {
    // FIXME adjust precision factor to correspond to formula/power/etc
    mpfr_prec_t prec = 24 + 2 * std::max(mpfr_get_prec(in.Cx), mpfr_get_prec(in.Cy));
    mpfr_set_prec(out.C.x.mpfr_ptr, prec);
    mpfr_set_prec(out.C.y.mpfr_ptr, prec);
    mpfr_t d;
    mpfr_init2(d, 53);
    mfpr_set_d(d, c.x.val, MPFR_RNDN);
    mpfr_mul_2exp(d, d, c.x.exp, MPFR_RNDN);
    mpfr_add(out.C.x, in.C.x, cx, MPFR_RNDN);
    mfpr_set_d(d, c.y.val, MPFR_RNDN);
    mpfr_mul_2exp(d, d, c.y.exp, MPFR_RNDN);
    mpfr_add(out.C.y, in.C.y, d, MPFR_RNDN);
    mpfr_clear(d);
    newton_ok = form->center(out.C, newton_period, &progress[1], running)
  }
  if (newton_action >= 2 && *running && newton_ok)
  {
    switch (newton_zoom_mode)
    {
      case nr_mode_relative_mini:
        {
          floatexp s (1);
          mat2<double> K (1);
          newton_ok = form->size(s, K, Zp, newton_period, &progress[3], running);
          newton_ok &= 1 / s < Zoom * Zoom * Zoom; // FIXME adjust size sanity check
          if (newton_ok)
          {
            floatexp l0 = log(newton_relative_start);
            floatexp l1 = log(1 / size);
            out.Zoom = exp(l0 + newton_relative_fold * (l1 - l0)) / newton_size_factor;
            mpfr_prec_t prec = 24 - out.Zoom.exp;
            mpfr_prec_round(out.Cx, prec);
            mpfr_prec_round(out.Cy, prec);
          }
        }
        break;
      case nr_mode_absolute_mini:
        {
          floatexp size = form->size(out.Cx, out.Cy, newton_period, &progress[3], running);
          newton_ok = 1 / size < Zoom * Zoom * Zoom; // FIXME adjust size sanity check
          if (newton_ok)
          {
            floatexp l1 = log(1 / size);
            out.Zoom = exp(l1 * newton_absolute_mini_power) / newton_size_factor;
            mpfr_prec_t prec = 24 - out.Zoom.exp;
            mpfr_prec_round(out.Cx, prec);
            mpfr_prec_round(out.Cy, prec);
          }
        }
        break;
      case nr_mode_absolute_domain:
        {
          floatexp size = form->domain_size(out.Cx, out.Cy, newton_period, &progress[3], running);
          newton_ok = 1 / size < Zoom * Zoom * Zoom; // FIXME adjust size sanity check
          if (newton_ok)
          {
            floatexp l1 = log(1 / size);
            out.Zoom = exp(l1 * newton_absolute_domain_power) / newton_size_factor;
            mpfr_prec_t prec = 24 - out.Zoom.exp;
            mpfr_prec_round(out.Cx, prec);
            mpfr_prec_round(out.Cy, prec);
          }
        }
        break;
    }
  }
  if (newton_action >= 3 && *running && ok)
  {
    out.K = form->skew(period, &progress[4], running);
  }
#endif
  *ended = true;
}

bool newton_zoom_enabled = false;
int newton_action = 0;
int newton_zoom_mode = 0;

int newton_relative_preset = 1;
floatexp newton_relative_start = 1;
std::string newton_relative_start_str = "1";
float newton_relative_fold = 0.5;

int newton_absolute_mini_preset = 1;
float newton_absolute_mini_power = 0.5;

int newton_absolute_domain_preset = 1;
float newton_absolute_domain_power = 1.0;

int newton_size_factor_preset = 2;
float newton_size_factor = 4;

bool newton_ball_method = true;

#define nr_action_period 0
#define nr_action_center 1
#define nr_action_size 2
#define nr_action_skew 3

#define nr_mode_relative_mini 0
#define nr_mode_absolute_mini 1
#define nr_mode_absolute_domain 2

void display_newton_window(param &par, bool *open)
{
  ImGui::Begin("Newton Zooming", open);
  ImGui::Checkbox("Activate", &newton_zoom_enabled);
  ImGui::Combo("Action", &newton_action, "Period\0" "Center\0" "Zoom\0" "Skew\0");
  ImGui::Combo("Zoom Mode", &newton_zoom_mode, "Minibrot Relative\0" "Minibrot Absolute\0" "Atom Domain Absolute\0");
  switch (newton_zoom_mode)
  {
    case nr_mode_relative_mini:
      if (InputFloatExp("Relative Start", &newton_relative_start, &newton_relative_start_str))
      {
        std::ostringstream s;
        s << newton_relative_start;
        newton_relative_start_str = s.str();
      }
      ImGui::SameLine();
      if (ImGui::Button("Capture"))
      {
        newton_relative_start = par.Zoom;
        std::ostringstream s;
        s << newton_relative_start;
        newton_relative_start_str = s.str();
      }
      if (ImGui::Combo("Relative Fold", &newton_relative_preset, "Custom\0" "0.5 (2x)\0" "0.75 (4x)\0" "0.875 (8x)\0" "0.9375 (16x)\0" "1.0 (Minibrot)\0"))
      {
        switch (newton_relative_preset)
        {
          case 1: newton_relative_fold = 0.5; break;
          case 2: newton_relative_fold = 0.75; break;
          case 3: newton_relative_fold = 0.875; break;
          case 4: newton_relative_fold = 0.9375; break;
          case 5: newton_relative_fold = 1.0; break;
        }
      }
      if (newton_relative_preset == 0)
      {
        ImGui::SameLine();
        ImGui::InputFloat("##RelativeFoldCustom", &newton_relative_fold);
      }
      break;
    case nr_mode_absolute_mini:
      if (ImGui::Combo("Absolute Power##MiniAbsolutePower", &newton_absolute_mini_preset, "Custom\0" "0.5 (2x)\0" "0.75 (4x)\0" "0.875 (8x)\0" "0.9375 (16x)\0" "1.0 (Minibrot)\0"))
      {
        switch (newton_absolute_mini_preset)
        {
          case 1: newton_absolute_mini_power = 0.5; break;
          case 2: newton_absolute_mini_power = 0.75; break;
          case 3: newton_absolute_mini_power = 0.875; break;
          case 4: newton_absolute_mini_power = 0.9375; break;
          case 5: newton_absolute_mini_power = 1.0; break;
        }
      }
      if (newton_absolute_mini_preset == 0)
      {
        ImGui::SameLine();
        ImGui::InputFloat("##MiniAbsolutePowerCustom", &newton_absolute_mini_power);
      }
      break;
    case nr_mode_absolute_domain:
      if (ImGui::Combo("Absolute Power##DomainAbsolutePower", &newton_absolute_domain_preset, "Custom\0" "1.0 (Domain)\0" "1.125 (Morph)\0"))
      {
        switch (newton_absolute_domain_preset)
        {
          case 1: newton_absolute_domain_power = 1.0; break;
          case 2: newton_absolute_domain_power = 1.125; break;
        }
      }
      if (newton_absolute_domain_preset == 0)
      {
        ImGui::SameLine();
        ImGui::InputFloat("##DomainAbsolutePowerCustom", &newton_absolute_domain_power);
      }
      break;
  }
  if (ImGui::Combo("Size Factor", &newton_size_factor_preset, "Custom\0" "10/1 (zoomed out)\0" "4/1\0" "1/1 (actual size)\0" "1/4\0" "1/10 (zoomed in)\0"))
  {
    switch (newton_size_factor_preset)
    {
      case 1: newton_size_factor = 10; break;
      case 2: newton_size_factor = 4; break;
      case 3: newton_size_factor = 1; break;
      case 4: newton_size_factor = 1./4; break;
      case 5: newton_size_factor = 1./10; break;
    }
  }
  if (newton_size_factor_preset == 0)
  {
    ImGui::SameLine();
    ImGui::InputFloat("##SizeFactorCustom", &newton_size_factor);
  }
  ImGui::Checkbox("Ball Method", &newton_ball_method);
  ImGui::End();
}
