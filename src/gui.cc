// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <chrono>
#include <cinttypes>
#include <thread>
#include <vector>

#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"
#include "imgui_stdlib.h"
#include <SDL.h>
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <SDL_opengles2.h>
#else
#include <SDL_opengl.h>
#endif
#include <mpreal.h>

#ifdef __EMSCRIPTEN__
#include "emscripten.h"
#endif

#ifdef HAVE_OMP
#include <omp.h>
#else
int omp_get_num_procs()
{
  return 1;
}
#endif

#include "colour.h"
#if defined(__EMSCRIPTEN__) || defined(__ANDROID__)
#include "display_web.h"
typedef display_web display_t;
#else
#include "display_gl.h"
typedef display_gl display_t;
#endif
#include "engine.h"
#include "floatexp.h"
#include "formula.h"
#include "main.h"
#include "map.h"
#include "param.h"
#include "stats.h"

#ifdef HAVE_GLDEBUG
#ifdef _WIN32
__attribute__((stdcall))
#endif
static void opengl_debug_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *message, const void *userParam)
{
#ifdef _WIN32
  // don't do anything in Windows
  // in Wine printing happens in multiple threads at once (?) and crashes
  return;
#endif
  (void) userParam;
  (void) length;
  if (source == 33350 && type == 33361 && id == 131185 && severity == 33387)
  {
    // silence extremely verbose message from NVIDIA driver about buffer memory
    return;
  }
  if (source == GL_DEBUG_SOURCE_SHADER_COMPILER && type == GL_DEBUG_TYPE_OTHER && severity == GL_DEBUG_SEVERITY_NOTIFICATION)
  {
    // silence verbose messages from AMDGPU driver about shader compilation
    return;
  }
  const char *source_name = "unknown";
  const char *type_name = "unknown";
  const char *severity_name = "unknown";
  switch (source)
  {
    case GL_DEBUG_SOURCE_API: source_name = "API"; break;
    case GL_DEBUG_SOURCE_WINDOW_SYSTEM: source_name = "Window System"; break;
    case GL_DEBUG_SOURCE_SHADER_COMPILER: source_name = "Shader Compiler"; break;
    case GL_DEBUG_SOURCE_THIRD_PARTY: source_name = "Third Party"; break;
    case GL_DEBUG_SOURCE_APPLICATION: source_name = "Application"; break;
    case GL_DEBUG_SOURCE_OTHER: source_name = "Other"; break;
  }
  switch (type)
  {
    case GL_DEBUG_TYPE_ERROR: type_name = "Error"; break;
    case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: type_name = "Deprecated Behaviour"; break;
    case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR: type_name = "Undefined Behaviour"; break;
    case GL_DEBUG_TYPE_PORTABILITY: type_name = "Portability"; break;
    case GL_DEBUG_TYPE_PERFORMANCE: type_name = "Performance"; break;
    case GL_DEBUG_TYPE_MARKER: type_name = "Marker"; break;
    case GL_DEBUG_TYPE_PUSH_GROUP: type_name = "Push Group"; break;
    case GL_DEBUG_TYPE_POP_GROUP: type_name = "Pop Group"; break;
    case GL_DEBUG_TYPE_OTHER: type_name = "Other"; break;
  }
  switch (severity)
  {
    case GL_DEBUG_SEVERITY_HIGH: severity_name = "High"; break;
    case GL_DEBUG_SEVERITY_MEDIUM: severity_name = "Medium"; break;
    case GL_DEBUG_SEVERITY_LOW: severity_name = "Low"; break;
    case GL_DEBUG_SEVERITY_NOTIFICATION: severity_name = "Notification"; break;
  }
  std::fprintf(stderr, "OpenGL : %s  %s  %s : %d : %s\n", severity_name, type_name, source_name, id, message);
}
#endif

// rendering state machine
progress_t progress[5];
bool quit = false;
bool running = false;
bool restart = false;
bool ended = true;
std::chrono::duration<double> duration = std::chrono::duration<double>::zero();
bool save = false;
bool save_exit = false;

// zoom by mouse drag
bool drag = false;
double drag_start_x = 0;
double drag_start_y = 0;

// last mouse coordinates
int mouse_x = 0;
int mouse_y = 0;

// imgui state
bool show_windows = true;
bool show_status_window = true;
bool show_location_window = true;
bool show_bailout_window = false;
bool show_information_window = true;
bool show_newton_window = false;
bool show_demo_window = false;

int mouse_action = 0;

void handle_event(SDL_Window *window, SDL_Event &e, param &par)
{
#define STOP \
  running = false; \
  while (! ended) \
    std::this_thread::sleep_for(std::chrono::milliseconds(1));

  int win_width = 0;
  int win_height = 0;
  SDL_GetWindowSize(window, &win_width, &win_height);
  SDL_Keymod mods = SDL_GetModState();
  bool shift = mods & KMOD_SHIFT;
  bool ctrl = mods & KMOD_CTRL;
  switch (e.type)
  {
    case SDL_QUIT:
      STOP
      quit = true;
      break;

    case SDL_MOUSEMOTION:
      mouse_x = e.motion.x;
      mouse_y = e.motion.y;
      break;

    case SDL_MOUSEBUTTONDOWN:
      switch (e.button.button)
      {
        case SDL_BUTTON_LEFT:
          drag = true;
          drag_start_x = e.button.x;
          drag_start_y = e.button.y;
          break;
        case SDL_BUTTON_RIGHT:
          if (drag)
          {
            drag = false;
          }
          else
          {
            STOP
            double cx = (e.button.x - win_width / 2.0) / (win_width / 2.0);
            double cy = (e.button.y - win_height / 2.0) / (win_height / 2.0);
            zoom(par, cx, cy, 0.5);
            restart = true;
          }
          break;
        case SDL_BUTTON_MIDDLE:
          {
            STOP
            double cx = (e.button.x - win_width / 2.0) / (win_width / 2.0);
            double cy = (e.button.y - win_height / 2.0) / (win_height / 2.0);
            zoom(par, cx, cy, 1, false);
            restart = true;
          }
        default:
          break;
      }
      break;

    case SDL_MOUSEBUTTONUP:
      switch (e.button.button)
      {
        case SDL_BUTTON_LEFT:
          if (drag)
          {
            double drag_end_x = e.button.x;
            double drag_end_y = e.button.y;
            double dx = (drag_end_x - drag_start_x) / (win_width / 2.0);
            double dy = (drag_end_y - drag_start_y) / (win_height / 2.0);
            double d = std::min(std::max(1 / std::hypot(dx, dy), 1.0/16.0), 16.0);
            double cx = (drag_start_x - win_width / 2.0) / (win_width / 2.0);
            double cy = (drag_start_y - win_height / 2.0) / (win_height / 2.0);
            drag = false;
            switch (mouse_action)
            {
              case 0:
                STOP
                zoom(par, cx, cy, d, false);
                restart = true;
                break;
              case 1:
                STOP
//                START_NEWTON(cx, cy, d)
                break;
            }
          }
          break;
        default:
          break;
      }
      break;

    case SDL_MOUSEWHEEL:
      if (e.wheel.y > 0)
      {
        double cx = (mouse_x - win_width / 2.0) / (win_width / 2.0);
        double cy = (mouse_y - win_height / 2.0) / (win_height / 2.0);
        STOP
        zoom(par, cx, cy, 2);
        restart = true;
      }
      if (e.wheel.y < 0)
      {
        double cx = (mouse_x - win_width / 2.0) / (win_width / 2.0);
        double cy = (mouse_y - win_height / 2.0) / (win_height / 2.0);
        STOP
        zoom(par, cx, cy, 0.5);
        restart = true;
      }
      break;

    case SDL_KEYDOWN:
      switch (e.key.keysym.sym)
      {
        case SDLK_ESCAPE:
          STOP
          break;

        case SDLK_F5:
          STOP
          restart = true;
          break;

        case SDLK_F10:
          show_windows = ! show_windows;
          break;

        case SDLK_KP_PLUS:
          STOP
          if (shift)
          {
            if (par.MaxPtbIters < count_t(1) << 48)
            {
              par.MaxPtbIters <<= 1;
            }
          }
          else
          {
            if (par.Iterations < count_t(1) << 48)
            {
              par.Iterations <<= 1;
              par.MaxRefIters <<= 1;
            }
          }
          restring(par);
          restart = true;
          break;
        case SDLK_KP_MINUS:
          STOP
          if (shift)
          {
            if (par.MaxPtbIters > count_t(1) << 8)
            {
              par.MaxPtbIters >>= 1;
            }
          }
          else
          {
            if (par.Iterations > count_t(1) << 8)
            {
              par.Iterations >>= 1;
              par.MaxRefIters >>= 1;
            }
          }
          restring(par);
          restart = true;
          break;

        case SDLK_KP_0:
          STOP
          zoom(par, 0, 0, 0.5);
          restart = true;
          break;
        case SDLK_KP_1:
          STOP
          zoom(par, -1, 1, 2);
          restart = true;
          break;
        case SDLK_KP_2:
          STOP
          zoom(par, 0, 1, 2);
          restart = true;
          break;
        case SDLK_KP_3:
          STOP
          zoom(par, 1, 1, 2);
          restart = true;
          break;
        case SDLK_KP_4:
          STOP
          zoom(par, -1, 0, 2);
          restart = true;
          break;
        case SDLK_KP_5:
          STOP
          zoom(par, 0, 0, 2);
          restart = true;
          break;
        case SDLK_KP_6:
          STOP
          zoom(par, 1, 0, 2);
          restart = true;
          break;
        case SDLK_KP_7:
          STOP
          zoom(par, -1, -1, 2);
          restart = true;
          break;
        case SDLK_KP_8:
          STOP
          zoom(par, 0, -1, 2);
          restart = true;
          break;
        case SDLK_KP_9:
          STOP
          zoom(par, 1, -1, 2);
          restart = true;
          break;

        case SDLK_HOME:
          STOP
          home(par);
          restart = true;
          break;

        case SDLK_s:
          if (ctrl)
          {
            save = true;
          }
          break;

        default:
          break;
      }
      break;
  }
}

void display_background(SDL_Window *window, display_t &dsp)
{
  int win_width = 0;
  int win_height = 0;
  SDL_GetWindowSize(window, &win_width, &win_height);
  double x0 = 1, y0 = 1, x1 = -1, y1 = -1;
  if (drag)
  {
    double cx = (drag_start_x - win_width / 2.0) / (win_width / 2.0);
    double cy = (drag_start_y - win_height / 2.0) / (win_height / 2.0);
    double r = std::hypot((mouse_x - drag_start_x) / (win_width / 2.0), (mouse_y - drag_start_y) / (win_height / 2.0));
    x0 = cx - r;
    x1 = cx + r;
    y0 = cy - r;
    y1 = cy + r;
  }
  int display_w = 0, display_h = 0;
  SDL_GL_GetDrawableSize(window, &display_w, &display_h);
  dsp.draw(display_w, display_h, x0, y0, x1, y1);
}

void display_window_window()
{
  ImGui::Begin("Windows");
  ImGui::Combo("Mouse Action", &mouse_action, "Navigate\0" "Newton\0");
  ImGui::Checkbox("Status", &show_status_window);
  ImGui::Checkbox("Location", &show_location_window);
  ImGui::Checkbox("Information", &show_information_window);
  ImGui::Checkbox("Newton Zooming", &show_newton_window);
  ImGui::Checkbox("ImGui Demo", &show_demo_window);
  ImGui::Text("Press F10 to toggle all");
  ImGui::End();
}

void display_status_window(bool *open)
{
  char ref[20], apx[20], sub[20], pix[20];
  float r = progress[0];
  float a = progress[2];
  float f = progress[3];
  float p = progress[4];
  std::snprintf(ref, sizeof(ref), "Ref: %3d%%", (int)(r * 100));
  std::snprintf(apx, sizeof(apx), "Apx: %3d%%", (int)(a * 100));
  std::snprintf(sub, sizeof(sub), "Sub: %3d%%", (int)(f * 100));
  std::snprintf(pix, sizeof(pix), "Pix: %3d%%", (int)(p * 100));
  const char *status = "Status: unknown";
  if (! running)
  {
    status = "Cancelled";
  }
  else if (ended)
  {
    status = "Completed";
  }
  else
  {
    status = "Working...";
  }
  ImGui::Begin("Status", open);
  ImGui::Text(status);
  ImGui::ProgressBar(r, ImVec2(-1.f, 0.f), ref);
  ImGui::ProgressBar(a, ImVec2(-1.f, 0.f), apx);
  ImGui::ProgressBar(f, ImVec2(-1.f, 0.f), sub);
  ImGui::ProgressBar(p, ImVec2(-1.f, 0.f), pix);
  count_t ms = std::ceil(1000 * duration.count());
  count_t s = ms / 1000;
  count_t m = s / 60;
  count_t h = m / 60;
  count_t d = h / 24;
  if (d > 0)
  {
    ImGui::Text("T: %dd%02dh%02dm%02ds%03dms", int(d), int(h % 24), int(m % 60), int(s % 60), int(ms % 1000));
  }
  else if (h > 0)
  {
    ImGui::Text("T: %dh%02dm%02ds%03dms", int(h), int(m % 60), int(s % 60), int(ms % 1000));
  }
  else if (m > 0)
  {
    ImGui::Text("T: %dm%02ds%03dms", int(m), int(s % 60), int(ms % 1000));
  }
  else if (s > 0)
  {
    ImGui::Text("T: %ds%03dms", int(s), int(ms % 1000));
  }
  else
  {
    ImGui::Text("T: %dms", int(ms));
  }
  ImGui::End();
}

bool InputFloatExp(const char *label, floatexp *x, std::string *str)
{

  if (ImGui::InputText(label, str, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    mpfr_t z;
    mpfr_init2(z, 53);
    mpfr_set_str(z, str->c_str(), 10, MPFR_RNDN);
    long e = 0;
    double m = mpfr_get_d_2exp(&e, z, MPFR_RNDN);
    mpfr_clear(z);
    *x = floatexp(m, e);
    return true;
  }
  return false;
}

void display_location_window(param &par, bool *open)
{
  ImGui::Begin("Location", open);
  ImGui::Text("Zoom");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Zoom", &par.sZoom, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    STOP
    mpfr_t zoom;
    mpfr_init2(zoom, 53);
    mpfr_set_str(zoom, par.sZoom.c_str(), 10, MPFR_RNDN);
    long e = 0;
    double m = mpfr_get_d_2exp(&e, zoom, MPFR_RNDN);
    mpfr_clear(zoom);
    par.Zoom = floatexp(m, e);
    restring(par);
    restart = true;
  }
  ImGui::PopItemWidth();
  ImGui::Text("Real");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Real", &par.sRe, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    STOP
    mpfr_set_str(par.C.x.mpfr_ptr(), par.sRe.c_str(), 10, MPFR_RNDN);
    restring(par);
    restart = true;
  }
  ImGui::PopItemWidth();
  ImGui::Text("Imag");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Imag", &par.sIm, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    STOP
    mpfr_set_str(par.C.y.mpfr_ptr(), par.sIm.c_str(), 10, MPFR_RNDN);
    restring(par);
    restart = true;
  }
  ImGui::PopItemWidth();
  ImGui::End();
}

void display_bailout_window(param &par, bool *open)
{
  ImGui::Begin("Bailout", open);
  ImGui::Text("Iterations   ");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Iterations", &par.sIterations, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsDecimal))
  {
    try
    {
      count_t tmp = std::stoll(par.sIterations);
      if (tmp > 0)
      {
        STOP
        par.Iterations = tmp;
        restring(par);
        restart = true;
      }
      else
      {
        restring(par);
      }
    }
    catch (std::invalid_argument &e)
    {
      restring(par);
    }
    catch (std::out_of_range &e)
    {
      restring(par);
    }
  }
  ImGui::PopItemWidth();
  ImGui::Text("Max Ref Iters");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##MaxRefIters", &par.sMaxRefIters, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsDecimal))
  {
    try
    {
      count_t tmp = std::stoll(par.sMaxRefIters);
      if (tmp > 0)
      {
        STOP
        par.MaxRefIters = tmp;
        restring(par);
        restart = true;
      }
      else
      {
        restring(par);
      }
    }
    catch (std::invalid_argument &e)
    {
      restring(par);
    }
    catch (std::out_of_range &e)
    {
      restring(par);
    }
  }
  ImGui::PopItemWidth();
  ImGui::Text("Max Ptb Iters");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##MaxPtbIters", &par.sMaxPtbIters, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsDecimal))
  {
    try
    {
      count_t tmp = std::stoll(par.sMaxPtbIters);
      if (tmp > 0)
      {
        STOP
        par.MaxPtbIters = tmp;
        restring(par);
        restart = true;
      }
      else
      {
        restring(par);
      }
    }
    catch (std::invalid_argument &e)
    {
      restring(par);
    }
    catch (std::out_of_range &e)
    {
      restring(par);
    }
  }
  ImGui::PopItemWidth();
  bool LockMaxRefItersToPeriod = par.LockMaxRefItersToPeriod;
  if (ImGui::Checkbox("Lock Max Ref Iters to Period", &LockMaxRefItersToPeriod))
  {
    STOP
    par.LockMaxRefItersToPeriod = LockMaxRefItersToPeriod;
    restart = true;
  }
  ImGui::Text("Escape Radius");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##EscapeRadius", &par.sEscapeRadius, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    try
    {
      double tmp = std::stod(par.sEscapeRadius);
      if (tmp >= 2) // FIXME TODO low bailout radius
      {
        STOP
        par.EscapeRadius = tmp;
        restring(par);
        restart = true;
      }
      else
      {
        restring(par);
      }
    }
    catch (std::invalid_argument &e)
    {
      restring(par);
    }
    catch (std::out_of_range &e)
    {
      restring(par);
    }
  }
  ImGui::PopItemWidth();
  ImGui::End();
}

void display_information_window(stats &sta, bool *open)
{
  ImGui::Begin("Information", open);
  ImGui::Text("Speedup            %.1fx", sta.iterations / (double) (sta.perturb_iterations + sta.bla_steps));
  ImGui::Text("Average Steps      %.1f", (sta.perturb_iterations + sta.bla_steps) / (double) sta.pixels);
  ImGui::Text("Average BLA Steps  %.1f", sta.bla_steps / (double) sta.pixels);
  ImGui::Text("Average Ptb Steps  %.1f", sta.perturb_iterations / (double) sta.pixels);
  ImGui::Text("Iterations Per BLA %.1f", sta.bla_iterations / (double) sta.bla_steps);
  ImGui::Text("Average Iterations %.1f", sta.iterations / (double) sta.pixels);
  ImGui::Text("Average BLA Iters. %.1f", sta.bla_iterations / (double) sta.pixels);
  ImGui::Text("Average Ptb Iters. %.1f", sta.perturb_iterations / (double) sta.pixels);
  ImGui::Text("Average Rebases    %.1f", sta.rebases / (double) sta.pixels);
  ImGui::Text("Minimum Iterations %" PRId64, sta.minimum_iterations);
  ImGui::Text("Maximum Iterations %" PRId64, sta.maximum_iterations);
  ImGui::End();
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

void display_gui(SDL_Window *window, display_t &dsp, param &par, stats &sta)
{
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplSDL2_NewFrame();
  ImGui::NewFrame();

  if (show_windows)
  {
    display_window_window();
    if (show_status_window)
    {
      display_status_window(&show_status_window);
    }
    if (show_location_window)
    {
      display_location_window(par, &show_location_window);
    }
    if (show_bailout_window)
    {
      display_bailout_window(par, &show_bailout_window);
    }
    if (show_information_window)
    {
      display_information_window(sta, &show_information_window);
    }
    if (show_newton_window)
    {
      display_newton_window(par, &show_newton_window);
    }
    if (show_demo_window)
    {
      ImGui::ShowDemoWindow(&show_demo_window);
    }
  }

  ImGui::Render();
  ImGuiIO& io = ImGui::GetIO();
  glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
  glClearColor(0.5, 0.5, 0.5, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  display_background(window, dsp);
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  SDL_GL_SwapWindow(window);
}

bool want_capture(int type)
{
  ImGuiIO& io = ImGui::GetIO();
  return
    (io.WantCaptureMouse && (type == SDL_MOUSEBUTTONDOWN || type == SDL_MOUSEBUTTONUP || type == SDL_MOUSEWHEEL || type == SDL_MOUSEMOTION)) ||
    (io.WantCaptureKeyboard && (type == SDL_KEYDOWN || type == SDL_KEYUP)) ;
}

const coord_t win_width = 1024;
const coord_t win_height = 576;
SDL_Window* window = nullptr;
param par;
stats sta;
formula *form = nullptr;
colour *clr = nullptr;
display_t *dsp = nullptr;
map *out = nullptr;
std::thread *bg = nullptr;
std::chrono::time_point<std::chrono::steady_clock> start_time;
count_t subframe = 0;

enum { st_start, st_reference, st_reference_end, st_subframe_start, st_subframe, st_subframe_end, st_idle, st_quit } state = st_start;

void main1()
{
  switch (state)
  {
    case st_start:
      if (quit)
      {
        state = st_quit;
      }
      else
      {
        progress[0] = 0;
        progress[1] = 0;
        progress[2] = 0;
        progress[3] = 0;
        progress[4] = 0;
        running = true;
        ended = false;
        restart = false;
        start_time = std::chrono::steady_clock::now();
        bg = new std::thread (reference_thread, std::ref(sta), form, std::cref(par), &progress[0], &running, &ended);
        state = st_reference;
      }
      break;
    case st_reference:
      if (quit)
      {
        state = st_quit;
      }
      else if (ended)
      {
        state = st_reference_end;
      }
      else
      {
        auto current_time = std::chrono::steady_clock::now();
        duration = current_time - start_time;
        display_gui(window, *dsp, par, sta);
        SDL_Event e;
        while (SDL_PollEvent(&e))
        {
          ImGui_ImplSDL2_ProcessEvent(&e);
          if (! want_capture(e.type))
          {
            handle_event(window, e, par);
          }
        }
      }
      break;
    case st_reference_end:
      {
        bg->join();
        delete bg;
        bg = nullptr;
        subframe = 0;
        state = st_subframe_start;
      }
      // fall-through
    case st_subframe_start:
      if (running && (par.MaxSubframes <= 0 || subframe < par.MaxSubframes))
      {
        ended = false;
        restart = false;
        progress[3] = par.MaxSubframes <= 0 ? 0 : subframe / progress_t(par.MaxSubframes);
        bg = new std::thread (subframe_thread, std::ref(*out), std::ref(sta), form, std::cref(par), subframe, &progress[4], &running, &ended);
        state = st_subframe;
      }
      else
      {
      }
      break;
    case st_subframe:
      if (! quit && ! ended && ! restart)
      {
        auto current_time = std::chrono::steady_clock::now();
        duration = current_time - start_time;
        display_gui(window, *dsp, par, sta);
        SDL_Event e;
        while (SDL_PollEvent(&e))
        {
          ImGui_ImplSDL2_ProcessEvent(&e);
          if (! want_capture(e.type))
          {
            handle_event(window, e, par);
          }
        }
      }
      else
      {
        state = st_subframe_end;
      }
      break;
    case st_subframe_end:
      {
        bg->join();
        delete bg;
        bg = nullptr;
        if (running) // was not interrupted
        {
          if (subframe == 0)
          {
            dsp->clear();
          }
          dsp->accumulate(*out);
          subframe++;
          if (subframe >= par.MaxSubframes)
          {
            progress[3] = 1;
            state = st_idle;
          }
          else
          {
            state = st_subframe_start;
          }
        }
        else
        {
          state = st_idle;
        }
      }
      // fall-through
    case st_idle:
      if (quit)
      {
        state = st_quit;
      }
      else if (restart)
      {
        state = st_start;
      }
      else
      {
        display_gui(window, *dsp, par, sta);
        if (save)
        {
          dsp->get_rgb(*out);
          out->saveEXR(par.Stem, (1 << Channel_R) | (1 << Channel_G) | (1 << Channel_B), omp_get_num_procs());
          save = false;
          if (save_exit)
          {
            quit = true;
          }
        }
        SDL_Event e;
        while (SDL_PollEvent(&e))
        {
          ImGui_ImplSDL2_ProcessEvent(&e);
          if (! want_capture(e.type))
          {
            handle_event(window, e, par);
          }
        }
      }
      break;
    case st_quit:
      {
      }
      break;
  }
}

int main(int argc, char **argv)
{
  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0)
  {
    std::cerr << argv[0] << ": error: SDL_Init: " << SDL_GetError() << std::endl;
    return 1;
  }

  // decide GL+GLSL versions
#if defined(IMGUI_IMPL_OPENGL_ES2)
  // GL ES 2.0 + GLSL 100
  const char* glsl_version = "#version 100";
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_ES);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
#elif defined(__APPLE__)
  // GL 3.2 Core + GLSL 150
  const char* glsl_version = "#version 150";
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG); // Always required on Mac
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
#else
  // GL 3.0 + GLSL 130
  const char* glsl_version = "#version 130";
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
#endif

  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
  SDL_WindowFlags window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL | /*SDL_WINDOW_RESIZABLE | */SDL_WINDOW_ALLOW_HIGHDPI);
  window = SDL_CreateWindow("Fraktaler 3", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, win_width, win_height, window_flags);
  if (! window)
  {
    std::cerr << argv[0] << ": error: SDL_CreateWindow: " << SDL_GetError() << std::endl;
    SDL_Quit();
    return 1;
  }
  SDL_GLContext gl_context = SDL_GL_CreateContext(window);
  if (! gl_context)
  {
    std::cerr << argv[0] << ": error: SDL_GL_CreateContext: " << SDL_GetError() << std::endl;
    SDL_Quit();
    return 1;
  }
  SDL_GL_MakeCurrent(window, gl_context);

#ifdef HAVE_GLEW
  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK)
  {
    std::cerr << argv[0] << ": error: glewInit" << std::endl;
    SDL_Quit();
    return 1;
  }
  glGetError(); // discard common error from glew
#endif

#ifdef HAVE_GLDEBUG
  if (glDebugMessageCallback)
  {
    glDebugMessageCallback(opengl_debug_callback, nullptr);
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
  }
#endif

  SDL_GL_SetSwapInterval(1);

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
  glClearColor(0.5, 0.5, 0.5, 1);
  glClear(GL_COLOR_BUFFER_BIT);

  GLint maximum_texture_size = 0;
  glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maximum_texture_size);

  // setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO(); (void)io;
  ImGui::StyleColorsDark();
  ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
  ImGui_ImplOpenGL3_Init(glsl_version);

  colours_init();
  formulas_init();

  par.ExponentialMap = false;
  par.ZoomOutSequence = false;
  par.Channels = Channels_default;
  par.Stem = "fraktaler-3.exr";
  par.Width = win_width;
  par.Height = win_height;
  par.EscapeRadius = 625;
  par.MaxSubframes = 16;
  home(par);
  if (argc == 7)
  {
    mpreal radius (2);
    radius.set_prec(53);
    radius = argv[3];
    par.Zoom = floatexp(2 / radius) / 1.5;
    mpfr_prec_t prec = 24 + par.Zoom.exp;
    par.C.x.set_prec(prec);
    par.C.y.set_prec(prec);
    par.C.x = argv[1];
    par.C.y = argv[2];
    par.K = rotation(std::atof(argv[4]));
    par.Iterations = atoll(argv[5]);
    par.MaxRefIters = par.Iterations;
    par.MaxPtbIters = atoll(argv[6]);
    restring(par);
    save = true;
    save_exit = true;
  }

  form = formulas[0]; // FIXME
  clr = colours[0]; // FIXME
  out = new map(par.Width, par.Height, par.Iterations);
  dsp = new display_t(clr);
  dsp->resize(out->width, out->height);
  reset(sta);

#ifdef __EMSCRIPTEN__
  emscripten_set_main_loop(main1, 0, true);
#else
  while (! quit)
  {
    main1();
  }
#endif

  // cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();
  SDL_GL_DeleteContext(gl_context);
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
