// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <chrono>
#include <cinttypes>
#include <map>
#include <thread>
#include <vector>

#include "types.h"

#include <SDL.h>
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <SDL_opengles2.h>
#else
#include <SDL_opengl.h>
#endif
#include <imgui.h>
#include <imgui_impl_sdl.h>
#include <imgui_impl_opengl3.h>
#include <imgui_stdlib.h>
#ifdef HAVE_FS
#include <imfilebrowser.h>
#endif
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/matrix_transform_2d.hpp>
#include <mpreal.h>
#include <toml.hpp>

#ifdef __EMSCRIPTEN__
#include "emscripten.h"
#include "emscripten/html5.h"
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
#include "map.h"
#include "param.h"
#include "stats.h"
#include "version.h"

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

const char *gl_version = "unknown";
int srgb_conversion = 0;

// global state
SDL_Window* window = nullptr;
param par;
stats sta;
colour *clr = nullptr;
display_t *dsp = nullptr;
map *out = nullptr;
std::thread *bg = nullptr;
std::chrono::time_point<std::chrono::steady_clock> start_time;
count_t subframe = 0;

// rendering state machine
std::vector<progress_t> progress;
progress_t newton_progress[4];
bool quit = false;
bool running = false;
bool restart = false;
bool continue_subframe_rendering = false;
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

// touch events
SDL_TouchID finger_device;
std::map<SDL_FingerID, std::pair<vec3, vec3>> fingers;
mat3 finger_transform(1.0f);
mat3 finger_transform_started(1.0f);

bool matrix_ok(const mat3 &m)
{
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j)
  {
    if (std::isnan(m[i][j]))
    {
      return false;
    }
    if (std::isinf(m[i][j]))
    {
      return false;
    }
  }
  if (std::abs(glm::determinant(m)) < 1.0e-6f)
  {
    return false;
  }
  return true;
}

void update_finger_transform()
{
  switch (fingers.size())
  {
    case 0: // identity
      {
        mat3 T = mat3(1.0f);
        if (matrix_ok(T))
        {
          finger_transform = T * finger_transform;
        }
      }
      break;
    case 1: // translate
      {
        const std::pair<vec3, vec3> &finger = (*fingers.begin()).second;
        const vec3 start = finger.first;
        const vec3 end = finger.second;
        mat3 T = mat3(1.0f);
        T = glm::translate(T, - vec2(start[0], start[1]) / start.z);
        T = glm::translate(T, vec2(end[0], end[1]) / end.z);
        if (matrix_ok(T))
        {
          finger_transform = T * finger_transform;
        }
      }
      break;
    case 2: // translate, rotate, scale
      {
        const std::pair<vec3, vec3> &finger1 = (*fingers.begin()).second;
        const std::pair<vec3, vec3> &finger2 = (*++fingers.begin()).second;
        const vec3 start1 = finger1.first;
        const vec3 end1 = finger1.second;
        const vec3 start2 = finger2.first;
        const vec3 end2 = finger2.second;
        const vec2 p1 = vec2(start1[0], start1[1]) / start1.z;
        const vec2 q1 = vec2(end1[0], end1[1]) / end1.z;
        const vec2 p2 = vec2(start2[0], start2[1]) / start2.z;
        const vec2 q2 = vec2(end2[0], end2[1]) / end2.z;
        const vec2 dp = p1 - p2;
        const vec2 dq = q1 - q2;
        const vec2 mp = (p1 + p2) * 0.5f;
        const vec2 mq = (q1 + q2) * 0.5f;
        const float ap = std::atan2(dp.y, dp.x);
        const float aq = std::atan2(dq.y, dq.x);
        const float sp = std::hypot(dp.y, dp.x);
        const float sq = std::hypot(dq.y, dq.x);
        const float t = aq - ap;
        const float s = sq / sp;
        const float co = s * std::cos(t);
        const float si = s * std::sin(t);
        const float px = mp[0];
        const float py = mp[1];
        const float qx = mq[0];
        const float qy = mq[1];
        mat3 N(1.0f, 0.0f, -px,   0.0f, 1.0f, -py,   0.0f, 0.0f, 1.0f);
        mat3 M(co,   -si,  0.0f,  si,   co,   0.0f,  0.0f, 0.0f, 1.0f);
        mat3 L(1.0f, 0.0f, qx,    0.0f, 1.0f, qy,    0.0f, 0.0f, 1.0f);
        mat3 T = transpose(L) * transpose(M) * transpose(N);
        if (matrix_ok(T))
        {
          finger_transform = T * finger_transform;
        }
      }
      break;
    default: // overdetermined system, just use first 3 fingers....
    case 3: // translate, rotate, scale, skew
      {
        const std::pair<vec3, vec3> &finger1 = (*fingers.begin()).second;
        const std::pair<vec3, vec3> &finger2 = (*++fingers.begin()).second;
        const std::pair<vec3, vec3> &finger3 = (*++++fingers.begin()).second;
        const mat3 start(finger1.first, finger2.first, finger3.first);
        const mat3 end(finger1.second, finger2.second, finger3.second);
        mat3 T = end * inverse(start);
        if (matrix_ok(T))
        {
          finger_transform = T * finger_transform;
        }
      }
      break;
  }
  for (auto & [k, finger] : fingers)
  {
    finger.first = finger.second;
  }
}

// imgui state

bool show_windows = true;

struct window
{
  bool show = false;
  int x = 64, y = 64, w = 256, h = 256;
};

struct windows
{
  struct window
    io = { false, -1, -1, 167, 54 },
    formula = { false, -1, -1, 225, 54 },
    colour = { false, -1, -1, 225, 54 },
    status = { false, -1, -1, 128, 160 },
    location = { false, -1, -1, -1, 104 },
    reference = { false, -1, -1, -1, 104 },
    algorithm = { false, -1, -1, 235, 254 },
    bailout = { false, -1, -1, 240, 123 },
    transform = { false, -1, -1, 274, 146 },
    information = { false, -1, -1, 442, 218 },
    quality = { false, -1, -1, 218, 77 },
    newton = { false, -1, -1, 384, 384 },
#ifdef HAVE_IMGUI_DEMO
    demo = { false, -1, -1, -1, -1 },
#endif
    about = { false, -1, -1, 576, -1 };
};

std::istream &operator>>(std::istream &ifs, windows &w)
{
  auto t = toml::parse(ifs);
#define LOAD(a) \
  w.a.show = toml::find_or(t, "window", #a, "show", w.a.show); \
  w.a.x = toml::find_or(t, "window", #a, "x", w.a.x); \
  w.a.y = toml::find_or(t, "window", #a, "y", w.a.y); \
  w.a.w = toml::find_or(t, "window", #a, "w", w.a.w); \
  w.a.h = toml::find_or(t, "window", #a, "h", w.a.h);
  LOAD(io)
  LOAD(formula)
  LOAD(colour)
  LOAD(status)
  LOAD(location)
  LOAD(reference)
  LOAD(algorithm)
  LOAD(bailout)
  LOAD(transform)
  LOAD(information)
  LOAD(quality)
  LOAD(newton)
#ifdef HAVE_IMGUI_DEMO
  LOAD(demo)
#endif
  LOAD(about)
#undef LOAD
  return ifs;
}

std::ostream &operator<<(std::ostream &ofs, const windows &p)
{
  const windows q;
#define SAVE2(a,b) if (p.a.b != q.a.b) { ofs << "window" << "." << #a << "." << #b << " = " << std::setw(70) << toml::value(p.a.b) << "\n"; }
#define SAVE(a) SAVE2(a, show) SAVE2(a, x) SAVE2(a, y) SAVE2(a, w) SAVE2(a, h)
  SAVE(io)
  SAVE(formula)
  SAVE(colour)
  SAVE(status)
  SAVE(location)
  SAVE(reference)
  SAVE(algorithm)
  SAVE(bailout)
  SAVE(transform)
  SAVE(information)
  SAVE(quality)
  SAVE(newton)
#ifdef HAVE_IMGUI_DEMO
  SAVE(demo)
#endif
  SAVE(about)
#undef SAVE
#undef SAVE2
  return ofs;
}

windows window_state;

enum user_event_codes
{
  code_persist = 1
};

Uint32 one_minute = 60 * 1000;
Uint32 persistence_timer_callback(Uint32 interval, void *p)
{
  (void) p;
  SDL_Event event;
  SDL_UserEvent userevent;
  userevent.type = SDL_USEREVENT;
  userevent.code = code_persist;
  userevent.data1 = 0;
  userevent.data2 = 0;
  event.type = SDL_USEREVENT;
  event.user = userevent;
  SDL_PushEvent(&event);
  return interval;
}

int mouse_action = 0;

const SDL_TouchID multitouch_device = 147;
SDL_FingerID multitouch_finger = 0;
std::map<SDL_FingerID, std::pair<coord_t, coord_t>> multitouch_fingers;

void multitouch_move_finger(SDL_FingerID finger, coord_t x, coord_t y)
{
  multitouch_fingers[finger] = std::pair<coord_t, coord_t>(x, y);
}

SDL_FingerID multitouch_add_finger(coord_t x, coord_t y)
{
  float md2 = 1.0f/0.0f;
  SDL_FingerID finger = 0;
  for (const auto & [id, p] : multitouch_fingers)
  {
    float dx = x - p.first;
    float dy = y - p.second;
    float d2 = dx * dx + dy * dy;
    if (d2 < md2)
    {
      md2 = d2;
      finger = id;
    }
  }
  if (md2 > 16 * 16) // FIXME hardcoded sensitivity
  {
    // none nearby, add one
    finger = 1;
    for (const auto & [id, v] : multitouch_fingers)
    {
      if (id == finger)
      {
        finger++;
      }
      else
      {
        break;
      }
    }
    multitouch_fingers[finger] = std::pair<coord_t, coord_t>(x, y);
  }
  else
  {
  }
  return finger;
}

SDL_FingerID multitouch_remove_finger(coord_t &x, coord_t &y)
{
  float md2 = 1.0f/0.0f;
  coord_t mx = x;
  coord_t my = y;
  SDL_FingerID finger = 0;
  for (const auto & [id, p] : multitouch_fingers)
  {
    float dx = x - p.first;
    float dy = y - p.second;
    float d2 = dx * dx + dy * dy;
    if (d2 < md2)
    {
      md2 = d2;
      finger = id;
      mx = p.first;
      my = p.second;
    }
  }
  if (finger)
  {
    x = mx;
    y = my;
    multitouch_fingers.erase(finger);
  }
  return finger;
}

int win_pixel_width = 0;
int win_pixel_height = 0;
void resize(coord_t super, coord_t sub)
{
  auto width = (win_pixel_width * super) / sub;
  auto height = (win_pixel_height * super) / sub;
  delete out;
  out = new map(width, height, par.p.bailout.iterations);
  dsp->resize(out->width, out->height);
}

void persist_state()
{
  try
  {
    par.save_toml(pref_path + "persistence.toml");
    syncfs();
  }
  catch (const std::exception &e)
  {
    SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "saving persistence %s", e.what());
  }
  try
  {
    std::ofstream ofs;
    ofs.exceptions(std::ofstream::badbit | std::ofstream::failbit);
    ofs.open(pref_path + "gui-settings.toml", std::ios_base::binary);
    ofs << window_state;
  }
  catch (const std::exception &e)
  {
    SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "saving GUI settings %s", e.what());
  }
}

complex<floatexp> newton_c(0);
floatexp newton_r = 1;
bool start_newton = false;
bool newton_running = false;
bool newton_ended = false;
param newton_par;
bool newton_ok = false;

void handle_event(SDL_Window *window, SDL_Event &e, param &par)
{
#ifdef __EMSCRIPTEN__
#define STOP \
  running = false; \
  while (! ended) \
    emscripten_sleep(1);
#else
#define STOP \
  running = false; \
  while (! ended) \
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
#endif

  int win_width = 0;
  int win_height = 0;
  SDL_GetWindowSize(window, &win_width, &win_height);
  SDL_Keymod mods = SDL_GetModState();
  bool shift = mods & KMOD_SHIFT;
  bool ctrl = mods & KMOD_CTRL;
  bool multitouch_emulation = ctrl && shift;
  switch (e.type)
  {
    case SDL_QUIT:
      STOP
      quit = true;
      persist_state();
      break;

    case SDL_MOUSEMOTION:
      if (e.motion.which != SDL_TOUCH_MOUSEID)
      {
        if (multitouch_emulation)
        {
          if (e.motion.state & SDL_BUTTON_LMASK)
          {
            multitouch_move_finger(multitouch_finger, e.motion.x, e.motion.y);
            SDL_Event t;
            t.type = SDL_FINGERMOTION;
            t.tfinger.type = SDL_FINGERMOTION;
            t.tfinger.timestamp = e.motion.timestamp;
            t.tfinger.touchId = multitouch_device;
            t.tfinger.fingerId = multitouch_finger;
            t.tfinger.x = e.motion.x / float(win_width);
            t.tfinger.y = e.motion.y / float(win_height);
            t.tfinger.dx = e.motion.xrel / float(win_width);
            t.tfinger.dy = e.motion.yrel / float(win_height);
            t.tfinger.pressure = 0.5f; // FIXME
#if 0
            t.tfinger.windowID = e.motion.windowID;
#endif
            SDL_PushEvent(&t);
          }
        }
        else
        {
          mouse_x = e.motion.x;
          mouse_y = e.motion.y;
        }
      }
      break;

    case SDL_MOUSEBUTTONDOWN:
      if (e.button.which != SDL_TOUCH_MOUSEID)
      {
        if (multitouch_emulation)
        {
          if (e.button.button == SDL_BUTTON_LEFT)
          {
            multitouch_finger = multitouch_add_finger(e.button.x, e.button.y);
            SDL_Event t;
            t.type = SDL_FINGERDOWN;
            t.tfinger.type = SDL_FINGERDOWN;
            t.tfinger.timestamp = e.button.timestamp;
            t.tfinger.touchId = multitouch_device;
            t.tfinger.fingerId = multitouch_finger;
            t.tfinger.x = e.button.x / float(win_width);
            t.tfinger.y = e.button.y / float(win_height);
            t.tfinger.dx = 0.0f;
            t.tfinger.dy = 0.0f;
            t.tfinger.pressure = 0.5f; // FIXME
#if 0
            t.tfinger.windowID = e.button.windowID;
#endif
            SDL_PushEvent(&t);
          }
          else if (e.button.button == SDL_BUTTON_RIGHT)
          {
            coord_t x = e.button.x;
            coord_t y = e.button.y;
            multitouch_finger = multitouch_remove_finger(x, y);
            if (multitouch_finger)
            {
              SDL_Event t;
              t.type = SDL_FINGERUP;
              t.tfinger.type = SDL_FINGERUP;
              t.tfinger.timestamp = e.button.timestamp;
              t.tfinger.touchId = multitouch_device;
              t.tfinger.fingerId = multitouch_finger;
              t.tfinger.x = x / float(win_width);
              t.tfinger.y = y / float(win_height);
              t.tfinger.dx = 0.0f;
              t.tfinger.dy = 0.0f;
              t.tfinger.pressure = 0.5f; // FIXME
#if 0
              t.tfinger.windowID = e.button.windowID;
#endif
              SDL_PushEvent(&t);
            }
          }
        }
        else
        {
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
                mat3 T = mat3(1.0f);
                T = glm::translate(T, vec2(float(e.button.x), float(win_height - e.button.y)));
                T = glm::scale(T, vec2(float(0.5f), float(0.5f)));
                T = glm::translate(T, -vec2(float(e.button.x), float(win_height - e.button.y)));
                finger_transform_started = T * finger_transform_started;
                restart = true;
              }
              break;
            case SDL_BUTTON_MIDDLE:
              {
                STOP
                double cx = (e.button.x - win_width / 2.0) / (win_width / 2.0);
                double cy = (e.button.y - win_height / 2.0) / (win_height / 2.0);
                zoom(par, cx, cy, 1, false);
                mat3 T = mat3(1.0f);
                T = glm::translate(T, -vec2(float(e.button.x - win_width / 2.0), float(win_height - e.button.y - win_height / 2.0)));
                finger_transform_started = T * finger_transform_started;
                restart = true;
              }
            default:
              break;
          }
        }
      }
      break;

    case SDL_MOUSEBUTTONUP:
      if (e.button.which != SDL_TOUCH_MOUSEID)
      {
        if (multitouch_emulation)
        {
        }
        else
        {
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
                    {
                      STOP
                      zoom(par, cx, cy, d, false);
                      mat3 T = mat3(1.0f);
                      T = glm::translate(T, vec2(float(win_width / 2.0), float(win_height / 2.0)));
                      T = glm::scale(T, vec2(float(d), float(d)));
                      T = glm::translate(T, -vec2(float(win_width / 2.0), float(win_height / 2.0)));
                      T = glm::translate(T, -vec2(float(drag_start_x - win_width / 2.0), float(win_height -drag_start_y - win_height / 2.0)));
                      finger_transform_started = T * finger_transform_started;
                      restart = true;
                    }
                    break;
                  case 1:
                    {
                      STOP
                      newton_c = get_delta_c(par, cx, cy);
                      newton_r = d / par.zoom;
                      start_newton = true;
                    }
                    break;
                }
              }
              break;
            default:
              break;
          }
        }
      }
      break;

    case SDL_MOUSEWHEEL:
      if (e.wheel.which != SDL_TOUCH_MOUSEID)
      {
        if (e.wheel.y > 0)
        {
          double cx = (mouse_x - win_width / 2.0) / (win_width / 2.0);
          double cy = (mouse_y - win_height / 2.0) / (win_height / 2.0);
          STOP
          zoom(par, cx, cy, 2);
          mat3 T = mat3(1.0f);
          T = glm::translate(T, vec2(float(mouse_x), float(win_height - mouse_y)));
          T = glm::scale(T, vec2(float(2), float(2)));
          T = glm::translate(T, -vec2(float(mouse_x), float(win_height - mouse_y)));
          finger_transform_started = T * finger_transform_started;
          restart = true;
        }
        if (e.wheel.y < 0)
        {
          double cx = (mouse_x - win_width / 2.0) / (win_width / 2.0);
          double cy = (mouse_y - win_height / 2.0) / (win_height / 2.0);
          STOP
          zoom(par, cx, cy, 0.5);
          mat3 T = mat3(1.0f);
          T = glm::translate(T, vec2(float(mouse_x), float(win_height - mouse_y)));
          T = glm::scale(T, vec2(float(0.5f), float(0.5f)));
          T = glm::translate(T, -vec2(float(mouse_x), float(win_height - mouse_y)));
          finger_transform_started = T * finger_transform_started;
          restart = true;
        }
      }
      break;

    case SDL_FINGERDOWN:
      switch (mouse_action)
      {
        case 0:
          {
            if (fingers.size() == 0)
            {
              finger_device = e.tfinger.touchId;
            }
            if (finger_device == e.tfinger.touchId)
            {
              vec3 f = vec3(e.tfinger.x * win_width, (1 - e.tfinger.y) * win_height, 1.0f);
              fingers[e.tfinger.fingerId] = std::pair<vec3, vec3>(f, f);
              update_finger_transform();
            }
          }
          break;
        case 1:
          {
            STOP
            newton_c = get_delta_c(par, (e.tfinger.x - 0.5) * 2, (0.5 - e.tfinger.y) * 2);
            newton_r = 0.1 / par.zoom;
            start_newton = true;
          }
          break;
      }
      break;
    case SDL_FINGERUP:
      if (finger_device == e.tfinger.touchId)
      {
        vec3 f = vec3(e.tfinger.x * win_width, (1 - e.tfinger.y) * win_height, 1.0f);
        fingers[e.tfinger.fingerId].second = f;
        update_finger_transform();
        fingers.erase(e.tfinger.fingerId);
        if (fingers.size() == 0)
        {
          STOP
          mat3 S = mat3(1.0f);
          // [0..w] x [0..h]
          S = glm::scale(S, vec2(float(win_width), float(win_height)));
          S = glm::scale(S, vec2(0.5f, 0.5f));
          S = glm::translate(S, vec2(1.0f));
          // [-1..1] x [-1..1]
          S = glm::inverse(S) * finger_transform * S;
          zoom(par, glm::inverse(S), finger_transform);
          finger_transform_started = finger_transform * finger_transform_started;
          finger_transform = mat3(1.0f);
          restart = true;
        }
      }
      break;

    case SDL_FINGERMOTION:
      if (finger_device == e.tfinger.touchId)
      {
        vec3 f = vec3(e.tfinger.x * win_width, (1 - e.tfinger.y) * win_height, 1.0f);
        fingers[e.tfinger.fingerId].second = f;
        update_finger_transform();
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
            if (par.p.bailout.maximum_perturb_iterations < count_t(1) << 48)
            {
              par.p.bailout.maximum_perturb_iterations <<= 1;
            }
          }
          else
          {
            if (par.p.bailout.iterations < count_t(1) << 48)
            {
              par.p.bailout.iterations <<= 1;
              par.p.bailout.maximum_reference_iterations <<= 1;
            }
          }
          restring_vals(par);
          restart = true;
          break;
        case SDLK_KP_MINUS:
          STOP
          if (shift)
          {
            if (par.p.bailout.maximum_perturb_iterations > count_t(1) << 8)
            {
              par.p.bailout.maximum_perturb_iterations >>= 1;
            }
          }
          else
          {
            if (par.p.bailout.iterations > count_t(1) << 8)
            {
              par.p.bailout.iterations >>= 1;
              par.p.bailout.maximum_reference_iterations >>= 1;
            }
          }
          restring_vals(par);
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

        case SDLK_q:
          if (ctrl)
          {
            quit = true;
          }
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

    case SDL_USEREVENT:
      switch (e.user.code)
      {
        case code_persist:
          persist_state();
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
  int display_w = 0, display_h = 0;
  SDL_GL_GetDrawableSize(window, &display_w, &display_h);
  dsp.draw(display_w, display_h, finger_transform * finger_transform_started, srgb_conversion);
  if (drag)
  {
    double cx = (drag_start_x - win_width / 2.0) / (win_width / 2.0);
    double cy = (drag_start_y - win_height / 2.0) / (win_height / 2.0);
    double r = std::hypot((mouse_x - drag_start_x) / (win_width / 2.0), (mouse_y - drag_start_y) / (win_height / 2.0));
    double x0 = cx - r;
    double x1 = cx + r;
    double y0 = cy - r;
    double y1 = cy + r;
    dsp.draw_rectangle(display_w, display_h, x0, y0, x1, y1, srgb_conversion);
  }
  if (fingers.size() > 0)
  {
    std::vector<glm::vec4> circles;
    float rx = 0.01;
    float ry = rx * win_width / win_height;
    for (const auto &f : fingers)
    {
      float cx = (f.second.second[0] - win_width / 2.0) / (win_width / 2.0);
      float cy = (f.second.second[1] - win_height / 2.0) / (win_height / 2.0);
      circles.push_back(glm::vec4(cx, cy, rx, ry));
    }
    dsp.draw_circles(display_w, display_h, circles, srgb_conversion);
  }
}

void display_window_window()
{
  ImGui::SetNextWindowPos(ImVec2(16, 16), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(177, 324), ImGuiCond_FirstUseEver);
  ImGui::Begin("Fraktaler 3");
//  ImGui::Combo("##MouseAction", &mouse_action, "Navigate\0");// "Newton\0");
  ImGui::Checkbox("Input/Ouput", &window_state.io.show);
  ImGui::Checkbox("Formula", &window_state.formula.show);
  ImGui::Checkbox("Colour", &window_state.colour.show);
  ImGui::Checkbox("Status", &window_state.status.show);
  ImGui::Checkbox("Location", &window_state.location.show);
  ImGui::Checkbox("Reference", &window_state.reference.show);
  ImGui::Checkbox("Bailout", &window_state.bailout.show);
  ImGui::Checkbox("Transform", &window_state.transform.show);
  ImGui::Checkbox("Algorithm", &window_state.algorithm.show);
  ImGui::Checkbox("Information", &window_state.information.show);
  ImGui::Checkbox("Quality", &window_state.quality.show);
  ImGui::Checkbox("Newton Zooming", &window_state.newton.show);
  ImGui::Checkbox("About", &window_state.about.show);
#ifdef HAVE_IMGUI_DEMO
  ImGui::Checkbox("ImGui Demo", &window_state.demo.show);
#endif
  ImGui::Text("Press F10 to toggle all");
  ImGui::End();
}

#ifdef HAVE_FS
ImGui::FileBrowser *load_dialog = nullptr;
ImGui::FileBrowser *save_dialog = nullptr;
#endif

bool reset_unlocked = false;

void display_set_window_dims(const struct window &w)
{
  const auto &io = ImGui::GetIO();
  int width = w.w > 0 ? w.w : io.DisplaySize.x - 16 * 2;
  int height = w.h > 0 ? w.h : io.DisplaySize.y - 16 * 2;
  int x = w.x > 0 ? w.x : (io.DisplaySize.x - width) / 2;
  int y = w.y > 0 ? w.y : (io.DisplaySize.y - height) / 2;
  ImGui::SetNextWindowPos(ImVec2(x, y), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(width, height), ImGuiCond_FirstUseEver);
}

void display_get_window_dims(struct window &w)
{
  ImVec2 p = ImGui::GetWindowPos();
  ImVec2 s = ImGui::GetWindowSize();
  w.x = p.x;
  w.y = p.y;
  w.w = s.x;
  w.h = s.y;
}

void display_io_window(bool *open)
{
  display_set_window_dims(window_state.io);
  ImGui::Begin("Input/Ouput", open);
  display_get_window_dims(window_state.io);
  ImGui::Checkbox("##ResetUnlocked", &reset_unlocked);
  ImGui::SameLine();
  if (ImGui::Button("Home") && reset_unlocked)
  {
    STOP
    reset_unlocked = false;
    home(par);
    restart = true;
  }
#ifdef HAVE_FS
  ImGui::SameLine();
  if (ImGui::Button("Load"))
  {
    load_dialog->Open();
  }
  ImGui::SameLine();
  if (ImGui::Button("Save"))
  {
    save_dialog->Open();
  }
#endif
  ImGui::End();
#ifdef HAVE_FS
  load_dialog->Display();
  save_dialog->Display();
  if (load_dialog->HasSelected())
  {
    std::string filename = load_dialog->GetSelected().string();
    try
    {
      STOP
      par.load_toml(filename);
      clr = colours[par.p.colour_id];
      dsp->set_colour(clr);
      par.p.image.subsampling = std::min(std::max(par.p.image.subsampling, 1), 32); // FIXME
      resize(1, par.p.image.subsampling);
      restart = true;
    }
    catch (std::exception &e)
    {
      SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "loading \"%s\": %s", filename.c_str(), e.what());
    }
    load_dialog->ClearSelected();
  }
  if (save_dialog->HasSelected())
  {
    std::string filename = save_dialog->GetSelected().string();
    try
    {
      par.save_toml(filename);
      syncfs();
    }
    catch (const std::exception &e)
    {
      SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "saving \"%s\": %s", filename.c_str(), e.what());
    }
    save_dialog->ClearSelected();
  }
#endif
}

void display_status_window(bool *open)
{
  const count_t count = par.p.formula.per.size();
  const char *status = "Status: unknown";
  if (! running)
  {
    status = "Cancelled";
  }
  else if (ended && (par.p.image.subframes > 0 && subframe >= par.p.image.subframes))
  {
    status = "Completed";
  }
  else
  {
    status = "Working...";
  }
  display_set_window_dims(window_state.status);
  ImGui::Begin("Status", open);
  display_get_window_dims(window_state.status);
  ImGui::Text("%s", status);
  float r = 0;
  for (int i = 0; i < count; ++i)
  {
    r += progress[i];
  }
  r /= count;
  char ref[20];
  std::snprintf(ref, sizeof(ref), "Ref: %3d%%", (int)(r * 100));
  ImGui::ProgressBar(r, ImVec2(-1.f, 0.f), ref);
  float a = 0;
  for (int i = 0; i < count; ++i)
  {
    a += progress[count + i];
  }
  a /= count;
  char apx[20];
  std::snprintf(apx, sizeof(apx), "Apx: %3d%%", (int)(a * 100));
  ImGui::ProgressBar(a, ImVec2(-1.f, 0.f), apx);
  char sub[20];
  float f = par.p.image.subframes == 0 ? 0 : glm::clamp(progress[2 * count], 0.0f, 1.0f);
  std::snprintf(sub, sizeof(sub), "Sub: %d/%d", (int) subframe, par.p.image.subframes);
  ImGui::ProgressBar(f, ImVec2(-1.f, 0.f), sub);
  char pix[20];
  float p = progress[2 * count + 1];
  std::snprintf(pix, sizeof(pix), "Pix: %3d%%", (int)(p * 100));
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

void display_formula_window(param &par, bool *open)
{
  display_set_window_dims(window_state.formula);
  ImGui::Begin("Formula", open);
  display_get_window_dims(window_state.formula);
  auto f = par.p.formula.per;
  const count_t count = f.size();
  bool changed = false;
  for (count_t i = 0; i < count; ++i)
  {
    ImGui::PushID(i);
    changed |= ImGui::Checkbox("|X|", &f[i].abs_x); ImGui::SameLine();
    changed |= ImGui::Checkbox("|Y|", &f[i].abs_y); ImGui::SameLine();
    changed |= ImGui::Checkbox("-X", &f[i].neg_x); ImGui::SameLine();
    changed |= ImGui::Checkbox("-Y", &f[i].neg_y); ImGui::SameLine();
    ImGui::PushItemWidth(100);
    changed |= ImGui::InputInt("P", &f[i].power, 1, 5); ImGui::SameLine();
    if (f[i].power < 2)
    {
      f[i].power = 2;
    }
    ImGui::PopItemWidth();
    if (ImGui::Button("+"))
    {
      f.insert(f.begin() + i, f[i]);
      changed |= true;
    }
    if (1 < count)
    {
      ImGui::SameLine();
      if (ImGui::Button("-"))
      {
        f.erase(f.begin() + i);
        changed |= true;
      }
    }
    ImGui::PopID();
  }
  if (changed)
  {
    STOP
    par.p.formula.per = f;
    restart = true;
  }
  ImGui::End();
}

bool colour_get_name(void *data, int n, const char **out_str)
{
  (void) data;
  *out_str = colours[n]->name();
  return true;
}

void display_colour_window(param &par, bool *open)
{
  display_set_window_dims(window_state.colour);
  ImGui::Begin("Colour", open);
  display_get_window_dims(window_state.colour);
  int colour_id = par.p.colour_id;
  if (ImGui::Combo("Colour", &colour_id, colour_get_name, &colours, colours.size()))
  {
    if (colour_id != par.p.colour_id)
    {
      STOP
      par.p.colour_id = colour_id;
      clr = colours[par.p.colour_id];
      dsp->set_colour(clr);
      restart = true;
    }
  }
  ImGui::End();
}

void display_location_window(param &par, bool *open)
{
  display_set_window_dims(window_state.location);
  ImGui::Begin("Location", open);
  display_get_window_dims(window_state.location);
  ImGui::Text("Zoom");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Zoom", &par.p.location.zoom, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    STOP
    mpfr_t zoom;
    mpfr_init2(zoom, 53);
    mpfr_set_str(zoom, par.p.location.zoom.c_str(), 10, MPFR_RNDN);
    long e = 0;
    double m = mpfr_get_d_2exp(&e, zoom, MPFR_RNDN);
    mpfr_clear(zoom);
    par.zoom = floatexp(m, e);
    restring_vals(par);
    restart = true;
  }
  ImGui::PopItemWidth();
  ImGui::Text("Real");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Real", &par.p.location.real, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    STOP
    mpfr_set_str(par.center.x.mpfr_ptr(), par.p.location.real.c_str(), 10, MPFR_RNDN);
    restring_locs(par);
    restart = true;
  }
  ImGui::PopItemWidth();
  ImGui::Text("Imag");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Imag", &par.p.location.imag, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    STOP
    mpfr_set_str(par.center.y.mpfr_ptr(), par.p.location.imag.c_str(), 10, MPFR_RNDN);
    restring_locs(par);
    restart = true;
  }
  ImGui::PopItemWidth();
  ImGui::End();
}

void display_reference_window(param &par, bool *open)
{
  display_set_window_dims(window_state.reference);
  ImGui::Begin("Reference", open);
  display_get_window_dims(window_state.reference);
  ImGui::Text("Period");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Period", &par.s_period, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    STOP
    unstring_vals(par);
    restart = true;
  }
  ImGui::PopItemWidth();
  ImGui::Text("Real");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Real", &par.p.reference.real, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    STOP
    unstring_locs(par);
    restart = true;
  }
  ImGui::PopItemWidth();
  ImGui::Text("Imag");
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Imag", &par.p.reference.imag, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    STOP
    unstring_locs(par);
    restart = true;
  }
  ImGui::PopItemWidth();
  ImGui::End();
}

void display_bailout_window(param &par, bool *open)
{
  display_set_window_dims(window_state.bailout);
  ImGui::Begin("Bailout", open);
  display_get_window_dims(window_state.bailout);
  ImGui::Text("Iterations   ");
  ImGui::SameLine();
  if (ImGui::Button("-##IterationsDown"))
  {
    STOP
    par.p.bailout.iterations >>= 1;
    par.p.bailout.iterations = std::max(par.p.bailout.iterations, count_t(1) << 6);
    restring_vals(par);
    restart = true;
  }
  ImGui::SameLine();
  if (ImGui::Button("+##IterationsUp"))
  {
    STOP
    par.p.bailout.iterations <<= 1;
    par.p.bailout.iterations = std::min(par.p.bailout.iterations, count_t(1) << 60);
    restring_vals(par);
    restart = true;
  }
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##Iterations", &par.s_iterations, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsDecimal))
  {
    try
    {
      count_t tmp = std::stoll(par.s_iterations);
      if (tmp > 0)
      {
        STOP
        par.p.bailout.iterations = tmp;
        restring_vals(par);
        restart = true;
      }
      else
      {
        restring_vals(par);
      }
    }
    catch (std::invalid_argument &e)
    {
      restring_vals(par);
    }
    catch (std::out_of_range &e)
    {
      restring_vals(par);
    }
  }
  ImGui::PopItemWidth();
  ImGui::Text("Max Ref Iters");
  ImGui::SameLine();
  if (ImGui::Button("-##MaxRefItersDown"))
  {
    STOP
    par.p.bailout.maximum_reference_iterations >>= 1;
    par.p.bailout.maximum_reference_iterations = std::max(par.p.bailout.maximum_reference_iterations, count_t(1) << 6);
    restring_vals(par);
    restart = true;
  }
  ImGui::SameLine();
  if (ImGui::Button("+##MaxRefItersUp"))
  {
    STOP
    par.p.bailout.maximum_reference_iterations <<= 1;
    par.p.bailout.maximum_reference_iterations = std::min(par.p.bailout.maximum_reference_iterations, count_t(1) << 60);
    restring_vals(par);
    restart = true;
  }
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##MaxRefIters", &par.s_maximum_reference_iterations, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsDecimal))
  {
    try
    {
      count_t tmp = std::stoll(par.s_maximum_reference_iterations);
      if (tmp > 0)
      {
        STOP
        par.p.bailout.maximum_reference_iterations = tmp;
        restring_vals(par);
        restart = true;
      }
      else
      {
        restring_vals(par);
      }
    }
    catch (std::invalid_argument &e)
    {
      restring_vals(par);
    }
    catch (std::out_of_range &e)
    {
      restring_vals(par);
    }
  }
  ImGui::PopItemWidth();
  ImGui::Text("Max Ptb Iters");
  ImGui::SameLine();
  if (ImGui::Button("-##MaxPtbItersDown"))
  {
    STOP
    par.p.bailout.maximum_perturb_iterations >>= 1;
    par.p.bailout.maximum_perturb_iterations = std::max(par.p.bailout.maximum_perturb_iterations, count_t(1) << 6);
    restring_vals(par);
    restart = true;
  }
  ImGui::SameLine();
  if (ImGui::Button("+##MaxPtbItersUp"))
  {
    STOP
    par.p.bailout.maximum_perturb_iterations <<= 1;
    par.p.bailout.maximum_perturb_iterations = std::min(par.p.bailout.maximum_perturb_iterations, count_t(1) << 60);
    restring_vals(par);
    restart = true;
  }
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##MaxPtbIters", &par.s_maximum_perturb_iterations, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsDecimal))
  {
    try
    {
      count_t tmp = std::stoll(par.s_maximum_perturb_iterations);
      if (tmp > 0)
      {
        STOP
        par.p.bailout.maximum_perturb_iterations = tmp;
        restring_vals(par);
        restart = true;
      }
      else
      {
        restring_vals(par);
      }
    }
    catch (std::invalid_argument &e)
    {
      restring_vals(par);
    }
    catch (std::out_of_range &e)
    {
      restring_vals(par);
    }
  }
  ImGui::PopItemWidth();
  ImGui::Text("Escape Radius");
  ImGui::SameLine();
  if (ImGui::Button("-##EscapeRadiusDown"))
  {
    STOP
    par.p.bailout.escape_radius /= 2;
    par.p.bailout.escape_radius = std::max(par.p.bailout.escape_radius, 2.0);
    restring_vals(par);
    restart = true;
  }
  ImGui::SameLine();
  if (ImGui::Button("+##EscapeRadiusUp"))
  {
    STOP
    par.p.bailout.escape_radius *= 2;
    par.p.bailout.escape_radius = std::min(par.p.bailout.escape_radius, 65536.0);
    restring_vals(par);
    restart = true;
  }
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##EscapeRadius", &par.s_escape_radius, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    try
    {
      double tmp = std::stod(par.s_escape_radius);
      if (tmp >= 2) // FIXME TODO low bailout radius
      {
        STOP
        par.p.bailout.escape_radius = tmp;
        restring_vals(par);
        restart = true;
      }
      else
      {
        restring_vals(par);
      }
    }
    catch (std::invalid_argument &e)
    {
      restring_vals(par);
    }
    catch (std::out_of_range &e)
    {
      restring_vals(par);
    }
  }
  ImGui::PopItemWidth();
  ImGui::Text("Inscape Radius");
  ImGui::SameLine();
  if (ImGui::Button("-##InscapeRadiusDown"))
  {
    STOP
    par.p.bailout.inscape_radius /= 2;
    par.p.bailout.inscape_radius = std::max(par.p.bailout.inscape_radius, 1.0 / (1 << 20));
    restring_vals(par);
    restart = true;
  }
  ImGui::SameLine();
  if (ImGui::Button("+##InscapeRadiusUp"))
  {
    STOP
    par.p.bailout.inscape_radius *= 2;
    par.p.bailout.inscape_radius = std::min(par.p.bailout.inscape_radius, 1.0 / (1 << 0));
    restring_vals(par);
    restart = true;
  }
  ImGui::SameLine();
  ImGui::PushItemWidth(-FLT_MIN);
  if (ImGui::InputText("##InscapeRadius", &par.s_inscape_radius, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific))
  {
    try
    {
      double tmp = std::stod(par.s_inscape_radius);
      if (tmp <= 1)
      {
        STOP
        par.p.bailout.inscape_radius = tmp;
        restring_vals(par);
        restart = true;
      }
      else
      {
        restring_vals(par);
      }
    }
    catch (std::invalid_argument &e)
    {
      restring_vals(par);
    }
    catch (std::out_of_range &e)
    {
      restring_vals(par);
    }
  }
  ImGui::PopItemWidth();
  ImGui::End();
}

void display_transform_window(param &par, bool *open)
{
  display_set_window_dims(window_state.transform);
  ImGui::Begin("Transform", open);
  display_get_window_dims(window_state.transform);
  bool reflect = par.p.transform.reflect;
  if (ImGui::Checkbox("Reflect", &reflect))
  {
    STOP
    par.p.transform.reflect = reflect;
    restring_vals(par);
    restart = true;
  }
  {
    float rotate = par.p.transform.rotate;
    bool changed = false;
    if (ImGui::Button("-##RotateDown"))
    {
      rotate -= 5;
      if (rotate < -360) rotate += 720;
      changed = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("0##RotateZero"))
    {
      rotate = 0;
      changed = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("+##RotateUp"))
    {
      rotate += 5;
      if (rotate > 360) rotate -= 720;
      changed = true;
    }
    ImGui::SameLine();
    ImGui::PushItemWidth(200);
    if (ImGui::SliderFloat("Rotate", &rotate, -360.f, 360.f, "%.2f") || changed)
    {
      STOP
      par.p.transform.rotate = rotate;
      restring_vals(par);
      restart = true;
    }
    ImGui::PopItemWidth();
  }
  if (ImGui::Button("Auto Stretch (DE)"))
  {
    if (subframe > 0)
    {
      STOP
      mat2<double> T = out->unskewDE();
      double d = determinant(T);
      if (! (std::isnan(d) || std::isinf(d) || d == 0))
      {
        par.transform = par.transform * T;
        unstring_vals(par);
        restart = true;
      }
    }
  }
  {
    float stretch_amount = par.p.transform.stretch_amount;
    bool changed = false;
    if (ImGui::Button("-##AmountDown"))
    {
      stretch_amount -= 5;
      if (stretch_amount < -1000) stretch_amount = -1000;
      changed = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("0##AmountZero"))
    {
      stretch_amount = 0;
      changed = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("+##AmountUp"))
    {
      stretch_amount += 5;
      if (stretch_amount > 1000) stretch_amount = 1000;
      changed = true;
    }
    ImGui::SameLine();
    ImGui::PushItemWidth(200);
    if (ImGui::SliderFloat("Stretch Amount", &stretch_amount, -1000.f, 1000.f, "%.2f") || changed) // FIXME
    {
      STOP
      par.p.transform.stretch_amount = stretch_amount;
      restring_vals(par);
      restart = true;
    }
    ImGui::PopItemWidth();
  }
  {
    float stretch_angle = par.p.transform.stretch_angle;
    bool changed = false;
    if (ImGui::Button("-##AngleDown"))
    {
      stretch_angle -= 5;
      if (stretch_angle < -360) stretch_angle += 720;
      changed = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("0##AngleZero"))
    {
      stretch_angle = 0;
      changed = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("+##AngleUp"))
    {
      stretch_angle += 5;
      if (stretch_angle > 360) stretch_angle -= 720;
      changed = true;
    }
    ImGui::SameLine();
    ImGui::PushItemWidth(200);
    if (ImGui::SliderFloat("Stretch Angle", &stretch_angle, -360.f, 360.f, "%.2f") || changed)
    {
      STOP
      par.p.transform.stretch_angle = stretch_angle;
      restring_vals(par);
      restart = true;
    }
    ImGui::PopItemWidth();
  }
  bool exponential_map = par.p.transform.exponential_map;
  if (ImGui::Checkbox("Exponential Map", &exponential_map))
  {
    STOP
    par.p.transform.exponential_map = exponential_map;
    restring_vals(par);
    restart = true;
  }
  ImGui::End();
}

void display_algorithm_window(param &par, bool *open)
{
  display_set_window_dims(window_state.algorithm);
  ImGui::Begin("Algorithm", open);
  display_get_window_dims(window_state.algorithm);
  bool lock_maximum_reference_iterations_to_period = par.p.algorithm.lock_maximum_reference_iterations_to_period;
  if (ImGui::Checkbox("Lock Max Ref Iters to Period", &lock_maximum_reference_iterations_to_period))
  {
    STOP
    par.p.algorithm.lock_maximum_reference_iterations_to_period = lock_maximum_reference_iterations_to_period;
    restart = true;
  }
#if 0
  bool reuse_reference = par.p.algorithm.reuse_reference;
  if (ImGui::Checkbox("Reuse Reference", &reuse_reference))
  {
    STOP
    par.p.algorithm.reuse_reference = reuse_reference;
    restart = true;
  }
  bool reuse_bilinear_approximation = par.p.algorithm.reuse_bilinear_approximation;
  if (ImGui::Checkbox("Reuse Bilinear Approximation", &reuse_bilinear_approximation))
  {
    STOP
    par.p.algorithm.reuse_bilinear_approximation = reuse_bilinear_approximation;
    restart = true;
  }
#endif
  ImGui::Text("Used Number Type: %s", nt_string[nt_current]);
  ImGui::Text("Number Type Selection");
  // number types drag-and-drop between two columns (left active, right unused)
  std::string names[6 * 2] = { "", "", "", "", "", "", "", "", "", "", "", "" };
  int i = 0;
  for (auto s : par.p.algorithm.number_types)
  {
    names[i] = s;
    i += 2;
  }
  i = 1;
  for (auto cp : nt_string)
  {
    std::string s(cp);
    if (s != "none" && std::find(par.p.algorithm.number_types.begin(), par.p.algorithm.number_types.end(), s) == par.p.algorithm.number_types.end())
    {
      names[i] = s;
      i += 2;
    }
  }
  for (int n = 0; n <
#ifdef HAVE_FLOAT128
  6
#else
  5
#endif
    * 2; n++)
  {
    ImGui::PushID(n);
    if ((n % 2) != 0)
    {
      ImGui::SameLine();
    }
    ImGui::Button(names[n].c_str(), ImVec2(100, 20));
    if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None))
    {
      ImGui::SetDragDropPayload("ALGORITHM_NUMBER_TYPE_CELL", &n, sizeof(int));
      ImGui::Text("Swap");
      ImGui::EndDragDropSource();
    }
    if (ImGui::BeginDragDropTarget())
    {
      if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("ALGORITHM_NUMBER_TYPE_CELL"))
      {
        IM_ASSERT(payload->DataSize == sizeof(int));
        int payload_n = *(const int*)payload->Data;
        std::swap(names[n], names[payload_n]);
      }
      ImGui::EndDragDropTarget();
    }
    ImGui::PopID();
  }
  std::vector<std::string> number_types;
  for (i = 0; i < 6 * 2; i += 2)
  {
    if (names[i] != "")
    {
      number_types.push_back(names[i]);
    }
  }
  if (number_types != par.p.algorithm.number_types)
  {
    STOP
    par.p.algorithm.number_types = number_types;
    restart = true;
  }
  ImGui::End();
}

void display_information_window(stats &sta, bool *open)
{
  display_set_window_dims(window_state.information);
  ImGui::Begin("Information", open);
  display_get_window_dims(window_state.information);
  double count = sta.iiters.s0 + sta.uiters.s0 + sta.iters.s0;
  ImGui::Text("Speedup       %.1fx", sta.iters.mean() / sta.steps.mean());
  ImGui::Text("Escaped       %.1f%%", 100.0 * sta.iters.s0 / count);
  ImGui::Text("Unscaped      %.1f%%", 100.0 * sta.uiters.s0 / count);
  ImGui::Text("Inscaped      %.1f%%", 100.0 * sta.iiters.s0 / count);
  ImGui::Text("Ex. Steps     %.1f (min %.1f, max %.1f, stddev %.1f)", sta.steps.mean(), sta.steps.mi, sta.steps.ma, sta.steps.stddev());
  ImGui::Text("Ex. Steps BLA %.1f (min %.1f, max %.1f, stddev %.1f)", sta.steps_bla.mean(), sta.steps_bla.mi, sta.steps_bla.ma, sta.steps_bla.stddev());
  ImGui::Text("Ex. Steps Ptb %.1f (min %.1f, max %.1f, stddev %.1f)", sta.steps_ptb.mean(), sta.steps_ptb.mi, sta.steps_ptb.ma, sta.steps_ptb.stddev());
  ImGui::Text("Ex. Iters     %.1f (min %.1f, max %.1f, stddev %.1f)", sta.iters.mean(), sta.iters.mi, sta.iters.ma, sta.iters.stddev());
  ImGui::Text("Ex. Iters BLA %.1f (min %.1f, max %.1f, stddev %.1f)", sta.iters_bla.mean(), sta.iters_bla.mi, sta.iters_bla.ma, sta.iters_bla.stddev());
  ImGui::Text("Ex. Iters Ptb %.1f (min %.1f, max %.1f, stddev %.1f)", sta.iters_ptb.mean(), sta.iters_ptb.mi, sta.iters_ptb.ma, sta.iters_ptb.stddev());
  ImGui::Text("Un. Steps     %.1f (min %.1f, max %.1f, stddev %.1f)", sta.usteps.mean(), sta.usteps.mi, sta.usteps.ma, sta.usteps.stddev());
  ImGui::Text("Un. Steps BLA %.1f (min %.1f, max %.1f, stddev %.1f)", sta.usteps_bla.mean(), sta.usteps_bla.mi, sta.usteps_bla.ma, sta.usteps_bla.stddev());
  ImGui::Text("Un. Steps Ptb %.1f (min %.1f, max %.1f, stddev %.1f)", sta.usteps_ptb.mean(), sta.usteps_ptb.mi, sta.usteps_ptb.ma, sta.usteps_ptb.stddev());
  ImGui::Text("Un. Iters     %.1f (min %.1f, max %.1f, stddev %.1f)", sta.uiters.mean(), sta.uiters.mi, sta.uiters.ma, sta.uiters.stddev());
  ImGui::Text("Un. Iters BLA %.1f (min %.1f, max %.1f, stddev %.1f)", sta.uiters_bla.mean(), sta.uiters_bla.mi, sta.uiters_bla.ma, sta.uiters_bla.stddev());
  ImGui::Text("Un. Iters Ptb %.1f (min %.1f, max %.1f, stddev %.1f)", sta.uiters_ptb.mean(), sta.uiters_ptb.mi, sta.uiters_ptb.ma, sta.uiters_ptb.stddev());
  ImGui::Text("In. Steps     %.1f (min %.1f, max %.1f, stddev %.1f)", sta.isteps.mean(), sta.isteps.mi, sta.isteps.ma, sta.isteps.stddev());
  ImGui::Text("In. Steps BLA %.1f (min %.1f, max %.1f, stddev %.1f)", sta.isteps_bla.mean(), sta.isteps_bla.mi, sta.isteps_bla.ma, sta.isteps_bla.stddev());
  ImGui::Text("In. Steps Ptb %.1f (min %.1f, max %.1f, stddev %.1f)", sta.isteps_ptb.mean(), sta.isteps_ptb.mi, sta.isteps_ptb.ma, sta.isteps_ptb.stddev());
  ImGui::Text("In. Iters     %.1f (min %.1f, max %.1f, stddev %.1f)", sta.iiters.mean(), sta.iiters.mi, sta.iiters.ma, sta.iiters.stddev());
  ImGui::Text("In. Iters BLA %.1f (min %.1f, max %.1f, stddev %.1f)", sta.iiters_bla.mean(), sta.iiters_bla.mi, sta.iiters_bla.ma, sta.iiters_bla.stddev());
  ImGui::Text("In. Iters Ptb %.1f (min %.1f, max %.1f, stddev %.1f)", sta.iiters_ptb.mean(), sta.iiters_ptb.mi, sta.iiters_ptb.ma, sta.iiters_ptb.stddev());
  ImGui::Text("Iters Ref Max %.1f (min %.1f, max %.1f, stddev %.1f)", sta.iters_ref.mean(), sta.iters_ref.mi, sta.iters_ref.ma, sta.iters_ref.stddev());
  ImGui::Text("Rebases       %.1f (min %.1f, max %.1f, stddev %.1f)", sta.rebases.mean(), sta.rebases.mi, sta.rebases.ma, sta.rebases.stddev());
  ImGui::Text("Rebases Small %.1f (min %.1f, max %.1f, stddev %.1f)", sta.rebases_small.mean(), sta.rebases_small.mi, sta.rebases_small.ma, sta.rebases_small.stddev());
  ImGui::Text("Rebases NoRef %.1f (min %.1f, max %.1f, stddev %.1f)", sta.rebases_noref.mean(), sta.rebases_noref.mi, sta.rebases_noref.ma, sta.rebases_noref.stddev());
  
  ImGui::End();
}

void display_quality_window(bool *open)
{
  display_set_window_dims(window_state.quality);
  ImGui::Begin("Quality", open);
  display_get_window_dims(window_state.quality);
  int subsampling = par.p.image.subsampling;
  if (ImGui::InputInt("Sub", &subsampling))
  {
    STOP
    par.p.image.subsampling = std::min(std::max(subsampling, 1), 32); // FIXME
    resize(1, par.p.image.subsampling);
    restart = true;
  }
  int subframes = par.p.image.subframes;
  if (ImGui::InputInt("Frames", &subframes))
  {
    par.p.image.subframes = std::min(std::max(subframes, 0), 65536); // FIXME
    continue_subframe_rendering = true;
  }
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

void display_newton_window(param &par, bool *open)
{
  display_set_window_dims(window_state.newton);
  ImGui::Begin("Newton Zooming", open);
  display_get_window_dims(window_state.newton);
  ImGui::Checkbox("Activate", &newton_zoom_enabled);
  mouse_action = newton_zoom_enabled ? 1 : 0;
  ImGui::Combo("Action", &par.p.newton.action, "Period\0" "Center\0" "Zoom\0" "Transform\0");
  int absolute = par.p.newton.absolute;
  ImGui::Combo("Zoom Mode", &absolute, "Relative\0" "Absolute\0");
  par.p.newton.absolute = absolute;
  if (! par.p.newton.absolute)
  {
      if (InputFloatExp("Relative Start", &newton_relative_start, &newton_relative_start_str))
      {
        std::ostringstream s;
        s << newton_relative_start;
        newton_relative_start_str = s.str();
      }
      ImGui::SameLine();
      if (ImGui::Button("Capture"))
      {
        newton_relative_start = par.zoom;
        std::ostringstream s;
        s << newton_relative_start;
        newton_relative_start_str = s.str();
      }
  }
  int power_preset = 0;
  float power_presets[6] = { par.p.newton.power, 0.5f, 0.75f, 0.875f, 0.9375f, 1.0f };
  for (int i = 1; i < 6; ++i)
  {
    if (par.p.newton.power == power_presets[i])
    {
      power_preset = i;
      break;
    }
  }
  if (ImGui::Combo("Power", &power_preset, "Custom\0" "0.5 (2x)\0" "0.75 (4x)\0" "0.875 (8x)\0" "0.9375 (16x)\0" "1.0 (Minibrot)\0"))
  {
    par.p.newton.power = power_presets[power_preset];
  }
  if (power_preset == 0)
  {
    ImGui::SameLine();
    ImGui::InputFloat("##PowerCustom", &par.p.newton.power);
  }
#if 0
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
#endif
  int factor_preset = 0;
  float factor_presets[6] = { par.p.newton.factor, 10.0f, 4.0f, 1.0f, 0.25f, 0.1f };
  for (int i = 1; i < 6; ++i)
  {
    if (par.p.newton.factor == factor_presets[i])
    {
      factor_preset = i;
      break;
    }
  }
  if (ImGui::Combo("Factor", &factor_preset, "Custom\0" "10/1 (zoomed out)\0" "4/1\0" "1/1 (actual size)\0" "1/4\0" "1/10 (zoomed in)\0"))
  {
    par.p.newton.factor = factor_presets[factor_preset];
    switch (newton_size_factor_preset)
    {
      case 1: newton_size_factor = 10; break;
      case 2: newton_size_factor = 4; break;
      case 3: newton_size_factor = 1; break;
      case 4: newton_size_factor = 1./4; break;
      case 5: newton_size_factor = 1./10; break;
    }
  }
  if (factor_preset == 0)
  {
    ImGui::SameLine();
    ImGui::InputFloat("##FactorCustom", &par.p.newton.factor);
  }
  ImGui::End();
}

bool newton_really_stop = false;
void display_newton_modal(bool &open)
{
  if (open)
  {
    ImGui::OpenPopup("Newton Zooming Progress");
  }
  ImVec2 center(ImGui::GetIO().DisplaySize.x * 0.5f, ImGui::GetIO().DisplaySize.y * 0.5f);
  ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
  if (ImGui::BeginPopupModal("Newton Zooming Progress", NULL, ImGuiWindowFlags_AlwaysAutoResize))
  {
    if (par.p.newton.action >= newton_action_period)
    {
      char str[30];
      std::snprintf(str, 30, "Period: %d%%", int(100 * newton_progress[0]));
      ImGui::ProgressBar(float(newton_progress[0]), ImVec2(200,0), &str[0]);
    }
    if (par.p.newton.action >= newton_action_center)
    {
      char str[30];
      std::snprintf(str, 30, "Steps: %d%%", int(100 * newton_progress[1]));
      ImGui::ProgressBar(float(newton_progress[1]), ImVec2(200,0), &str[0]);
      std::snprintf(str, 30, "Center: %d%%", int(100 * newton_progress[2]));
      ImGui::ProgressBar(float(newton_progress[2]), ImVec2(200,0), &str[0]);
    }
    if (par.p.newton.action >= newton_action_zoom)
    {
      char str[30];
      std::snprintf(str, 30, "Size: %d%%", int(100 * newton_progress[3]));
      ImGui::ProgressBar(float(newton_progress[3]), ImVec2(200,0), &str[0]);
    }
    ImGui::Checkbox("##ReallyStop", &newton_really_stop);
    ImGui::SameLine();
    if (ImGui::Button("Stop") && newton_really_stop)
    {
      newton_running = false;
    }
    if (newton_ended)
    {
      open = false;
      ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
  }
}

std::string about_text = "";

void display_about_window(bool *open)
{
  if (about_text == "")
  {
    about_text = version(gl_version) + "\n\n\n\n" + license();
  }
  display_set_window_dims(window_state.about);
  ImGui::Begin("About", open);
  display_get_window_dims(window_state.about);
  ImGui::TextUnformatted(about_text.c_str());
  ImGui::End();
}

void display_gui(SDL_Window *window, display_t &dsp, param &par, stats &sta, bool newton_modal = false)
{
  int win_screen_width = 0;
  int win_screen_height = 0;
  SDL_GetWindowSize(window, &win_screen_width, &win_screen_height);
  SDL_GL_GetDrawableSize(window, &win_pixel_width, &win_pixel_height);
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplSDL2_NewFrame();
  ImGui::NewFrame();
  if (show_windows)
  {
    display_window_window();
    if (window_state.io.show)
    {
      display_io_window(&window_state.io.show);
    }
    if (window_state.status.show)
    {
      display_status_window(&window_state.status.show);
    }
    if (window_state.formula.show)
    {
      display_formula_window(par, &window_state.formula.show);
    }
    if (window_state.colour.show)
    {
      display_colour_window(par, &window_state.colour.show);
    }
    if (window_state.location.show)
    {
      display_location_window(par, &window_state.location.show);
    }
    if (window_state.reference.show)
    {
      display_reference_window(par, &window_state.reference.show);
    }
    if (window_state.bailout.show)
    {
      display_bailout_window(par, &window_state.bailout.show);
    }
    if (window_state.transform.show)
    {
      display_transform_window(par, &window_state.transform.show);
    }
    if (window_state.algorithm.show)
    {
      display_algorithm_window(par, &window_state.algorithm.show);
    }
    if (window_state.information.show)
    {
      display_information_window(sta, &window_state.information.show);
    }
    if (window_state.quality.show)
    {
      display_quality_window(&window_state.quality.show);
    }
    if (window_state.newton.show)
    {
      display_newton_window(par, &window_state.newton.show);
    }
    if (window_state.about.show)
    {
      display_about_window(&window_state.about.show);
    }
#ifdef HAVE_IMGUI_DEMO
    if (window_state.demo.show)
    {
      ImGui::ShowDemoWindow(&window_state.demo.show);
    }
#endif
  }
  display_newton_modal(newton_modal);

  ImGui::Render();
  glViewport(0, 0, win_pixel_width, win_pixel_height);
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
    (io.WantCaptureMouse && (
      type == SDL_MOUSEBUTTONDOWN ||
      type == SDL_MOUSEBUTTONUP ||
      type == SDL_MOUSEWHEEL ||
      type == SDL_MOUSEMOTION ||
      type == SDL_FINGERDOWN ||
      type == SDL_FINGERUP ||
      type == SDL_FINGERMOTION)) ||
    (io.WantCaptureKeyboard && (
      type == SDL_KEYDOWN ||
      type == SDL_KEYUP)) ;
}

enum { st_start, st_reference, st_reference_end, st_subframe_start, st_subframe, st_subframe_end, st_idle, st_quit, st_newton_start, st_newton, st_newton_end } state = st_start;

int gui_busy = 2;
bool just_did_newton = false;
void main1()
{
  const count_t count = par.p.formula.per.size();
  switch (state)
  {
    case st_start:
      if (quit)
      {
        state = st_quit;
      }
      else
      {
        progress.resize(2 * count + 2);
        for (int i = 0; i < 2 * count + 2; ++i)
        {
          progress[i] = 0;
        }
        running = true;
        ended = false;
        restart = false;
        continue_subframe_rendering = false;
        start_time = std::chrono::steady_clock::now();
        bg = new std::thread (reference_thread, std::ref(sta), std::ref(par), just_did_newton, &progress[0], &running, &ended);
        state = st_reference;
        just_did_newton = false;
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
          gui_busy = 2;
        }
      }
      break;
    case st_reference_end:
      {
        bg->join();
        delete bg;
        bg = nullptr;
        restring_locs(par); // reference may have been updated to image center
        subframe = 0;
        state = st_subframe_start;
      }
      // fall-through
    case st_subframe_start:
      if (running && (par.p.image.subframes <= 0 || subframe < par.p.image.subframes))
      {
        ended = false;
        restart = false;
        progress[2 * count] = par.p.image.subframes <= 0 ? 0 : subframe / progress_t(par.p.image.subframes);
        bg = new std::thread (subframe_thread, 0, std::ref(*out), std::ref(sta), std::cref(par), subframe, &progress[2 * count + 1], &running, &ended);
        state = st_subframe;
      }
      else
      {
        state = st_idle;
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
          gui_busy = 2;
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
            finger_transform_started = mat3(1.0f);
            dsp->clear();
          }
          dsp->accumulate(*out);
          subframe++;
          if (par.p.image.subframes > 0 && subframe >= par.p.image.subframes)
          {
            progress[2 * count] = 1;
            state = st_idle;
          }
          else
          {
            continue_subframe_rendering = true;
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
      else if (start_newton)
      {
        state = st_newton_start;
      }
      else if (restart)
      {
        state = st_start;
      }
      else if (continue_subframe_rendering)
      {
        continue_subframe_rendering = false;
        state = st_subframe_start;
      }
      else
      {
        display_gui(window, *dsp, par, sta);
        if (save)
        {
          dsp->get_rgb(*out);
          out->saveEXR(par.p.render.filename + ".exr", Channels_RGB, std::thread::hardware_concurrency());
          save = false;
          if (save_exit)
          {
            quit = true;
          }
        }
        bool got_event = false;
        do
        {
          got_event = false;
          gui_busy--;
          if (gui_busy < 0)
          {
            gui_busy = 0;
          }
          SDL_Event e;
          if (gui_busy ? (got_event = SDL_PollEvent(&e)) : SDL_WaitEvent(&e))
          {
            ImGui_ImplSDL2_ProcessEvent(&e);
            if (! want_capture(e.type))
            {
              handle_event(window, e, par);
            }
            gui_busy = 2;
          }
        }
        while (got_event);
      }
      break;
    case st_quit:
      {
      }
      break;
    case st_newton_start:
      {
        newton_running = true;
        newton_ended = false;
        start_newton = false;
        newton_progress[0] = 0;
        newton_progress[1] = 0;
        newton_progress[2] = 0;
        newton_progress[3] = 0;
        newton_c.x = newton_c.x - floatexp(par.reference.x - par.center.x);
        newton_c.y = newton_c.y - floatexp(par.reference.y - par.center.y);
        bg = new std::thread (newton_thread, std::ref(newton_par), std::ref(newton_ok), std::cref(par), std::cref(newton_c), std::cref(newton_r), &newton_progress[0], &newton_running, &newton_ended);
        state = st_newton;
      }
      break;
    case st_newton:
      if (! quit && ! newton_ended)
      {
        display_gui(window, *dsp, par, sta, true);
        SDL_Event e;
        while (SDL_PollEvent(&e))
        {
          ImGui_ImplSDL2_ProcessEvent(&e);
          if (! want_capture(e.type))
          {
            handle_event(window, e, par);
          }
        }
        gui_busy = 2;
      }
      else
      {
        state = st_newton_end;
      }
      break;
    case st_newton_end:
      {
        bg->join();
        delete bg;
        bg = nullptr;
        if (newton_ok)
        {
          par = newton_par;
          just_did_newton = true;
          state = st_start;
        }
        else
        {
          state = st_idle;
        }
      }
      break;
  }

}

#ifdef __EMSCRIPTEN__
int main0(int argc, char **argv);
extern "C" int EMSCRIPTEN_KEEPALIVE main00(void)
{
  int argc = 1;
  char *argv[] = { "fraktaler-3", nullptr };
  return main0(argc, argv);
}
#endif

int main0(int argc, char **argv)
{
  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0)
  {
    const std::string message = "SDL_Init: " + std::string(SDL_GetError());
    if (0 != SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Fraktaler 3", message.c_str(), nullptr))
    {
      std::cerr << argv[0] << ": " << message << std::endl;
    }
    return 1;
  }

#ifdef __ANDROID__
  SDL_DisplayMode mode;
  if (SDL_GetDesktopDisplayMode(0, &mode) != 0)
  {
    const std::string message = "SDL_GetDesktopDisplayMode: " + std::string(SDL_GetError());
    if (0 != SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Fraktaler 3", message.c_str(), nullptr))
    {
      SDL_LogCritical(SDL_LOG_CATEGORY_APPLICATION, "%s", message.c_str());
    }
    SDL_Quit();
    return 1;
  }
  int win_screen_width = mode.w;
  int win_screen_height = mode.h;
#else
  int win_screen_width = 1024;
  int win_screen_height = 576;
#endif

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
  SDL_WindowFlags window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL /* | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI */);
  window = SDL_CreateWindow("Fraktaler 3", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, win_screen_width, win_screen_height, window_flags);
  if (! window)
  {
    const std::string message = "SDL_CreateWindow: " + std::string(SDL_GetError());
    if (0 != SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Fraktaler 3", message.c_str(), nullptr))
    {
      SDL_LogCritical(SDL_LOG_CATEGORY_APPLICATION, "%s", message.c_str());
    }
    SDL_Quit();
    return 1;
  }
  SDL_GLContext gl_context = SDL_GL_CreateContext(window);
  if (! gl_context)
  {
    const std::string message = "SDL_GL_CreateContext: " + std::string(SDL_GetError());
    if (0 != SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Fraktaler 3", message.c_str(), window))
    {
      SDL_LogCritical(SDL_LOG_CATEGORY_APPLICATION, "%s", message.c_str());
    }
    SDL_Quit();
    return 1;
  }
  SDL_GL_MakeCurrent(window, gl_context);
#ifdef __EMSCRIPTEN__
  bool EXT_sRGB = emscripten_webgl_enable_extension(emscripten_webgl_get_current_context(), "EXT_sRGB");
#endif

#ifdef HAVE_GLEW
  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK)
  {
    const std::string message = "glewInit: " + std::string(SDL_GetError());
    if (0 != SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Fraktaler 3", message.c_str(), window))
    {
      SDL_LogCritical(SDL_LOG_CATEGORY_APPLICATION, "%s", message.c_str());
    }
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

  gl_version = (const char *) glGetString(GL_VERSION);
#ifdef __EMSCRIPTEN__
  if (! EXT_sRGB)
  {
    if (is_webgl_1(gl_version))
    {
      const std::string message = "could not enable WebGL 1.0 EXT_sRGB";
      if (0 != SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Fraktaler 3", message.c_str(), window))
      {
        SDL_LogCritical(SDL_LOG_CATEGORY_APPLICATION, "%s", message.c_str());
      }
      SDL_Quit();
      return 1;
    }
  }
  srgb_conversion = 1; // FIXME, should check if framebuffer is linear or sRGB
#endif

  SDL_GL_SetSwapInterval(1);
  SDL_GL_GetDrawableSize(window, &win_pixel_width, &win_pixel_height);

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
  glClearColor(0.5, 0.5, 0.5, 1);
  glClear(GL_COLOR_BUFFER_BIT);

  GLint maximum_texture_size = 0;
  glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maximum_texture_size);

  // setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui::StyleColorsDark();
  ImGui::GetIO().IniFilename = nullptr;
  ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
  ImGui_ImplOpenGL3_Init(glsl_version);

  colours_init();

#ifdef HAVE_FS
  load_dialog = new ImGui::FileBrowser(ImGuiFileBrowserFlags_CloseOnEsc);
  load_dialog->SetTitle("Load...");
  load_dialog->SetTypeFilters({ ".toml", ".exr" });
  save_dialog = new ImGui::FileBrowser(ImGuiFileBrowserFlags_CloseOnEsc | ImGuiFileBrowserFlags_EnterNewFilename | ImGuiFileBrowserFlags_CreateNewDir);
  save_dialog->SetTitle("Save...");
  save_dialog->SetTypeFilters({ ".toml", ".exr" });
#endif

  {
#ifdef __EMSCRIPTEN__
    pref_path = "/fraktaler-3/";
#else
    char *dir = SDL_GetPrefPath("uk.co.mathr", "fraktaler-3");
    if (dir)
    {
      pref_path = dir;
      SDL_free(dir);
    }
    else
    {
      pref_path = "";
    }
#endif
    try
    {
      std::ifstream ifs(pref_path + "gui-settings.toml", std::ios_base::binary);
      ifs.exceptions(std::ifstream::badbit);
      ifs >> window_state;
    }
    catch (const std::exception &e)
    {
      SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "loading GUI settings: %s", e.what());
    }
    try
    {
      if (argc > 1)
      {
        par.load_toml(argv[1]);
      }
      else
      {
        par.load_toml(pref_path + "persistence.toml");
      }
      // FIXME overrides just-loaded parameter
      par.p.image.width = win_pixel_width;
      par.p.image.height = win_pixel_height;
      par.p.image.subsampling = std::min(std::max(par.p.image.subsampling, 1), 32); // FIXME
      clr = colours[par.p.colour_id];
      out = new map((par.p.image.width + par.p.image.subsampling - 1) / par.p.image.subsampling, (par.p.image.height + par.p.image.subsampling - 1) / par.p.image.subsampling, par.p.bailout.iterations);
      dsp = new display_t(clr);
      resize(1, par.p.image.subsampling);
      reset(sta);
    }
    catch (const std::exception &e)
    {
      SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "loading parameter: %s", e.what());
    }
  }
  SDL_AddTimer(one_minute, persistence_timer_callback, nullptr);

  populate_number_type_wisdom();
#ifdef __EMSCRIPTEN__
  emscripten_set_main_loop(main1, 0, false);
  return 0;
#else
  while (! quit)
  {
    main1();
  }
#endif

#ifdef HAVE_FS
  if (load_dialog)
  {
    delete load_dialog;
    load_dialog = nullptr;
  }
  if (save_dialog)
  {
    delete save_dialog;
    save_dialog = nullptr;
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

int main(int argc, char **argv)
{
#ifdef __EMSCRIPTEN__
  EM_ASM(
    FS.mkdir('/fraktaler-3');
    FS.mount(IDBFS, {}, '/fraktaler-3');
    FS.syncfs(true, function (err) {
      assert(! err);
      ccall('main00', 'number');
    });
  );
#else
  return main0(argc, argv);
#endif
}
