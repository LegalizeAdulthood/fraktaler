// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <chrono>
#include <thread>

#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <SDL2/SDL_opengles2.h>
#else
#include <SDL2/SDL_opengl.h>
#endif
#include <mpfr.h>

#include "display.h"
#include "floatexp.h"
#include "main.h"
#include "map.h"
#include "param.h"
#include "reference.h"
#include "render.h"

void main_window_thread(map &out, const param &par, progress_t *progress, bool *running, bool *ended)
{
  floatexp Zoom = par.Zoom;
  if (Zoom > e10(1, 4900))
  {
std::cerr << "floatexp" << std::endl;
    complex<floatexp> *Zfe = new complex<floatexp>[par.MaximumReferenceIterations];
    const count_t M = reference(Zfe, par.MaximumReferenceIterations, par.Cx, par.Cy, &progress[0], running);
    progress[1] = 1;
    render(out, par, Zoom, M, Zfe, &progress[2], running);
    delete[] Zfe;
  }
  else if (Zoom > e10(1, 300))
  {
std::cerr << "long double" << std::endl;
    complex<long double> *Zld = new complex<long double>[par.MaximumReferenceIterations];
    const count_t M = reference(Zld, par.MaximumReferenceIterations, par.Cx, par.Cy, &progress[0], running);
    progress[1] = 1;
    render(out, par, (long double)(Zoom), M, Zld, &progress[2], running);
    delete[] Zld;
  }
  else if (Zoom > e10(1, 30))
  {
std::cerr << "double" << std::endl;
    complex<double> *Zd = new complex<double>[par.MaximumReferenceIterations];
    const count_t M = reference(Zd, par.MaximumReferenceIterations, par.Cx, par.Cy, &progress[0], running);
    progress[1] = 1;
    render(out, par, double(Zoom), M, Zd, &progress[2], running);
    delete[] Zd;
  }
  else
  {
std::cerr << "float" << std::endl;
    complex<float> *Zf = new complex<float>[par.MaximumReferenceIterations];
    const count_t M = reference(Zf, par.MaximumReferenceIterations, par.Cx, par.Cy, &progress[0], running);
    progress[1] = 1;
    render(out, par, float(Zoom), M, Zf, &progress[2], running);
    delete[] Zf;
  }
  *ended = true;
#if 0
  SDL_Event e;
  e.type = SDL_USEREVENT;
  e.user.code = 1;
  e.user.data1 = nullptr;
  e.user.data2 = nullptr;
  SDL_PushEvent(&e);
#endif
}

static void glfw_error_callback(int error, const char* description)
{
  std::fprintf(stderr, "glfw error %d: %s\n", error, description);
}

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

// rendering state machine
bool quit = false;
bool running = false;
bool restart = false;
bool ended = true;

// zoom by mouse drag
bool drag = false;
double drag_start_x = 0;
double drag_start_y = 0;

// last mouse coordinates
int mouse_x = 0;
int mouse_y = 0;

// imgui state
bool show_demo_window = true;
bool show_another_window = false;
ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

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
            STOP
            zoom(par, cx, cy, d, false);
            restart = true;
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

        case SDLK_KP_PLUS:
          STOP
          if (shift)
          {
            if (par.PerturbIterations < count_t(1) << 48)
            {
              par.PerturbIterations <<= 1;
            }
          }
          else
          {
            if (par.Iterations < count_t(1) << 48)
            {
              par.Iterations <<= 1;
              par.MaximumReferenceIterations <<= 1;
            }
          }
          restart = true;
          break;
        case SDLK_KP_MINUS:
          STOP
          if (shift)
          {
            if (par.PerturbIterations > count_t(1) << 8)
            {
              par.PerturbIterations >>= 1;
            }
          }
          else
          {
            if (par.Iterations > count_t(1) << 8)
            {
              par.Iterations >>= 1;
              par.MaximumReferenceIterations >>= 1;
            }
          }
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

        default:
          break;
      }
      break;
  }
}

void display_background(SDL_Window *window, display &dsp)
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

void display_gui(SDL_Window *window, display &dsp)
{
  // Start the Dear ImGui frame
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplSDL2_NewFrame();
  ImGui::NewFrame();

  // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
  if (show_demo_window)
    ImGui::ShowDemoWindow(&show_demo_window);

  // 2. Show a simple window that we create ourselves. We use a Begin/End pair to created a named window.
  {
    static float f = 0.0f;
    static int counter = 0;

    ImGui::Begin("Hello, world!");                          // Create a window called "Hello, world!" and append into it.

    ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
    ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
    ImGui::Checkbox("Another Window", &show_another_window);

    ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
    ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color

    if (ImGui::Button("Button"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
      counter++;
    ImGui::SameLine();
    ImGui::Text("counter = %d", counter);

    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
    ImGui::End();
  }

  // 3. Show another simple window.
  if (show_another_window)
  {
    ImGui::Begin("Another Window", &show_another_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
    ImGui::Text("Hello from another window!");
    if (ImGui::Button("Close Me"))
      show_another_window = false;
    ImGui::End();
  }

  // rendering
  ImGui::Render();
  ImGuiIO& io = ImGui::GetIO();
  glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
  glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
  glClear(GL_COLOR_BUFFER_BIT);
  display_background(window, dsp);
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  SDL_GL_SwapWindow(window);
}

int main_window(int argc, char **argv)
{
  const coord_t win_width = 1024;
  const coord_t win_height = 576;

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
  SDL_Window* window = SDL_CreateWindow("Fraktaler 3", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, win_width, win_height, window_flags);
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

  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK)
  {
    std::cerr << argv[0] << ": error: glewInit" << std::endl;
    SDL_Quit();
    return 1;
  }
  glGetError(); // discard common error from glew
  if (glDebugMessageCallback)
  {
    glDebugMessageCallback(opengl_debug_callback, nullptr);
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
  }

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

  param par;
  mpfr_init2(par.Cx, 24);
  mpfr_init2(par.Cy, 24);
  mpfr_set_d(par.Cx, 0, MPFR_RNDN);
  mpfr_set_d(par.Cy, 0, MPFR_RNDN);
  par.Zoom = 1;
  par.Iterations = 1024;
  par.ReferencePeriod = 0;
  par.MaximumReferenceIterations = par.Iterations;
  par.PerturbIterations = 1024;
  par.ExponentialMap = false;
  par.ZoomOutSequence = false;
  par.Channels = Channels_default;
  par.Stem = "fraktaler-3.exr";
  par.Width = 1920;
  par.Height = 1080;

  map out(par.Width, par.Height, par.Iterations);

  display dsp;
  dsp.resize(out);

  while (! quit)
  {
    std::cerr << "render" << std::endl;
    {
      progress_t progress[4] = { 0, 0, 0, 0 };
      running = true;
      ended = false;
      restart = false;
      std::thread bg(main_window_thread, std::ref(out), std::cref(par), &progress[0], &running, &ended);
      int tick = 0;
      while (! quit && ! ended)
      {
        display_gui(window, dsp);
        SDL_Event e;
        while (SDL_PollEvent(&e))
        {
          ImGui_ImplSDL2_ProcessEvent(&e);
          if (! (io.WantCaptureMouse || io.WantCaptureKeyboard))
          {
            handle_event(window, e, par);
          }
        }
        if (! quit && ! ended)
        {
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
          if (tick == 0)
          {
            std::cerr
              << "Reference["     << std::setw(3) << int(progress[0] * 100) << "%] "
              << "Frame["         << std::setw(3) << int(progress[1] * 100) << "%] "
              << "Approximation[" << std::setw(3) << int(progress[2] * 100) << "%] "
              << "Pixels["        << std::setw(3) << int(progress[3] * 100) << "%] "
              << "\r";
          }
          tick = (tick + 1) % 500;
        }
      }
      bg.join();
    }
    if (running) // was not interrupted
    {
      std::cerr << "display" << std::endl;
      dsp.upload_raw(out);
      dsp.colourize();
    }
    while (! quit && ! restart)
    {
      display_gui(window, dsp);
      SDL_Event e;
      if (SDL_WaitEvent(&e))
      {
        ImGui_ImplSDL2_ProcessEvent(&e);
        if (! (io.WantCaptureMouse || io.WantCaptureKeyboard))
        {
          handle_event(window, e, par);
        }
      }
    }
    std::cerr << "quit or restart" << std::endl;
  }

  // cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();
  SDL_GL_DeleteContext(gl_context);
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
