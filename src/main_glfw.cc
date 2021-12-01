// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <chrono>
#include <thread>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
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
    complex<floatexp> *Zfe = new complex<floatexp>[par.MaximumReferenceIterations];
    const count_t M = reference(Zfe, par.MaximumReferenceIterations, par.Cx, par.Cy, &progress[0], running);
    progress[1] = 1;
    render(out, par, Zoom, M, Zfe, &progress[2], running);
    delete[] Zfe;
  }
  else if (Zoom > e10(1, 300))
  {
    complex<long double> *Zld = new complex<long double>[par.MaximumReferenceIterations];
    const count_t M = reference(Zld, par.MaximumReferenceIterations, par.Cx, par.Cy, &progress[0], running);
    progress[1] = 1;
    render(out, par, (long double)(Zoom), M, Zld, &progress[2], running);
    delete[] Zld;
  }
  else if (Zoom > e10(1, 30))
  {
    complex<double> *Zd = new complex<double>[par.MaximumReferenceIterations];
    const count_t M = reference(Zd, par.MaximumReferenceIterations, par.Cx, par.Cy, &progress[0], running);
    progress[1] = 1;
    render(out, par, double(Zoom), M, Zd, &progress[2], running);
    delete[] Zd;
  }
  else
  {
    complex<float> *Zf = new complex<float>[par.MaximumReferenceIterations];
    const count_t M = reference(Zf, par.MaximumReferenceIterations, par.Cx, par.Cy, &progress[0], running);
    progress[1] = 1;
    render(out, par, float(Zoom), M, Zf, &progress[2], running);
    delete[] Zf;
  }
  *ended = true;
}


static void glfw_error_callback(int error, const char* description)
{
  std::fprintf(stderr, "glfw error %d: %s\n", error, description);
}

static void APIENTRY opengl_debug_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *message, const void *userParam)
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

static GLFWwindow *create_window(int major, int minor)
{
  const size_t output_window_width = 1024;
  const size_t output_window_height = 576;
  glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, major);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, minor);
  if ((major > 3) || (major == 3 && minor >= 2))
  {
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  }
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
  glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
#if (GLFW_VERSION_MAJOR > 3) || (GLFW_VERSION_MAJOR == 3 && GLFW_VERSION_MINOR >= 3)
  glfwWindowHint(GLFW_SCALE_TO_MONITOR, GLFW_TRUE);
#endif
  return glfwCreateWindow(output_window_width, output_window_height, "Fraktaler 3", 0, 0);
}

int main_window(int argc, char **argv)
{
  // initialize OpenGL
  glfwSetErrorCallback(glfw_error_callback);
  if (! glfwInit())
  {
    return 1;
  }
  // try OpenGL 4 first
  const char *glsl_version = "#version 400 core";
  GLFWwindow *window = create_window(4, 0);
  if (! window)
  {
    // try OpenGL 3
    glsl_version = "#version 330 core";
    window = create_window(3, 3);
    if (! window)
    {
      // give up
      return 1;
    }
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);
  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK)
  {
    std::fprintf(stderr, "glew error\n");
    return 1;
  }
  glGetError(); // discard common error from glew
  if (glDebugMessageCallback)
  {
    glDebugMessageCallback(opengl_debug_callback, nullptr);
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
  }

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
  glClearColor(0.5, 0.5, 0.5, 1);
  glClear(GL_COLOR_BUFFER_BIT);

  GLint maximum_texture_size = 0;
  glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maximum_texture_size);

  param par;
  mpfr_init2(par.Cx, 24);
  mpfr_init2(par.Cy, 24);
  mpfr_set_d(par.Cx, 0, MPFR_RNDN);
  mpfr_set_d(par.Cy, 0, MPFR_RNDN);
  par.Zoom = 1;
  par.Iterations = 1000;
  par.ReferencePeriod = 1;
  par.MaximumReferenceIterations = par.ReferencePeriod + 1;
  par.PerturbIterations = 1000;
  par.ExponentialMap = false;
  par.ZoomOutSequence = false;
  par.Channels = Channels_default;
  par.Stem = "fraktaler-3.exr";
  par.Width = 1920;
  par.Height = 1080;

  map out(par.Width, par.Height, par.Iterations);

  display dsp;

  bool quit = false;
  bool restart = false;
  while (! quit)
  {
    {
      progress_t progress[4] = { 0, 0, 0, 0 };
      bool running = true;
      bool ended = false;
      std::thread bg(main_window_thread, std::ref(out), std::cref(par), &progress[0], &running, &ended);
      int tick = 0;
      while (! quit && ! ended)
      {
        glfwPollEvents();
        if (glfwWindowShouldClose(window))
        {
          quit = true;
          running = false;
        }
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
      bg.join();
    }
    dsp.resize(out);
    dsp.upload_raw(out);
    dsp.colourize();
    while (! quit && ! restart)
    {
      glfwWaitEvents();
      if (glfwWindowShouldClose(window))
      {
        quit = true;
      }
      int display_w, display_h;
      glfwGetFramebufferSize(window, &display_w, &display_h);
      dsp.draw(display_w, display_h);
      glfwSwapBuffers(window);
    }
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  return 0;
}
