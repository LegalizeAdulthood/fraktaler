#pragma once
#include "gles2.h"
#define IMGUI_IMPL_OPENGL_LOADER_CUSTOM

#if defined(_WIN32)

#define OPENGL_GLSL_VERSION "#version 330 core"
#define OPENGL_MAJOR 3
#define OPENGL_MINOR 3
#define OPENGL_PROFILE SDL_GL_CONTEXT_PROFILE_CORE
#define OPENGL_FLAGS SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG
// define this while our GL loader does not have glPolygonMode
#define IMGUI_IMPL_OPENGL_ES3
// desktop OpenGL has float buffers in core
//#define WANT_EXT_COLOR_BUFFER_FLOAT
//#define NEED_EXT_COLOR_BUFFER_FLOAT

#else

#define OPENGL_GLSL_VERSION "#version 300 es"
#define OPENGL_MAJOR 3
#define OPENGL_MINOR 0
#define OPENGL_PROFILE SDL_GL_CONTEXT_PROFILE_ES
#define OPENGL_FLAGS 0
#define IMGUI_IMPL_OPENGL_ES3
// we can try to do without, losing per-tile dynamic range and fidelity
// a proper fix might be OpenGL ES 3.1 compute shaders...
#define WANT_EXT_COLOR_BUFFER_FLOAT
//#define NEED_EXT_COLOR_BUFFER_FLOAT

#endif

#define ImDrawIdx unsigned int

