LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE := main

SDL_PATH := ../SDL

LOCAL_C_INCLUDES := $(LOCAL_PATH)/$(SDL_PATH)/include

# Add your application source files here...
LOCAL_SRC_FILES :=  bla.cc bla.h colour.cc colour.h colour_monochrome.h colour_rainbow.h complex.h display.h display_gl.cc display_gl.h dual.h engine.cc engine.h floatexp.h formula.cc formula.h formula_burningship.h formula_mandelbrot.h fraktaler-3-source.7z.cc glutil.cc glutil.h gui.cc gui.h map.cc map.h gui.cc matrix.h param.cc param.h source.cc source.h stats.cc stats.h types.h version.cc version.h

LOCAL_SHARED_LIBRARIES := SDL2

LOCAL_LDLIBS := -lGLESv1_CM -lGLESv2 -llog

include $(BUILD_SHARED_LIBRARY)
