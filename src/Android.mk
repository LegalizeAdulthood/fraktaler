# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021-2024 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)
LOCAL_MODULE    := libgmp
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libgmp.a
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libgmpxx
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libgmpxx.a
LOCAL_STATIC_LIBRARIES := libgmp
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libmpfr
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libmpfr.a
LOCAL_STATIC_LIBRARIES := libgmp
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libdeflate
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libdeflate.a
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libImath
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libImath-3_1.a
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libIex
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libIex-3_2.a
LOCAL_STATIC_LIBRARIES := libImath
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libIlmThread
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libIlmThread-3_2.a
LOCAL_STATIC_LIBRARIES := libIex
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libOpenEXRCore
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libOpenEXRCore-3_2.a
LOCAL_STATIC_LIBRARIES := libdeflate
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libOpenEXR
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libOpenEXR-3_2.a
LOCAL_STATIC_LIBRARIES := libOpenEXRCore libIex libImath libdeflate
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libOpenEXRUtil
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libOpenEXRUtil-3_2.a
LOCAL_STATIC_LIBRARIES := libOpenEXR libOpenEXRCore libIex libImath libdeflate
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libjpeg
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libjpeg.a
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libz
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libz.a
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE    := libpng
LOCAL_SRC_FILES := $(LOCAL_PATH)/$(TARGET_ARCH_ABI)/lib/libpng16.a
LOCAL_STATIC_LIBRARIES := libz
include $(PREBUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE := main
SDL_PATH := ../SDL
LOCAL_C_INCLUDES := \
$(LOCAL_PATH)/$(SDL_PATH)/include \
$(LOCAL_PATH)/$(TARGET_ARCH_ABI)/include \
$(LOCAL_PATH)/$(TARGET_ARCH_ABI)/include/Imath \
$(LOCAL_PATH)/$(TARGET_ARCH_ABI)/include/OpenEXR \
$(LOCAL_PATH)/imgui \
$(LOCAL_PATH)/imgui/backends \
$(LOCAL_PATH)/imgui/misc/cpp \
$(LOCAL_PATH)/imgui-filebrowser \
$(LOCAL_PATH)/implot \

VERSION ?= $(shell test -d $(LOCAL_PATH)/../.git && git describe --always --dirty=+ || (cat $(LOCAL_PATH)/../VERSION.txt | head -n 1))

LOCAL_CFLAGS := \
-fPIC \
-DHAVE_GUI \
-DHAVE_FS \
-DHAVE_EXR=3 \
-DIMGUI_USER_CONFIG="\"f3imconfig.h\"" \
-DFRAKTALER_3_VERSION_STRING="\"$(VERSION)\"" \
-DIMGUI_GIT_VERSION_STRING="\"$(shell cd $(LOCAL_PATH)/imgui && git describe --tags --always --dirty=+)\"" \
-DIMGUI_FILEBROWSER_GIT_VERSION_STRING="\"$(shell cd $(LOCAL_PATH)/imgui-filebrowser && git describe --tags --always --dirty=+)\"" \
-DIMPLOT_GIT_VERSION_STRING="\"$(shell cd $(LOCAL_PATH)/implot && git describe --tags --always --dirty=+)\"" \
-DTOML11_GIT_VERSION_STRING="\"+git\"" \

LOCAL_CPPFLAGS := -std=c++2a
LOCAL_CPP_FEATURES := exceptions rtti
LOCAL_CPP_EXTENSION := .cpp .cc
LOCAL_SRC_FILES := \
batch.cc \
bla.cc \
colour.cc \
display_gles.cc \
engine.cc \
gles2.cc \
glutil.cc \
gui.cc \
histogram.cc \
hybrid.cc \
image_raw.cc \
image_rgb.cc \
jpeg.cc \
main.cc \
opencl.cc \
param.cc \
png.cc \
render.cc \
source.cc \
version.cc \
wisdom.cc \
imgui/imgui.cpp \
imgui/imgui_demo.cpp \
imgui/imgui_draw.cpp \
imgui/imgui_tables.cpp \
imgui/imgui_widgets.cpp \
imgui/backends/imgui_impl_sdl2.cpp \
imgui/backends/imgui_impl_opengl3.cpp \
imgui/misc/cpp/imgui_stdlib.cpp \
implot/implot.cpp \
implot/implot_items.cpp \

LOCAL_SHARED_LIBRARIES := SDL2
LOCAL_STATIC_LIBRARIES := libjpeg libpng libz libOpenEXR libIlmThread libIex libImath libdeflate libmpfr libgmpxx libgmp
LOCAL_LDLIBS := -lGLESv1_CM -lGLESv2 -lGLESv3 -llog
#ifneq ($(TARGET_ARCH_ABI), arm64-v8a)
#LOCAL_LDFLAGS := -Wl,--no-warn-shared-textrel
#endif
include $(BUILD_SHARED_LIBRARY)
