// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2024 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "gles2.h"

struct colour;
struct patom;
struct tile;

colour *colour_new();
void colour_delete(colour *c);

void colour_set_image_size(colour *c, int width, int height);
void colour_set_zoom_log_2(colour *c, float zoom_log_2);
void colour_set_time(colour *c, float time);

void colour_set_program(struct colour *u, GLuint program);

bool colour_display(struct colour *u, bool show_gui);
bool colour_display_late(struct colour *u);

std::string colour_get_shader(colour *u);
bool colour_set_shader(colour *u, std::string source);

std::vector<std::map<std::string, patom>> colour_get_uniforms(colour *u);
bool colour_set_uniforms(colour *u, std::vector<std::map<std::string, patom>> &r);

void colour_upload(const colour *u);
void colour_tile(colour *c, int x, int y, int subframe, tile *data);
