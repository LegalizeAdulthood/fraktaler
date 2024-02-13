// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2024 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

precision highp float;

in vec2 Internal_coord;
out vec4 Internal_frag_colour;

// input tile
uniform highp usampler2D Internal_N0;
uniform highp usampler2D Internal_N1;
uniform        sampler2D Internal_NF;
uniform        sampler2D Internal_T;
uniform        sampler2D Internal_DEX;
uniform        sampler2D Internal_DEY;
uniform highp usampler2D Internal_BLA;
uniform highp usampler2D Internal_PTB;
uniform bool Internal_have_N0;
uniform bool Internal_have_N1;
uniform bool Internal_have_NF;
uniform bool Internal_have_T;
uniform bool Internal_have_DEX;
uniform bool Internal_have_DEY;
uniform bool Internal_have_BLA;
uniform bool Internal_have_PTB;

// coordinates
uniform uvec3 Internal_tile;
uniform uvec2 Internal_tile_size;
uniform uvec2 Internal_image_size;
uniform uint Internal_frame;
uniform float Internal_zoom_log_2;
uniform float Internal_time;

// BEGIN public API
vec3  colour(void); // implement this
uint  getN0 (void) { return Internal_have_N0  ? texture(Internal_N0,  Internal_coord).x : 0u; }
uint  getN1 (void) { return Internal_have_N1  ? texture(Internal_N1,  Internal_coord).x : 0u; }
float getNF (void) { return Internal_have_NF  ? texture(Internal_NF,  Internal_coord).x : 0.0; }
float getT  (void) { return Internal_have_T   ? texture(Internal_T,   Internal_coord).x : 0.0; }
float getDEX(void) { return Internal_have_DEX ? texture(Internal_DEX, Internal_coord).x : 0.0; }
float getDEY(void) { return Internal_have_DEY ? texture(Internal_DEY, Internal_coord).x : 0.0; }
uint  getBLA(void) { return Internal_have_BLA ? texture(Internal_BLA, Internal_coord).x : 0u; }
uint  getPTB(void) { return Internal_have_PTB ? texture(Internal_PTB, Internal_coord).x : 0u; }
vec2  getDE (void) { return vec2(getDEX(), getDEY()); }
vec2  getCoord(void);
ivec2 getImageSize(void) { return ivec2(Internal_image_size); }
float getTime(void) { return Internal_time; }
float getZoomLog2(void) { return Internal_zoom_log_2; }
float linear2sRGB(float c);
vec3 linear2sRGB(vec3 c);
float sRGB2linear(float c);
vec3 sRGB2linear(vec3 c);
vec3 hsv2sRGB(vec3 c);
vec3 hsv2rgb(vec3 c);
// END public API

// https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
float linear2sRGB(float c)
{
  c = clamp(c, 0.0, 1.0);
  const float a = 0.055;
  if (c <= 0.0031308)
    return 12.92 * c;
  else
    return (1.0 + a) * float(pow(c, 1.0 / 2.4)) - a;
}

vec3 linear2sRGB(vec3 c)
{
  return vec3(linear2sRGB(c.x), linear2sRGB(c.y), linear2sRGB(c.z));
}

// https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
float sRGB2linear(float c)
{
  c = clamp(c, 0.0, 1.0);
  const float a = 0.055;
  if (c <= 0.04045)
    return c / 12.92;
  else
    return float(pow((c + a) / (1.0 + a), 2.4));
}

vec3 sRGB2linear(vec3 c)
{
  return vec3(sRGB2linear(c.x), sRGB2linear(c.y), sRGB2linear(c.z));
}

// http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl
vec3 hsv2sRGB(vec3 c)
{
  vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  vec3 rgb = c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
  return rgb;
}

vec3 hsv2rgb(vec3 c)
{
  return sRGB2linear(hsv2sRGB(c));
}

// http://www.burtleburtle.net/bob/hash/integer.html
uint burtle_hash(uint a)
{
  a = (a+0x7ed55d16u) + (a<<12);
  a = (a^0xc761c23cu) ^ (a>>19);
  a = (a+0x165667b1u) + (a<<5);
  a = (a+0xd3a2646cu) ^ (a<<9);
  a = (a+0xfd7046c5u) + (a<<3);
  a = (a^0xb55a4f09u) ^ (a>>16);
  return a;
}

float radical_inverse(uint a, uint base)
{
  float one_minus_epsilon = 0.9999999;
  float base1 = 1.0 / float(base);
  uint reversed = 0u;
  float base1n = 1.0;
  while (a > 0u)
  {
    uint next  = a / base;
    uint digit = a - base * next;
    reversed = reversed * base + digit;
    base1n = base1n * base1;
    a = next;
  }
  return min(float(reversed) * base1n, one_minus_epsilon);
}

float wrap(float v)
{
  return v - floor(v);
}

float triangle(float a)
{
  float b = 2.0 * a - 1.0;
  float c = sqrt(abs(b));
  float e = b > 0.0 ? c - 1.0 : 1.0 - c;
  return e;
}

vec2 jitter(uvec3 p)
{
  uint ix = (Internal_frame * Internal_image_size.y + p.y) * Internal_image_size.x + p.x;
  float h = float(burtle_hash(ix)) * 2.3283064365387e-10;
  return vec2
    ( triangle(wrap(radical_inverse(p.z, 2u) + h))
    , triangle(wrap(radical_inverse(p.z, 3u) + h))
    );
}

vec2 getCoord(void) { vec2 p = (vec2(Internal_tile.xy) + Internal_coord) * vec2(Internal_tile_size); return p + jitter(uvec3(p, Internal_tile.z)); }

void main(void)
{
  Internal_frag_colour = vec4(colour(), 1.0);
}

#line 0 1
