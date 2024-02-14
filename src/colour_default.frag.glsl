uniform float brightness;
uniform float contrast;
uniform float exposure;
uniform float gamma;

vec3 colour(void)
{
  vec2 de = getDE();
  float v = 0.75 + 0.125 * 0.25 * log(4.0 * 4.0 * dot(de, de));
  v = exp2(exposure) * pow(clamp((v + brightness - 0.5) * exp2(contrast) + 0.5, 0.0, 1.0), 1.0 / gamma);
  return vec3(v);
}
