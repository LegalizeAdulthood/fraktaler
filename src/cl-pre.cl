#ifndef NUMBER_TYPE
#error NUMBER_TYPE not defined
#endif

#if NUMBER_TYPE == 2
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

#define MAX_PHASES 64
#define MAX_LEVELS 64

#if NUMBER_TYPE == 1 // float
typedef float real;
#else
#if NUMBER_TYPE == 2 // double
typedef double real;
#else

#if 0 && NUMBER_TYPE == 4 // floatexp

typedef float mantissa;
typedef int exponent;
struct floatexp
{
  mantissa val;
  exponent exp;
};
typedef struct floatexp real;

#else
#if 0 && NUMBER_TYPE == 5 // softfloat

struct softfloat
{
  uint se;
  uint m;
};
typedef struct softfloat real;

#else
#error unsupported NUMBER_TYPE (can handle 1,2)
#endif
#endif
#endif
#endif
#if NUMBER_TYPE == 1 || NUMBER_TYPE == 2

#define float_from_real(x) ((float)(x))
#define real_from_float(x) ((real)(x))
#define real_from_int(x) ((real)(x))
#define real_from_long(x) ((real)(x))
#define real_neg_real(x) (-(x))
#define real_sqr_real(x) ((x)*(x))
#define real_mul2_real(x) ((x)*2.0)
#define real_div2_real(x) ((x)*0.5)
#define real_1div_real(x) (1.0/(x))
#define real_add_real_real(x,y) ((x)+(y))
#define real_sub_real_real(x,y) ((x)-(y))
#define real_mul_real_real(x,y) ((x)*(y))
#define real_div_real_real(x,y) ((x)/(y))
#define real_sqrt_real(x) (sqrt((x)))
#define real_abs_real(x) (fabs((x)))
#define real_exp_real(x) (exp((x)))
#define real_log_real(x) (log((x)))
#define real_sin_real(x) (sin((x)))
#define real_cos_real(x) (cos((x)))
#define real_floor_real(x) (floor((x)))
#define real_atan2_real_real(y,x) (atan2((y),(x)))
#define real_hypot_real_real(x,y) (hypot((x),(y)))
#define real_nextafter_real_real(x,y) (nextafter((x),(y)))
#define real_min_real_real(x,y) (min((x),(y)))
#define bool_gtzero_real(x) ((x)>0.0)
#define bool_gezero_real(x) ((x)>=0.0)
#define bool_ltzero_real(x) ((x)<0.0)
#define bool_lt_real_real(x,y) ((x)<(y))
#define bool_isnan_real(x) (isnan((x)))
#define bool_isinf_real(x) (isinf((x)))
#define real_twopi() (6.283185307179586)

#endif

struct mat2
{
  real a, b, c, d;
};

struct complex
{
  real x, y;
};

struct dual
{
  real x; real dx[2];
};

struct complexdual
{
  struct dual x, y;
};

struct blaR2
{
  struct mat2 A, B;
  real r2;
  long l;
};

struct complex complex_mul_complex_complex(struct complex a, struct complex b)
{
  struct complex r;
  r.x = real_sub_real_real(real_mul_real_real(a.x, b.x), real_mul_real_real(a.y, b.y));
  r.y = real_add_real_real(real_mul_real_real(a.x, b.y), real_mul_real_real(a.y, b.x));
  return r;
}

struct complex complex_div_real_complex(real a, struct complex b)
{
  struct complex r;
  const real d = real_add_real_real(real_sqr_real(b.x), real_sqr_real(b.y));
  r.x = real_mul_real_real(a, real_div_real_real(b.x, d));
  r.y = real_neg_real(real_mul_real_real(a, real_div_real_real(b.y, d)));
  return r;
}

struct dual dual_neg_dual(struct dual a)
{
  struct dual r;
  r.x = real_neg_real(a.x);
  r.dx[0] = real_neg_real(a.dx[0]);
  r.dx[1] = real_neg_real(a.dx[1]);
  return r;
}

struct dual dual_abs_dual(struct dual a)
{
  return bool_ltzero_real(a.x) ? dual_neg_dual(a) : a;
}

struct dual dual_mul_real_dual(real a, struct dual b)
{
  struct dual r;
  r.x = real_mul_real_real(a, b.x);
  r.dx[0] = real_mul_real_real(a, b.dx[0]);
  r.dx[1] = real_mul_real_real(a, b.dx[1]);
  return r;
}

struct dual dual_mul_dual_real(struct dual b, real a)
{
  struct dual r;
  r.x = real_mul_real_real(a, b.x);
  r.dx[0] = real_mul_real_real(a, b.dx[0]);
  r.dx[1] = real_mul_real_real(a, b.dx[1]);
  return r;
}

struct dual dual_add_dual_real(struct dual a, real b)
{
  struct dual r = a;
  r.x = real_add_real_real(r.x, b);
  return r;
}

struct dual dual_add_real_dual(real b, struct dual a)
{
  struct dual r = a;
  r.x = real_add_real_real(r.x, b);
  return r;
}

struct dual dual_sub_dual_real(struct dual a, real b)
{
  struct dual r = a;
  r.x = real_sub_real_real(r.x, b);
  return r;
}

struct dual dual_add_dual_dual(struct dual a, struct dual b)
{
  struct dual r;
  r.x = real_add_real_real(a.x, b.x);
  r.dx[0] = real_add_real_real(a.dx[0], b.dx[0]);
  r.dx[1] = real_add_real_real(a.dx[1], b.dx[1]);
  return r;
}

struct dual dual_sub_dual_dual(struct dual a, struct dual b)
{
  struct dual r;
  r.x = real_sub_real_real(a.x, b.x);
  r.dx[0] = real_sub_real_real(a.dx[0], b.dx[0]);
  r.dx[1] = real_sub_real_real(a.dx[1], b.dx[1]);
  return r;
}

struct dual dual_mul_dual_dual(struct dual a, struct dual b)
{
  struct dual r;
  r.x = real_mul_real_real(a.x, b.x);
  r.dx[0] = real_add_real_real(real_mul_real_real(a.x, b.dx[0]), real_mul_real_real(a.dx[0], b.x));
  r.dx[1] = real_add_real_real(real_mul_real_real(a.x, b.dx[1]), real_mul_real_real(a.dx[1], b.x));
  return r;
}

struct dual dual_exp_dual(struct dual a)
{
  struct dual r;
  r.x = real_exp_real(a.x);
  r.dx[0] = real_mul_real_real(r.x, a.dx[0]);
  r.dx[1] = real_mul_real_real(r.x, a.dx[1]);
  return r;
}

struct dual dual_cos_dual(struct dual a)
{
  struct dual r;
  r.x = real_cos_real(a.x);
  const real d = real_neg_real(real_sin_real(a.x));
  r.dx[0] = real_mul_real_real(d, a.dx[0]);
  r.dx[1] = real_mul_real_real(d, a.dx[1]);
  return r;
}

struct dual dual_sin_dual(struct dual a)
{
  struct dual r;
  r.x = real_sin_real(a.x);
  const real d = real_cos_real(a.x);
  r.dx[0] = real_mul_real_real(d, a.dx[0]);
  r.dx[1] = real_mul_real_real(d, a.dx[1]);
  return r;
}

struct dual dual_diffabs_real_dual(real c, struct dual d)
{
  const real cd = real_add_real_real(c, d.x);
  const struct dual c2d = dual_add_real_dual(real_mul2_real(c), d);
  return bool_gezero_real(c) ? bool_gezero_real(cd) ? d : dual_neg_dual(c2d) : bool_gtzero_real(cd) ? c2d : dual_neg_dual(d);
}

struct complexdual complexdual_add_complex_complexdual(struct complex a, struct complexdual b)
{
  struct complexdual r;
  r.x = dual_add_real_dual(a.x, b.x);
  r.y = dual_add_real_dual(a.y, b.y);
  return r;
}

struct complexdual complexdual_add_complexdual_complexdual(struct complexdual a, struct complexdual b)
{
  struct complexdual r;
  r.x = dual_add_dual_dual(a.x, b.x);
  r.y = dual_add_dual_dual(a.y, b.y);
  return r;
}

struct complexdual complexdual_mul_complexdual_complex(struct complexdual a, struct complex b)
{
  struct complexdual r;
  r.x = dual_sub_dual_dual(dual_mul_dual_real(a.x, b.x), dual_mul_dual_real(a.y, b.y));
  r.y = dual_add_dual_dual(dual_mul_dual_real(a.y, b.x), dual_mul_dual_real(a.x, b.y));
  return r;
}

struct complexdual complexdual_mul_complexdual_complexdual(struct complexdual a, struct complexdual b)
{
  struct complexdual r;
  r.x = dual_sub_dual_dual(dual_mul_dual_dual(a.x, b.x), dual_mul_dual_dual(a.y, b.y));
  r.y = dual_add_dual_dual(dual_mul_dual_dual(a.y, b.x), dual_mul_dual_dual(a.x, b.y));
  return r;
}

struct complexdual complexdual_mul_mat2_complexdual(struct mat2 a, struct complexdual b)
{
  struct complexdual r;
  r.x = dual_add_dual_dual(dual_mul_real_dual(a.a, b.x), dual_mul_real_dual(a.b, b.y));
  r.y = dual_add_dual_dual(dual_mul_real_dual(a.c, b.x), dual_mul_real_dual(a.d, b.y));
  return r;
}

struct complex complex_mul_complex_mat2(struct complex a, struct mat2 b)
{
  struct complex r;
  r.x = real_add_real_real(real_mul_real_real(a.x, b.a), real_mul_real_real(a.y, b.c));
  r.y = real_add_real_real(real_mul_real_real(a.x, b.b), real_mul_real_real(a.y, b.d));
  return r;
}

real real_norm_complex(struct complex a)
{
  return real_add_real_real(real_sqr_real(a.x), real_sqr_real(a.y));
}

real real_norm_complexdual(struct complexdual a)
{
  return real_add_real_real(real_sqr_real(a.x.x), real_sqr_real(a.y.x));
}

real real_arg_complex(struct complex a)
{
  return real_atan2_real_real(a.y, a.x);
}

struct config
{
  long config_size;
  long number_type;
  /* shape */
  long height;
  long width;
  long subframes;
  long frame;
  /* bailout */
  long Iterations;
  real ER2;
  long PerturbIterations;
  /* transform */
  long transform_exponential_map;
  struct mat2 transform_K;
  real pixel_spacing;
  real offset_x;
  real offset_y;
  /* ref layout */
  long number_of_phases;
  long ref_size[MAX_PHASES];
  long ref_start[MAX_PHASES];
  /* bla layout */
  long bla_size[MAX_PHASES];
  long bla_levels[MAX_PHASES];
  long bla_start[MAX_PHASES][MAX_LEVELS];
};

// http://www.burtleburtle.net/bob/hash/integer.html
uint burtle_hash(uint a)
{
  a = (a+0x7ed55d16) + (a<<12);
  a = (a^0xc761c23c) ^ (a>>19);
  a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);
  a = (a+0xfd7046c5) + (a<<3);
  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}

real radical_inverse(long a, const long base)
{
  const real one_minus_epsilon = real_nextafter_real_real(real_from_int(1), real_from_int(0));
  const real base1 = real_1div_real(base);
  long reversed = 0;
  real base1n = real_from_int(1);
  while (a)
  {
    const long next  = a / base;
    const long digit = a - base * next;
    reversed = reversed * base + digit;
    base1n = real_mul_real_real(base1n, base1);
    a = next;
  }
  return real_min_real_real(real_mul_real_real(real_from_long(reversed), base1n), one_minus_epsilon);
}

real wrap(const real v)
{
  return real_sub_real_real(v, real_floor_real(v));
}

real triangle(const real a)
{
  const real b = real_sub_real_real(real_mul2_real(a), real_from_int(1));
  const real c = real_sqrt_real(real_abs_real(b));
  const real e = bool_gtzero_real(b) ? real_sub_real_real(c, real_from_int(1)) : real_sub_real_real(real_from_int(1), c);
  return e;
}

void jitter(const long width, const long height, const long frame, const long i, const long j, const long k, real *x, real *y)
{
  long ix = (frame * height + j) * width + i;
  real h = real_div_real_real(real_from_long(burtle_hash(ix)), real_from_long(0x100000000L));
  *x = triangle(wrap(real_add_real_real(radical_inverse(k, 2), h)));
  *y = triangle(wrap(real_add_real_real(radical_inverse(k, 3), h)));
}

__global const struct blaR2 *lookup_bla(__constant const struct config *config, __global const struct blaR2 *bla, long phase, long m, real z2)
{
  if (m <= 0)
  {
    return 0;
  }
  if (! (m < config->bla_size[phase]))
  {
    return 0;
  }
  __global const struct blaR2 *ret = 0;
  long ix = m - 1;
  for (long level = 0; level < config->bla_levels[phase]; ++level)
  {
    long ixm = (ix << level) + 1;
    long start = config->bla_start[phase][level];
    if (m == ixm && bool_lt_real_real(z2, bla[start + ix].r2))
    {
      ret = &bla[start + ix];
    }
    else
    {
      break;
    }
    ix = ix >> 1;
  }
  return ret;
}

__kernel void fraktaler3
( __constant const struct config *config
, __global const real *ref
, __global const struct blaR2 *bla
/* accumulate linear RGB */
, __global float *RGB
/* output raw data */
, __global uint *N0
, __global uint *N1
, __global float *NF
, __global float *T
, __global float *DEX
, __global float *DEY
, const long subframe
)
{
  // sanity check
  if (config->config_size != sizeof(struct config) || config->number_type != NUMBER_TYPE)
  {
    return;
  }
  const long j = get_global_id(0);
  const long i = get_global_id(1);
  const float degree = 2; // FIXME
  {
    real di, dj;
    jitter(config->width, config->height, config->frame, i, j, subframe, &di, &dj);
    struct dual u0 = { real_add_real_real(real_add_real_real(real_from_long(i), real_from_float(0.5)), di), { real_from_int(1), real_from_int(0) } };
    struct dual v0 = { real_add_real_real(real_add_real_real(real_from_long(j), real_from_float(0.5)), dj), { real_from_int(0), real_from_int(1) } };
    if (config->transform_exponential_map)
    {
      struct dual re = dual_mul_real_dual(real_div_real_real(real_neg_real(real_log_real(real_from_int(2))), real_from_long(config->height)), v0);
      struct dual im = dual_mul_real_dual(real_div_real_real(real_twopi(), real_from_long(config->width)), u0);
      real R = real_div2_real(real_hypot_real_real(real_from_long(config->width), real_from_long(config->height)));
      struct dual r = dual_exp_dual(re);
      struct dual c = dual_cos_dual(im);
      struct dual s = dual_sin_dual(im);
      u0 = dual_mul_real_dual(R, dual_mul_dual_dual(r, c));
      v0 = dual_mul_real_dual(R, dual_mul_dual_dual(r, s));
    }
    else
    {
      u0 = dual_sub_dual_real(u0, real_div2_real(real_from_long(config->width)));
      v0 = dual_sub_dual_real(v0, real_div2_real(real_from_long(config->height)));
    }
    // FIXME should K multiply offset?
    const struct complex C = { ref[config->ref_start[0] + 2], ref[config->ref_start[0] + 3] }; // FIXME
    struct dual cx = dual_add_dual_real(dual_mul_dual_real(u0, config->pixel_spacing), config->offset_x);
    struct dual cy = dual_add_dual_real(dual_mul_dual_real(v0, config->pixel_spacing), config->offset_y);
    struct complexdual c = { cx, cy };
    c = complexdual_mul_mat2_complexdual(config->transform_K, c);
    long phase = 0;
    long m = 0;
    long n = 0;
    long iters_ptb = 0;
    struct complex Z = { ref[config->ref_start[phase] + 0], ref[config->ref_start[phase] + 1] };
    struct complexdual z = { { real_from_int(0), { real_from_int(0), real_from_int(0) } }, { real_from_int(0), { real_from_int(0), real_from_int(0) } } };
    real z2 = real_norm_complexdual(z);
    struct complexdual Zz = complexdual_add_complex_complexdual(Z, z);
    real Zz2 = real_norm_complexdual(Zz);
    while (n < config->Iterations && Zz2 < config->ER2 && iters_ptb < config->PerturbIterations)
    {
      // bla steps
      __global const struct blaR2 *b = 0;
      while (n < config->Iterations && bool_lt_real_real(Zz2, config->ER2) && (b = lookup_bla(config, bla, phase, m, z2)))
      {
        z = complexdual_add_complexdual_complexdual(complexdual_mul_mat2_complexdual(b->A, z), complexdual_mul_mat2_complexdual(b->B, c));
        z2 = real_norm_complexdual(z);
        n += b->l;
        m += b->l;

        // rebase
        if (! (n < config->Iterations && bool_lt_real_real(Zz2, config->ER2) && iters_ptb < config->PerturbIterations))
        {
          break;
        }
        if (! (m < config->ref_size[phase]))
        {
          break;
        }
        struct complex Z = { ref[config->ref_start[phase] + 2 * m], ref[config->ref_start[phase] + 2 * m + 1] };
        Zz = complexdual_add_complex_complexdual(Z, z);
        Zz2 = real_norm_complexdual(Zz);
        if (bool_lt_real_real(Zz2, z2) || m + 1 == config->ref_size[phase])
        {
          z = Zz;
          phase = (phase + m) % config->number_of_phases;
          m = 0;
        }
      }

      // perturbation iteration
      {
        if (! (n < config->Iterations && bool_lt_real_real(Zz2, config->ER2) && iters_ptb < config->PerturbIterations))
        {
          break;
        }
        if (! (m < config->ref_size[phase]))
        {
          break;
        }

        // z = f(C, Z, c, z)
{
