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

#if NUMBER_TYPE == 4 // floatexp

typedef float mantissa;
typedef int exponent;
struct floatexp
{
  mantissa val;
  exponent exp;
};
typedef struct floatexp real;
__constant static const exponent LARGE_EXPONENT = sizeof(mantissa) == sizeof(float) ? 126 : 1022;
__constant static const exponent EXP_MIN = sizeof(exponent) == sizeof(int) ? (exponent)(-0x00800000) : (exponent)(-0x0080000000000000L);
__constant static const exponent EXP_MAX = sizeof(exponent) == sizeof(int) ? (exponent)( 0x00800000) : (exponent)( 0x0080000000000000L);

struct floatexp floatexp_from_mantissa_exponent(const mantissa aval, const exponent aexp)
{
  struct floatexp r;
  if (aval == 0)
  {
    r.val = aval;
    r.exp = EXP_MIN;
  }
  else if (isnan(aval))
  {
    r.val = aval;
    r.exp = EXP_MIN;
  }
  else if (isinf(aval))
  {
    r.val = aval;
    r.exp = EXP_MAX;
  }
  else
  {
    int e = 0;
    mantissa f_val = frexp(aval, &e);
    exponent f_exp = e + aexp;
    if (f_exp >= EXP_MAX)
    {
      r.val = f_val / (mantissa)(0);
      r.exp = EXP_MAX;
    }
    else if (f_exp <= EXP_MIN)
    {
      r.val = f_val * (mantissa)(0);
      r.exp = EXP_MIN;
    }
    else
    {
      r.val = f_val;
      r.exp = f_exp;
    }
  }
  return r;
}

struct floatexp floatexp_from_float(const float aval)
{
  return floatexp_from_mantissa_exponent(aval, 0);
}

struct floatexp floatexp_from_int(const int aval)
{
  return floatexp_from_mantissa_exponent(aval, 0);
}

struct floatexp floatexp_from_long(const long aval)
{
  return floatexp_from_mantissa_exponent(aval, 0);
}

mantissa mantissa_from_floatexp(const struct floatexp a)
{
  if (a.exp < -126)
  {
    return a.val * (mantissa)(0);
  }
  if (a.exp > 126)
  {
    return a.val / (mantissa)(0);
  }
  return ldexp((mantissa)(a.val), a.exp);
}

float float_from_floatexp(const struct floatexp a)
{
  return mantissa_from_floatexp(a);
}

struct floatexp floatexp_abs_floatexp(const struct floatexp f)
{
  struct floatexp fe = { fabs(f.val), f.exp };
  return fe;
}

struct floatexp floatexp_neg_floatexp(const struct floatexp f)
{
  struct floatexp fe = { -f.val, f.exp };
  return fe;
}

struct floatexp floatexp_sqr_floatexp(const struct floatexp a)
{
  return floatexp_from_mantissa_exponent(a.val * a.val, a.exp << 1);
}

struct floatexp floatexp_mul_floatexp_floatexp(const struct floatexp a, const struct floatexp b)
{
  return floatexp_from_mantissa_exponent(a.val * b.val, a.exp + b.exp);
}

struct floatexp floatexp_mul2exp_floatexp_exponent(const struct floatexp a, const exponent b)
{
  struct floatexp fe = { a.val, a.exp + b };
  return fe;
}

struct floatexp floatexp_div_floatexp_floatexp(const struct floatexp a, const struct floatexp b)
{
  return floatexp_from_mantissa_exponent(a.val / b.val, a.exp - b.exp);
}

struct floatexp floatexp_div2exp_floatexp_exponent(const struct floatexp a, const exponent b)
{
  struct floatexp fe = { a.val, a.exp - b };
  return fe;
}

struct floatexp floatexp_add_floatexp_floatexp(const struct floatexp a, const struct floatexp b)
{
  if (a.exp > b.exp)
  {
    struct floatexp c = { b.val, b.exp - a.exp };
    return floatexp_from_mantissa_exponent(a.val + mantissa_from_floatexp(c), a.exp);
  }
  else
  {
    struct floatexp c = { a.val, a.exp - b.exp };
    return floatexp_from_mantissa_exponent(mantissa_from_floatexp(c) + b.val, b.exp);
  }
}

struct floatexp floatexp_sub_floatexp_floatexp(const struct floatexp a, const struct floatexp b)
{
  return floatexp_add_floatexp_floatexp(a, floatexp_neg_floatexp(b));
}

int int_cmp_mantissa_mantissa(const mantissa a, const mantissa b)
{
  return (a > b) - (b > a);
}

int int_cmp_floatexp_floatexp(const struct floatexp a, const struct floatexp b)
{
  if (a.val > 0)
  {
    if (b.val <= 0)
    {
      return 1;
    }
    else if (a.exp > b.exp)
    {
      return 1;
    }
    else if (a.exp < b.exp)
    {
      return -1;
    }
    else
    {
      return int_cmp_mantissa_mantissa(a.val, b.val);
    }
  }
  else
  {
    if (b.val > 0)
    {
      return -1;
    }
    else if (a.exp > b.exp)
    {
      return -1;
    }
    else if (a.exp < b.exp)
    {
      return 1;
    }
    else
    {
      return int_cmp_mantissa_mantissa(a.val, b.val);
    }
  }
}

struct floatexp floatexp_min_floatexp_floatexp(const struct floatexp a, const struct floatexp b)
{
  return int_cmp_floatexp_floatexp(a, b) <= 0 ? a : b;
}

bool bool_gezero_floatexp(const struct floatexp a)
{
  return a.val >= 0.0;
}

bool bool_gtzero_floatexp(const struct floatexp a)
{
  return a.val > 0.0;
}

bool bool_ltzero_floatexp(const struct floatexp a)
{
  return a.val < 0;
}

bool bool_lt_floatexp_floatexp(const struct floatexp a, const struct floatexp b)
{
  return int_cmp_floatexp_floatexp(a, b) < 0;
}

struct floatexp floatexp_floor_floatexp(const struct floatexp a)
{
  if (a.exp >= LARGE_EXPONENT)
  {
    // already an integer
    return a;
  }
  else if (a.exp <= -LARGE_EXPONENT)
  {
    // very small
    if (a.val < 0)
    {
      return floatexp_from_int(-1);
    }
    else
    {
      return floatexp_from_int(0);
    }
  }
  else
  {
    // won't under/overflow
    return floatexp_from_mantissa_exponent(floor(mantissa_from_floatexp(a)), 0);
  }
}

struct floatexp floatexp_sqrt_floatexp(const struct floatexp a)
{
  return floatexp_from_mantissa_exponent
    ( sqrt((a.exp & 1) ? 2.0 * a.val : a.val)
    , (a.exp & 1) ? (a.exp - 1) / 2 : a.exp / 2
    );
}

struct floatexp floatexp_log_floatexp(const struct floatexp a)
{
  return floatexp_from_mantissa_exponent(log(a.val) + log(2.0) * a.exp, 0);
}

struct floatexp floatexp_pow_floatexp_ulong(struct floatexp x, ulong n)
{
  switch (n)
  {
    case 0: return floatexp_from_int(1);
    case 1: return x;
    case 2: return floatexp_sqr_floatexp(x);
    case 3: return floatexp_mul_floatexp_floatexp(x, floatexp_sqr_floatexp(x));
    case 4: return floatexp_sqr_floatexp(floatexp_sqr_floatexp(x));
    case 5: return floatexp_mul_floatexp_floatexp(x, floatexp_sqr_floatexp(floatexp_sqr_floatexp(x)));
    case 6: return floatexp_sqr_floatexp(floatexp_mul_floatexp_floatexp(x, floatexp_sqr_floatexp(x)));
    case 7: return floatexp_mul_floatexp_floatexp(x, floatexp_sqr_floatexp(floatexp_mul_floatexp_floatexp(x, floatexp_sqr_floatexp(x))));
    case 8: return floatexp_sqr_floatexp(floatexp_sqr_floatexp(floatexp_sqr_floatexp(x)));
    default:
    {
      struct floatexp y = floatexp_from_int(1);
      while (n > 1)
      {
        if (n & 1)
          y = floatexp_mul_floatexp_floatexp(y, x);
        x = floatexp_sqr_floatexp(x);
        n >>= 1;
      }
      return floatexp_mul_floatexp_floatexp(x, y);
    }
  }
}

struct floatexp floatexp_exp_floatexp(const struct floatexp a)
{
  if (-53 <= a.exp && a.exp <= 8) return floatexp_from_mantissa_exponent(exp(mantissa_from_floatexp(a)), 0);
  if (61 <= a.exp) return a.val > 0.0 ? floatexp_from_mantissa_exponent(a.val / 0.0, 0) : floatexp_from_mantissa_exponent(0.0, 0);
  if (a.exp < -53) return floatexp_from_mantissa_exponent(1.0, 0);
  return floatexp_pow_floatexp_ulong(floatexp_from_mantissa_exponent(exp(a.val), 0), ((ulong)(1)) << a.exp);
}

struct floatexp floatexp_sin_floatexp(const struct floatexp a)
{
  mantissa y = mantissa_from_floatexp(a);
  if (isinf(y))
  {
    return floatexp_from_mantissa_exponent(0.0/0.0, 0);
  }
  if (isnan(y))
  {
    return floatexp_from_mantissa_exponent(y, 0);
  }
  if (y == 0) // FIXME denormalized numbers lose precision
  {
    return a;
  }
  return floatexp_from_mantissa_exponent(sin(y), 0);
}

struct floatexp floatexp_cos_floatexp(const struct floatexp a)
{
  mantissa y = mantissa_from_floatexp(a);
  if (isinf(y))
  {
    return floatexp_from_mantissa_exponent(0.0/0.0, 0);
  }
  if (isnan(y))
  {
    return a;
  }
  return floatexp_from_mantissa_exponent(cos(y), 0);
}

struct floatexp floatexp_diffabs_floatexp_floatexp(const struct floatexp c, const struct floatexp d)
{
  const struct floatexp cd = floatexp_add_floatexp_floatexp(c, d);
  const struct floatexp c2d = floatexp_add_floatexp_floatexp(floatexp_mul2exp_floatexp_exponent(c, 1), d);
  return c.val >= 0.0 ? cd.val >= 0.0 ? d : floatexp_neg_floatexp(c2d) : cd.val > 0.0 ? c2d : floatexp_neg_floatexp(d);
}

struct floatexp floatexp_hypot_floatexp_floatexp(const struct floatexp x, const struct floatexp y)
{
  return floatexp_sqrt_floatexp(floatexp_add_floatexp_floatexp(floatexp_sqr_floatexp(x), floatexp_sqr_floatexp(y)));
}

struct floatexp floatexp_atan2_floatexp_floatexp(struct floatexp y, struct floatexp x)
{
  struct floatexp z = floatexp_hypot_floatexp_floatexp(y, x);
  x = floatexp_div_floatexp_floatexp(x, z);
  y = floatexp_div_floatexp_floatexp(y, z);
  return floatexp_from_mantissa_exponent(atan2(mantissa_from_floatexp(y), mantissa_from_floatexp(x)), 0);
}

struct floatexp floatexp_nextafter_floatexp_floatexp(const struct floatexp x, const struct floatexp y)
{
  return floatexp_from_mantissa_exponent(nextafter(x.val, mantissa_from_floatexp(floatexp_div2exp_floatexp_exponent(y, x.exp))), x.exp);
}

bool bool_isinf_floatexp(const struct floatexp x)
{
  return isinf(x.val);
}

bool bool_isnan_floatexp(const struct floatexp x)
{
  return isnan(x.val);
}

#define float_from_real(x) (float_from_floatexp((x)))
#define real_from_float(x) (floatexp_from_float((x)))
#define real_from_int(x) (floatexp_from_int((x)))
#define real_from_long(x) (floatexp_from_long((x)))
#define real_neg_real(x) (floatexp_neg_floatexp((x)))
#define real_sqr_real(x) (floatexp_sqr_floatexp((x)))
#define real_mul2_real(x) (floatexp_mul2exp_floatexp_exponent((x), 1))
#define real_div2_real(x) (floatexp_div2exp_floatexp_exponent((x), 1))
#define real_1div_real(x) (floatexp_div_floatexp_floatexp(floatexp_from_int(1), (x)))
#define real_add_real_real(x,y) (floatexp_add_floatexp_floatexp((x),(y)))
#define real_sub_real_real(x,y) (floatexp_sub_floatexp_floatexp((x),(y)))
#define real_mul_real_real(x,y) (floatexp_mul_floatexp_floatexp((x),(y)))
#define real_div_real_real(x,y) (floatexp_div_floatexp_floatexp((x),(y)))
#define real_sqrt_real(x) (floatexp_sqrt_floatexp((x)))
#define real_abs_real(x) (floatexp_abs_floatexp((x)))
#define real_exp_real(x) (floatexp_exp_floatexp((x)))
#define real_log_real(x) (floatexp_log_floatexp((x)))
#define real_sin_real(x) (floatexp_sin_floatexp((x)))
#define real_cos_real(x) (floatexp_cos_floatexp((x)))
#define real_floor_real(x) (floatexp_floor_floatexp((x)))
#define real_atan2_real_real(y,x) (floatexp_atan2_floatexp_floatexp((y),(x)))
#define real_hypot_real_real(x,y) (floatexp_hypot_floatexp_floatexp((x),(y)))
#define real_nextafter_real_real(x,y) (floatexp_nextafter_floatexp_floatexp((x),(y)))
#define real_min_real_real(x,y) (floatexp_min_floatexp_floatexp((x),(y)))
#define bool_gtzero_real(x) (bool_gtzero_floatexp((x)))
#define bool_gezero_real(x) (bool_gezero_floatexp((x)))
#define bool_ltzero_real(x) (bool_ltzero_floatexp((x)))
#define bool_lt_real_real(x,y) (bool_lt_floatexp_floatexp((x),(y)))
#define bool_isnan_real(x) (bool_isnan_floatexp((x)))
#define bool_isinf_real(x) (bool_isinf_floatexp((x)))
#define real_twopi() (floatexp_from_mantissa_exponent(6.283185307179586, 0))

#else
#if 0 && NUMBER_TYPE == 5 // softfloat

struct softfloat
{
  uint se;
  uint m;
};
typedef struct softfloat real;

#else
#error unsupported NUMBER_TYPE (can handle 1,2,4)
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
  const real base1 = real_1div_real(real_from_long(base));
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
  const long j = get_global_id(0);
  const long i = get_global_id(1);
  // sanity check
  if (config->config_size != sizeof(struct config) || config->number_type != NUMBER_TYPE)
  {
    const long k = j * config->width + i;
    if (RGB)
    {
      float v = 254.0 / 255.0;
      RGB[3*k+0] = v;
      RGB[3*k+1] = v;
      RGB[3*k+2] = v;
    }
    return;
  }
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
    while (n < config->Iterations && bool_lt_real_real(Zz2, config->ER2) && iters_ptb < config->PerturbIterations)
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
