#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define MAX_PHASES 64
#define MAX_LEVELS 64

struct mat2
{
  double a, b, c, d;
};

struct complex
{
  double x, y;
};

struct dual
{
  double x; double dx[2];
};

struct complexdual
{
  struct dual x, y;
};

struct blaR2
{
  struct mat2 A, B;
  double r2;
  long l;
};

struct complex complex_mul_complex_complex(struct complex a, struct complex b)
{
  struct complex r;
  r.x = a.x * b.x - a.y * b.y;
  r.y = a.x * b.y + a.y * b.x;
  return r;
}

struct complex complex_div_double_complex(double a, struct complex b)
{
  struct complex r;
  const double d = b.x * b.x + b.y * b.y;
  r.x = a * b.x / d;
  r.y = -a * b.y / d;
  return r;
}

struct dual dual_neg_dual(struct dual a)
{
  struct dual r;
  r.x = -a.x;
  r.dx[0] = -a.dx[0];
  r.dx[1] = -a.dx[1];
  return r;
}

struct dual dual_abs_dual(struct dual a)
{
  return a.x < 0 ? dual_neg_dual(a) : a;
}

struct dual dual_mul_double_dual(double a, struct dual b)
{
  struct dual r;
  r.x = a * b.x;
  r.dx[0] = a * b.dx[0];
  r.dx[1] = a * b.dx[1];
  return r;
}

struct dual dual_mul_dual_double(struct dual b, double a)
{
  struct dual r;
  r.x = a * b.x;
  r.dx[0] = a * b.dx[0];
  r.dx[1] = a * b.dx[1];
  return r;
}

struct dual dual_add_dual_double(struct dual a, double b)
{
  struct dual r = a;
  r.x += b;
  return r;
}

struct dual dual_add_double_dual(double b, struct dual a)
{
  struct dual r = a;
  r.x += b;
  return r;
}

struct dual dual_sub_dual_double(struct dual a, double b)
{
  struct dual r = a;
  r.x -= b;
  return r;
}

struct dual dual_add_dual_dual(struct dual a, struct dual b)
{
  struct dual r;
  r.x = a.x + b.x;
  r.dx[0] = a.dx[0] + b.dx[0];
  r.dx[1] = a.dx[1] + b.dx[1];
  return r;
}

struct dual dual_sub_dual_dual(struct dual a, struct dual b)
{
  struct dual r;
  r.x = a.x - b.x;
  r.dx[0] = a.dx[0] - b.dx[0];
  r.dx[1] = a.dx[1] - b.dx[1];
  return r;
}

struct dual dual_mul_dual_dual(struct dual a, struct dual b)
{
  struct dual r;
  r.x = a.x *  b.x;
  r.dx[0] = a.x * b.dx[0] + a.dx[0] * b.x;
  r.dx[1] = a.x * b.dx[1] + a.dx[1] * b.x;
  return r;
}

struct dual dual_exp_dual(struct dual a)
{
  struct dual r;
  r.x = exp(a.x);
  r.dx[0] = r.x * a.dx[0];
  r.dx[1] = r.x * a.dx[1];
  return r;
}

struct dual dual_cos_dual(struct dual a)
{
  struct dual r;
  r.x = cos(a.x);
  const double d = -sin(a.x);
  r.dx[0] = d * a.dx[0];
  r.dx[1] = d * a.dx[1];
  return r;
}

struct dual dual_sin_dual(struct dual a)
{
  struct dual r;
  r.x = sin(a.x);
  const double d = cos(a.x);
  r.dx[0] = d * a.dx[0];
  r.dx[1] = d * a.dx[1];
  return r;
}

struct dual dual_diffabs_double_dual(double c, struct dual d)
{
  const double cd = c + d.x;
  const struct dual c2d = dual_add_double_dual(2 * c, d);
  return c >= 0 ? cd >= 0 ? d : dual_neg_dual(c2d) : cd > 0 ? c2d : dual_neg_dual(d);
}

struct complexdual complexdual_add_complex_complexdual(struct complex a, struct complexdual b)
{
  struct complexdual r;
  r.x = dual_add_double_dual(a.x, b.x);
  r.y = dual_add_double_dual(a.y, b.y);
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
  r.x = dual_sub_dual_dual(dual_mul_dual_double(a.x, b.x), dual_mul_dual_double(a.y, b.y));
  r.y = dual_add_dual_dual(dual_mul_dual_double(a.y, b.x), dual_mul_dual_double(a.x, b.y));
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
  r.x = dual_add_dual_dual(dual_mul_double_dual(a.a, b.x), dual_mul_double_dual(a.b, b.y));
  r.y = dual_add_dual_dual(dual_mul_double_dual(a.c, b.x), dual_mul_double_dual(a.d, b.y));
  return r;
}

struct complex complex_mul_complex_mat2(struct complex a, struct mat2 b)
{
  struct complex r;
  r.x = a.x * b.a + a.y * b.c;
  r.y = a.x * b.b + a.y * b.d;
  return r;
}

double double_norm_complex(struct complex a)
{
  return a.x * a.x + a.y * a.y;
}

double double_norm_complexdual(struct complexdual a)
{
  return a.x.x * a.x.x + a.y.x * a.y.x;
}

struct config
{
  /* shape */
  long height;
  long width;
  long subframes;
  long frame;
  /* bailout */
  long Iterations;
  double ER2;
  long PerturbIterations;
  /* transform */
  long transform_exponential_map;
  struct mat2 transform_K;
  double pixel_spacing;
  double offset_x;
  double offset_y;
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

double radical_inverse(long a, const long base)
{
  const double one_minus_epsilon = 0.99999999999999989;
  const double base1 = 1.0 / base;
  long reversed = 0;
  double base1n = 1;
  while (a)
  {
    const long next  = a / base;
    const long digit = a - base * next;
    reversed = reversed * base + digit;
    base1n *= base1;
    a = next;
  }
  return min(reversed * base1n, one_minus_epsilon);
}

double wrap(const double v)
{
  return v - floor(v);
}

double triangle(const double a)
{
  const double b = a * 2 - 1;
  const double c = sqrt(fabs(b));
  const double e = b > 0 ? c - 1 : 1 - c;
  return e;
}

void jitter(const long width, const long height, const long frame, const long i, const long j, const long k, double *x, double *y)
{
  long ix = (frame * height + j) * width + i;
  double h = burtle_hash(ix) / (double) 0x100000000L;
  *x = triangle(wrap(radical_inverse(k, 2) + h));
  *y = triangle(wrap(radical_inverse(k, 3) + h));
}

__global const struct blaR2 *lookup_bla(__constant const struct config *config, __global const struct blaR2 *bla, long phase, long m, double z2)
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
    if (m == ixm && z2 < bla[start + ix].r2)
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
, __global const double *ref
, __global const struct blaR2 *bla
, __global float *grey
, const long subframe
)
{
  const long j = get_global_id(0);
  const long i = get_global_id(1);
  const float degree = 2; // FIXME
  {
    double di, dj;
    jitter(config->width, config->height, config->frame, i, j, subframe, &di, &dj);
    struct dual u0 = { i + di, { 1, 0 } };
    struct dual v0 = { j + dj, { 0, 1 } };
    if (config->transform_exponential_map)
    {
      struct dual re = dual_mul_double_dual(-0.6931471805599453 / config->height, v0); // log 2
      struct dual im = dual_mul_double_dual(6.283185307179586 / config->width, u0); // 2 pi
      double R = 0.5 * hypot((double) config->width, (double) config->height);
      struct dual r = dual_exp_dual(re);
      struct dual c = dual_cos_dual(im);
      struct dual s = dual_sin_dual(im);
      u0 = dual_mul_double_dual(R, dual_mul_dual_dual(r, c));
      v0 = dual_mul_double_dual(R, dual_mul_dual_dual(r, s));
    }
    else
    {
      u0 = dual_sub_dual_double(u0, config->width / 2.0);
      v0 = dual_sub_dual_double(v0, config->height / 2.0);
    }
    // FIXME should K multiply offset?
    const struct complex C = { ref[config->ref_start[0] + 2], ref[config->ref_start[0] + 3] }; // FIXME
    struct dual cx = dual_add_dual_double(dual_mul_dual_double(u0, config->pixel_spacing), config->offset_x);
    struct dual cy = dual_add_dual_double(dual_mul_dual_double(v0, config->pixel_spacing), config->offset_y);
    struct complexdual c = { cx, cy };
    c = complexdual_mul_mat2_complexdual(config->transform_K, c);
    long phase = 0;
    long m = 0;
    long n = 0;
    long iters_ptb = 0;
    struct complex Z = { ref[config->ref_start[phase] + 0], ref[config->ref_start[phase] + 1] };
    struct complexdual z = { { 0, { 0, 0 } }, { 0, { 0, 0 } } };
    double z2 = double_norm_complexdual(z);
    struct complexdual Zz = complexdual_add_complex_complexdual(Z, z);
    double Zz2 = double_norm_complexdual(Zz);
    while (n < config->Iterations && Zz2 < config->ER2 && iters_ptb < config->PerturbIterations)
    {
      // bla steps
      __global const struct blaR2 *b = 0;
      while (n < config->Iterations && Zz2 < config->ER2 && (b = lookup_bla(config, bla, phase, m, z2)))
      {
        z = complexdual_add_complexdual_complexdual(complexdual_mul_mat2_complexdual(b->A, z), complexdual_mul_mat2_complexdual(b->B, c));
        z2 = double_norm_complexdual(z);
        n += b->l;
        m += b->l;

        // rebase
        if (! (n < config->Iterations && Zz2 < config->ER2 && iters_ptb < config->PerturbIterations))
        {
          break;
        }
        if (! (m < config->ref_size[phase]))
        {
          break;
        }
        struct complex Z = { ref[config->ref_start[phase] + 2 * m], ref[config->ref_start[phase] + 2 * m + 1] };
        Zz = complexdual_add_complex_complexdual(Z, z);
        Zz2 = double_norm_complexdual(Zz);
        if (Zz2 < z2 || m + 1 == config->ref_size[phase])
        {
          z = Zz;
          phase = (phase + m) % config->number_of_phases;
          m = 0;
        }
      }

      // perturbation iteration
      {
        if (! (n < config->Iterations && Zz2 < config->ER2 && iters_ptb < config->PerturbIterations))
        {
          break;
        }
        if (! (m < config->ref_size[phase]))
        {
          break;
        }

        // z = f(C, Z, c, z)
{
