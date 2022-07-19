// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

        // z = f(C, Z, c, z)
}

        z2 = real_norm_complexdual(z);
        n++;
        m++;
        iters_ptb++;
      }

      {
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
    }

    // compute output
    const struct complex Z1 = { Zz.x.x, Zz.y.x };
    const struct mat2 J = { Zz.x.dx[0], Zz.x.dx[1], Zz.y.dx[0], Zz.y.dx[1] };
    const struct complex dC = complex_mul_complex_mat2(complex_mul_complex_mat2(Z1, J), config->transform_K);
    const real Z1norm = real_norm_complex(Z1);
    struct complex de = complex_div_real_complex(real_mul_real_real(Z1norm, real_div2_real(real_log_real(Z1norm))), dC);
    float nf = clamp(1.0f - log(log(float_from_real(Z1norm)) / log(float_from_real(config->ER2))) / log(degree), 0.0f, 1.0f);
    float t = atan2(float_from_real(Z1.y), float_from_real(Z1.x)) / 6.283185307179586f;
    t -= floor(t);
    if (bool_lt_real_real(Zz2, config->ER2) || bool_isnan_real(de.x) || bool_isinf_real(de.x) || bool_isnan_real(de.y) || bool_isinf_real(de.y))
    {
      n = config->Iterations;
      nf = 0;
      t = 0;
      de.x = real_from_int(0);
      de.y = real_from_int(0);
    }
    const long k = (j - y0) * config->tile_width + (i - x0);
    /* accumulate colour */
    if (RGB)
    {
      /* colouring algorithm FIXME */
      const float v = clamp(0.75f + 0.125f * 0.5f * log(4.0f * 4.0f * float_from_real(real_norm_complex(de))), 0.0f, 1.0f);
      RGB[3*k+0] = v;
      RGB[3*k+1] = v;
      RGB[3*k+2] = v;
    }
    /* output raw */
    const long Nbias = 1024;
    ulong nn = n + Nbias;
    if (n >= config->Iterations)
    {
      nn = ~((ulong)(0));
    }
    if (N0)
    {
      N0[k] = nn;
    }
    if (N1)
    {
      N1[k] = nn >> 32;
    }
    if (NF)
    {
      NF[k] = nf;
    }
    if (T)
    {
      T[k] = t;
    }
    if (DEX)
    {
      DEX[k] = float_from_real(de.x);
    }
    if (DEY)
    {
      DEY[k] = float_from_real(de.y);
    }
  }
}
