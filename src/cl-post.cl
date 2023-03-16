// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

        // z = f(C, Z, c, z)
}

        z2 = real_norm_complexdual(z);
        n++;
        m++;
        iters_ptb++;
        last_degree = next_degree;
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
        next_degree = config->degree[n % config->number_of_phases];
        if (bool_lt_real_real(Zz2, real_mul_real_real(real_from_int(next_degree), z2)) || m + 1 == config->ref_size[phase])
        {
          z = Zz;
          z2 = Zz2;
          phase = (phase + m) % config->number_of_phases;
          m = 0;
        }
      }
    }

    // compute output
    const struct complex Z1 = { Zz.x.x, Zz.y.x };
    const struct mat2 J = { Zz.x.dx[0], Zz.x.dx[1], Zz.y.dx[0], Zz.y.dx[1] };
    const struct complex dC = complex_mul_complex_mat2(Z1, J);
    const real Z1norm = real_norm_complex(Z1);
    struct complex de = complex_div_real_complex(real_mul_real_real(Z1norm, real_div2_real(real_log_real(Z1norm))), dC);
    float nf = clamp(1.0f - log(float_from_real(real_log_real(Z1norm)) / float_from_real(real_log_real(config->ER2))) / log((float) last_degree), 0.0f, 1.0f);
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
      double nn = ((double)(n)) + ((double)(nf));
      nn /= 1243.0f * 0.5f;
      nn -= floor(nn);
      const float hue = ((float)(nn));
      const float nde = float_from_real(real_norm_complex(de));
      const float chroma_x = cos(2.0f * 3.141592653f * hue) + (nde > 0.0f ? 0.5f * float_from_real(de.x) / sqrt(nde) : 0.0f);
      const float chroma_y = sin(2.0f * 3.141592653f * hue) + (nde > 0.0f ? 0.5f * float_from_real(de.y) / sqrt(nde) : 0.0f);
      const float h = atan2(chroma_y, chroma_x) / (2.0f * 3.141592653f);
      const float s = tanh(0.25f * (chroma_x * chroma_x + chroma_y * chroma_y));
      const float v = clamp(0.75f + 0.125f * 0.5f * log(4.0f * 4.0f * nde), 0.0f, 1.0f);
      float r = 0.0f, g = 0.0f, b = 0.0f;
      if (v > 0.0f)
      {
        hsv2rgb(h, 0.0 * s, v, &r, &g, &b);
      }
      RGB[3*k+0] = r;
      RGB[3*k+1] = g;
      RGB[3*k+2] = b;
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
