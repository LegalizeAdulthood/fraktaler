        // z = f(C, Z, c, z)
}

        z2 = double_norm_complexdual(z);
        n++;
        m++;
        iters_ptb++;
      }

      {
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
    }

    // compute output
    const struct complex Z1 = { Zz.x.x, Zz.y.x };
    const struct mat2 J = { Zz.x.dx[0], Zz.x.dx[1], Zz.y.dx[0], Zz.y.dx[1] };
    const struct complex dC = complex_mul_complex_mat2(complex_mul_complex_mat2(Z1, J), config->transform_K);
    const double Z1norm = double_norm_complex(Z1);
    struct complex de = complex_div_double_complex(0.5 * Z1norm * log(Z1norm), dC);
    double nf = 0;
    double t = 0;
    float nf = min(max(1 - log(log(Z1norm) / log(config->ER2)) / log(degree), 0), 1);
    float t = double_complex_arg(Z1) / 6.283185307179586;
    t -= floor(t);
    if (Zz2 < config->ER2 || isnan(de.x) || isinf(de.x) || isnan(de.y) || isinf(de.y))
    {
      n = config->Iterations;
      nf = 0;
      t = 0;
      de.x = 0;
      de.y = 0;
    }
    const long k = j * config->width + i;
    /* accumulate colour */
    if (RGB)
    {
      /* colouring algorithm FIXME */
      const float v = clamp(0.75 + 0.125 * 0.5 * log(4.0 * 4.0 * double_norm_complex(de)), 0.0, 1.0);
      if (subframe == 0)
      {
        RGB[3*k+0] = 0;
        RGB[3*k+1] = 0;
        RGB[3*k+2] = 0;
      }
      RGB[3*k+0] += v;
      RGB[3*k+1] += v;
      RGB[3*k+2] += v;
      if (subframe == config->subframes - 1)
      {
        RGB[3*k+0] /= config->subframes;
        RGB[3*k+1] /= config->subframes;
        RGB[3*k+2] /= config->subframes;
      }
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
      DEX[k] = de.x;
    }
    if (DEY)
    {
      DEY[k] = de.y;
    }
  }
}
