#ifndef STIFFLIBS_H
#define STIFFLIBS_H

#define CHECK_CUDA(call) { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << ": " \
                  << cudaGetErrorString(err) << std::endl; \
        exit(EXIT_FAILURE); \
    } \
}

//----- acquire_slot
__device__ int acquire_slot(int* slot_status, int N) {
    for (int i = 0; i < N; ++i) {
        if (atomicCAS(&slot_status[i], 0, 1) == 0) {
            return i; // got slot i
        }
    }
    return -1; // no slot available right now
}


//----- release_slot
__device__ void release_slot(int* slot_status, int slot_id) {
    atomicExch(&slot_status[slot_id], 0);
}


//----- ludcmp
__device__ void ludcmp(float *a, int *indx, int n)
{
    const float TINY = 1.0e-20f;
    int i, imax = 0, j, k;
    float big, dum, sum, temp;

    // Allocate vv on stack (private to thread)
    float vv[10];  // change 10 to appropriate value !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    float d = 1.0f;

    for (i = 0; i < n; i++) 
    {
        big = 0.0f;
        for (j = 0; j < n; j++) 
        {
            temp = fabsf(a[i * n + j]);
            if (temp > big) big = temp;
        }
        if (big == 0.0f) 
        {
            // No printf in production device code, consider setting a flag or return NaN later!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            return;
        }
        vv[i] = 1.0f / big;
    }

    for (j = 0; j < n; j++) 
    {
        for (i = 0; i < j; i++) 
        {
            sum = a[i * n + j];
            for (k = 0; k < i; k++)
                sum -= a[i * n + k] * a[k * n + j];
            a[i * n + j] = sum;
        }

        big = 0.0f;
        for (i = j; i < n; i++) 
        {
            sum = a[i * n + j];
            for (k = 0; k < j; k++)
                sum -= a[i * n + k] * a[k * n + j];
            a[i * n + j] = sum;
            dum = vv[i] * fabsf(sum);
            if (dum >= big) 
            {
                big = dum;
                imax = i;
            }
        }

        if (j != imax) 
        {
            for (k = 0; k < n; k++) 
            {
                dum = a[imax * n + k];
                a[imax * n + k] = a[j * n + k];
                a[j * n + k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }

        indx[j] = imax;
        if (a[j * n + j] == 0.0f) a[j * n + j] = TINY;

        if (j != n - 1) 
        {
            dum = 1.0f / a[j * n + j];
            for (i = j + 1; i < n; i++)
                a[i * n + j] *= dum;
        }
    }
}


//----- lubksb
__device__ void lubksb(float *a, int *indx, float *b, int n)
{
    int i, ii = 0, ip, j;
    float sum;

    for (i = 0; i < n; i++) 
    {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];

        if (ii != 0) {
            for (j = ii - 1; j < i; j++) 
            {
                sum -= a[i * n + j] * b[j];
            }
        } 
        else if (sum != 0.0f) 
        {
            ii = i + 1;
        }

        b[i] = sum;
    }

    for (i = n - 1; i >= 0; i--) 
    {
        sum = b[i];
        for (j = i + 1; j < n; j++) 
        {
            sum -= a[i * n + j] * b[j];
        }
        b[i] = sum / a[i * n + i];
    }
}


//----- stiff_d
__device__ void stiff_d(float *y, float *dydx, float &x, const float htry, // we use &x because x will b updated step by step !
                        const float eps, float *yscal, float &hnext, int n, // n is nvar !
                        float *a, float *dfdy, int *indx, float *dfdx, float *dysav,
                        float *err, float *ysav, float *g1, float *g2, float *g3, float *g4, int idx) // idx is particle index !
{
    const float SAFETY=0.9f,GROW=1.5f,PGROW=-0.25f,SHRNK=0.5f;
    const float PSHRNK=-1.0f/3.0f, ERRCON=0.1296f;
    const int MAXTRY=40;
    const float GAM=0.5f, A21=2.0f, A31=48.0f/25.0f, A32=6.0f/25.0f;
    const float C21=-8.0f, C31=372.0f/25.0f, C32=12.0f/5.0f;
    const float C41=-112.0f/125.0f, C42=-54.0f/125.0f, C43=-2.0f/5.0f;
    const float B1=19.0f/9.0f, B2=0.5f, B3=25.0f/108.0f, B4=125.0f/108.0f;
    const float E1=17.0f/54.0f, E2=7.0f/36.0f, E3=0.0f, E4=125.0f/108.0f;
    const float C1X=0.5f, C2X=-1.5f, C3X=121.0f/50.0f, C4X=29.0f/250.0f;
    const float A2X=1.0f, A3X=3.0f/5.0f;

    int i, j, jtry;
    float errmax, h;
    float xsav = x;

    for (i = 0; i < n; i++)
    {
        ysav[i] = y[i];
        dysav[i] = dydx[i];
    }

    jacobn_s(ysav, dfdx, dfdy, n); // n is nvar.
    h = htry;

    for (jtry = 0; jtry < MAXTRY; jtry++) // Note: ysav is never overwritten so in every loop it has the initial values.
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++) a[i * n + j] = -dfdy[i * n + j];
            a[i * n + i] += 1.0f / (GAM * h);
        }

        ludcmp(a, indx, n);

        for (i = 0; i < n; i++) g1[i] = dysav[i] + h * C1X * dfdx[i];
        lubksb(a, indx, g1, n);
        for (i = 0; i < n; i++) y[i] = ysav[i] + A21 * g1[i];

        x = xsav + A2X * h;
        derivs(y, dydx);

        for (i = 0; i < n; i++) g2[i] = dydx[i] + h * C2X * dfdx[i] + C21 * g1[i] / h;
        lubksb(a, indx, g2, n);
        for (i = 0; i < n; i++) y[i] = ysav[i] + A31 * g1[i] + A32 * g2[i];

        x = xsav + A3X * h;
        derivs(y, dydx);

        for (i = 0; i < n; i++) g3[i] = dydx[i] + h * C3X * dfdx[i] + (C31 * g1[i] + C32 * g2[i]) / h;
        lubksb(a, indx, g3, n);

        for (i = 0; i < n; i++) g4[i] = dydx[i] + h * C4X * dfdx[i] + (C41 * g1[i] + C42 * g2[i] + C43 * g3[i]) / h;
        lubksb(a, indx, g4, n);

        for (i = 0; i < n; i++) {
            y[i] = ysav[i] + B1 * g1[i] + B2 * g2[i] + B3 * g3[i] + B4 * g4[i];
            err[i] = E1 * g1[i] + E2 * g2[i] + E3 * g3[i] + E4 * g4[i];
        }

        x = xsav + h;
        errmax = 0.0f;
        for (i = 0; i < n; i++)
        {
            float e = fabsf(err[i] / yscal[i]);
            if (e > errmax) errmax = e;
        }

        errmax /= eps;
        if (errmax <= 1.0f)
        {
            hnext = (errmax > ERRCON ? SAFETY * h * powf(errmax, PGROW) : GROW * h);
            return;
        } else
        {
            hnext = SAFETY * h * powf(errmax, PSHRNK);
            h = (h >= 0.0f ? fmaxf(hnext, SHRNK * h) : fminf(hnext, SHRNK * h));
        }
    }
    // Instead of printf (not allowed in __device__ unless using printf from device), set a flag or use debug mode
}


//----- sign_d
__device__ float sign_d(float a, float b)
{
    if (b >= 0.f)
    {
        return (a >= 0.f) ? a : -a;
    }
    else
    {
        return (a >= 0.f) ? -a : a;
    }
}


//----- odeint_d
// I can use d_y as both ystart in odeint and y in stiff. They are both the same.
__global__ void odeint_d(float *y_c, float x1, float x2, float eps, float htry, float hmin, int nvar, float *dydx_c, float *yscal_c, // y was ystart !
                              float *a_c, float *dfdy_c, int *indx_c, float *dfdx_c, float *dysav_c, float *err_c, float *ysav_c,
                              float *g1_c, float *g2_c, float *g3_c, float *g4_c, int *slot_status, int N_batch)
{

  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  
  if (idx < N_batch)
  {
    
    float *y = &y_c[idx * nvar]; // Good idea, so that we do not need to do idx * nvar everywhere we have y ! Now y points to idx * nvar location !
  
    //********************************************************
    //********************************************************
    // We need to find a free slot for the particle idx !!!!!!
    //********************************************************
    //********************************************************
    //-------- find a free slot -----------
    int slot_id = idx; //-1;
    /*
    while (slot_id == -1) {
        slot_id = acquire_slot(slot_status, N_con); // N_con --> N_concurrent.
    }
    */
    // now all of them are pointers of the start of the slot so working with them is extremely simpler now !
    float *a     = &a_c[slot_id * nvar * nvar]; // 3D
    float *dfdy  = &dfdy_c[slot_id * nvar * nvar]; // 3D
    float *yscal = &yscal_c[slot_id * nvar];
    float *dydx  = &dydx_c[slot_id * nvar];
    int *indx  = &indx_c[slot_id * nvar];
    float *dfdx  = &dfdx_c[slot_id * nvar];
    float *dysav = &dysav_c[slot_id * nvar];
    float *err   = &err_c[slot_id * nvar];
    float *ysav  = &ysav_c[slot_id * nvar];
    float *g1    = &g1_c[slot_id * nvar];
    float *g2    = &g2_c[slot_id * nvar];
    float *g3    = &g3_c[slot_id * nvar];
    float *g4    = &g4_c[slot_id * nvar];
    //-------------------------------------

    const int MAXSTP = 10000; // May need to decrease it !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    const float TINY = 1.0e-30f;

    float x = x1;
    float h = sign_d(htry, x2 - x1);
    float hnext;

    // in y each row represents one SPH particle and each column represent an ionization fraction or ionic abundance !
    for (int nstp = 0; nstp < MAXSTP; nstp++) 
    {
      derivs(y, dydx); // x--> scalar,  y--> N_part * nvar,   dydx--> N_part * nvar. Note: dydx is the starting pointer of the free slot!
      for (int i = 0; i < nvar; i++)
        yscal[i] = fabsf(y[i]) + fabsf(dydx[i] * h) + TINY; // Only y contains the whole 1,000,000 particles !

      if ((x + h - x2) * (x + h - x1) > 0.0f)
        h = x2 - x;

      stiff_d(y, dydx, x, h, eps, yscal, hnext, nvar,
              a, dfdy, indx, dfdx, dysav, err, ysav, g1, g2, g3, g4, idx); // n_s is n_free_slot. It will be used by 'a' and 'dfdy' in stiff !

      if ((x - x2) * (x2 - x1) >= 0.0f)
      {
        // We should free the slot here.
        //release_slot(slot_status, slot_id);
        return;
      }

      if (fabsf(hnext) <= hmin) 
      {
        for (int i = 0; i < nvar; i++)
          y[i] = nanf("small_step"); // it used to be y[idx * nvar + i] but since we defined *y = &y_c[idx * nvar] we can only use y[i] instead. Great !!!
        return;
      }

      h = hnext;
    }

    for (int i = 0; i < nvar; i++)
      y[i] = nanf("too_many_steps"); // it used to be y[idx * nvar + i] but since we defined *y = &y_c[idx * nvar] we can only use y[i] instead. Great !!!
  }
}



#endif
