

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
                              float *g1_c, float *g2_c, float *g3_c, float *g4_c, int *slot_status, int N_con, int Npart)
{

  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  
  if (idx < Npart)
  {
    
    float *y = &y_c[idx * nvar]; // Good idea, so that we do not need to do idx * nvar everywhere we have y ! Now y points to idx * nvar location !
  
    //********************************************************
    //********************************************************
    // We need to find a free slot for the particle idx !!!!!!
    //********************************************************
    //********************************************************
    //-------- find a free slot -----------
    int slot_id = -1;
    while (slot_id == -1) {
        slot_id = acquire_slot(slot_status, N_con); // N_con --> N_concurrent.
    }
    // now all of them are pointers of the start of the slot so working with them is extremely simpler now !
    float *a     = &a_c[slot_id * nvar * nvar]; // 3D
    float *dfdy  = &dfdy_c[slot_id * nvar * nvar]; // 3D
    float *yscal = &yscal_c[slot_id * nvar];
    float *dydx  = &dydx_c[slot_id * nvar];
    float *indx  = &indx_c[slot_id * nvar];
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
      derivs(x, y, dydx, nvar, idx); // x--> scalar,  y--> N_part * nvar,   dydx--> N_part * nvar. Note: dydx is the starting pointer of the free slot!
      for (int i = 0; i < nvar; i++)
        yscal[i] = fabsf(y[i]) + fabsf(dydx[i] * h) + TINY; // Only y contains the whole 1,000,000 particles !

      if ((x + h - x2) * (x + h - x1) > 0.0f)
        h = x2 - x;

      stiff_d(y, dydx, x, h, eps, yscal, hnext, nvar,
              a, dfdy, indx, dfdx, dysav, err, ysav, g1, g2, g3, g4, idx); // n_s is n_free_slot. It will be used by 'a' and 'dfdy' in stiff !

      if ((x - x2) * (x2 - x1) >= 0.0f)
      {
        return;
      }

      if (fabsf(hnext) <= hmin) 
      {
        for (int i = 0; i < nvar; i++)
          y[i] = nanf("small_step"); // it used to be y[idx * nvar + i] but since we defined *y = &y_c[idx * nvar] we can only use y[i] instead. Great !!!
          // We should free the slot here.
          release_slot(slot_status, slot_id);
        return;
      }

      h = hnext;
    }

    for (int i = 0; i < nvar; i++)
      y[i] = nanf("too_many_steps"); // it used to be y[idx * nvar + i] but since we defined *y = &y_c[idx * nvar] we can only use y[i] instead. Great !!!
  }
}


