#ifndef DIFF_JAC_H
#define DIFF_JAC_H

//----- jacobn_s
__device__ void jacobn_s(float *y, float *dfdx, float *dfdy, int n)
{
    // Assumes system size is exactly 3 (based on hardcoded entries)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for (int i = 0; i < n; i++)
        dfdx[i] = 0.0f;

    // Fill Jacobian matrix dfdy (row-major order)
    dfdy[0 * n + 0] = -0.013f - 1000.0f * y[2];
    dfdy[0 * n + 1] = 0.0f;
    dfdy[0 * n + 2] = -1000.0f * y[0];

    dfdy[1 * n + 0] = 0.0f;
    dfdy[1 * n + 1] = -2500.0f * y[2];
    dfdy[1 * n + 2] = -2500.0f * y[1];

    dfdy[2 * n + 0] = -0.013f - 1000.0f * y[2];
    dfdy[2 * n + 1] = -2500.0f * y[2];
    dfdy[2 * n + 2] = -1000.0f * y[0] - 2500.0f * y[1];
}


//----- derivs
__device__ void derivs(float *y, float *dydx)
{
    dydx[0] = -0.013f * y[0] - 1000.0f * y[0] * y[2];
    dydx[1] = -2500.0f * y[1] * y[2];
    dydx[2] = -0.013f * y[0] - 1000.0f * y[0] * y[2] - 2500.0f * y[1] * y[2];
}

#endif
