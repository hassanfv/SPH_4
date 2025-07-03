#ifndef DIFF_JAC_H
#define DIFF_JAC_H


//----- derivs
__device__ void derivs(float *y, float *dydx)
{
    dydx[0] = -10.0f * y[0] + y[1] * y[2];
    dydx[1] = -0.1f * y[1] + y[0] * y[2];
    dydx[2] = -100.0f * y[2] + y[3];
    dydx[3] = -0.5f * y[3] + y[4] * y[0];
    dydx[4] = -y[4] + y[5] * y[6];
    dydx[5] = -0.01f * y[5] + y[6];
    dydx[6] = -20.0f * y[6] + y[7];
    dydx[7] = -0.2f * y[7] + y[8] * y[1];
    dydx[8] = -y[8] + y[9] * y[0];
    dydx[9] = -5.0f * y[9] + y[0] * y[2];
}


//----- jacobn_s
__device__ void jacobn_s(float *y, float *dfdx, float *dfdy, int n)
{
    for (int i = 0; i < n; i++)
        dfdx[i] = 0.0f;

    for (int i = 0; i < n * n; i++)
        dfdy[i] = 0.0f;

    dfdy[0 * n + 0] = -10.0f;
    dfdy[0 * n + 1] = y[2];
    dfdy[0 * n + 2] = y[1];

    dfdy[1 * n + 0] = y[2];
    dfdy[1 * n + 1] = -0.1f;
    dfdy[1 * n + 2] = y[0];

    dfdy[2 * n + 2] = -100.0f;
    dfdy[2 * n + 3] = 1.0f;

    dfdy[3 * n + 0] = y[4];
    dfdy[3 * n + 3] = -0.5f;
    dfdy[3 * n + 4] = y[0];

    dfdy[4 * n + 4] = -1.0f;
    dfdy[4 * n + 5] = y[6];
    dfdy[4 * n + 6] = y[5];

    dfdy[5 * n + 5] = -0.01f;
    dfdy[5 * n + 6] = 1.0f;

    dfdy[6 * n + 6] = -20.0f;
    dfdy[6 * n + 7] = 1.0f;

    dfdy[7 * n + 1] = y[8];
    dfdy[7 * n + 7] = -0.2f;
    dfdy[7 * n + 8] = y[1];

    dfdy[8 * n + 0] = y[9];
    dfdy[8 * n + 8] = -1.0f;
    dfdy[8 * n + 9] = y[0];

    dfdy[9 * n + 0] = y[2];
    dfdy[9 * n + 2] = y[0];
    dfdy[9 * n + 9] = -5.0f;
}



#endif
