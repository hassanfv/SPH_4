#ifndef DIFF_JAC_H
#define DIFF_JAC_H


__device__ void derivs(float *y, float *dydx, int n = 50)
{
    for (int i = 0; i < n; i++) {
        float a = 0.1f + 0.01f * i;
        float b = 0.2f + 0.005f * i;
        int ip1 = (i + 1) % n;
        int ip2 = (i + 2) % n;
        dydx[i] = -a * y[i] + b * y[ip1] * __sinf(y[ip2]);
    }
}

__device__ void jacobn_s(float *y, float *dfdx, float *dfdy, int n = 50)
{
    for (int i = 0; i < n; i++)
        dfdx[i] = 0.0f;

    for (int i = 0; i < n * n; i++)
        dfdy[i] = 0.0f;

    for (int i = 0; i < n; i++) {
        float a = 0.1f + 0.01f * i;
        float b = 0.2f + 0.005f * i;
        int ip1 = (i + 1) % n;
        int ip2 = (i + 2) % n;

        dfdy[i * n + i] = -a;
        dfdy[i * n + ip1] = b * __sinf(y[ip2]);
        dfdy[i * n + ip2] = b * y[ip1] * __cosf(y[ip2]);
    }
}




#endif
