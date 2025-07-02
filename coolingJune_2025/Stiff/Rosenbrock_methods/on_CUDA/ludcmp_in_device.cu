


__device__ void ludcmp(float *a, int *indx, int n)
{
    const float TINY = 1.0e-20f;
    int i, imax = 0, j, k;
    float big, dum, sum, temp;

    // Allocate vv on stack (private to thread)
    float vv[3];  // change 3 to appropriate value !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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





