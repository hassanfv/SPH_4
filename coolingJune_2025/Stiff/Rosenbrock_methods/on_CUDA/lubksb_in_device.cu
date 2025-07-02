


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

