#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>

using namespace std;

const int nvar = 3; // Lorenz chaotic problem !
const int kmax = 4000;
const float dxsav = 0.0001;

const int MAXSTP = 10000; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
const float TINY = 1.0e-30;

float sigma = 10.0;
float betax = 8.0 / 3.0;
float rho = 28.0;

int kount;

float *xx_p = new float[kmax + 1];
float *y_p = new float[nvar * (kmax + 1)];


//----- sign_h
float sign_h(float a, float b)
{
    if (b >= 0.f)
	{
        return (a >= 0.f) ? a : -a;
    } else
	{
        return (a >= 0.f) ? -a : a;
    }
}

//----- derivs
void derivs(const float x, float *y, float *dydx)
{
  dydx[0] = sigma * (y[1] - y[0]);
  dydx[1] = y[0] * (rho - y[2]) - y[1];
  dydx[2] = y[0] * y[1] - betax * y[2];
}



//----- rkck_h
void rkck_h(float *y, float *dydx, const float &x, const float &h, float *yout, float *yerr, int &n,
	void derivs(const float, float*, float*))
{
  const float a2=0.2f, a3=0.3f, a4=0.6f, a5=1.0f, a6=0.875f,
	b21=0.2f, b31=3.0f/40.0f, b32=9.0f/40.0f, b41=0.3f, b42 = -0.9f,
	b43=1.2f, b51 = -11.0f/54.0f, b52=2.5f, b53 = -70.0f/27.0f,
	b54=35.0f/27.0f, b61=1631.0f/55296.0f, b62=175.0f/512.0f,
	b63=575.0f/13824.0f, b64=44275.0f/110592.0f, b65=253.0f/4096.0f,
	c1=37.0f/378.0f, c3=250.0f/621.0f, c4=125.0f/594.0f, c6=512.0f/1771.0f,
	dc1=c1-2825.0f/27648.0f, dc3=c3-18575.0f/48384.0f,
	dc4=c4-13525.0f/55296.0f, dc5 = -277.00f/14336.0f, dc6=c6-0.25f;
		
	int i;

	//int n = y.size();
	float ak2[n], ak3[n], ak4[n], ak5[n], ak6[n], ytemp[n];
	
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	derivs(x+a2*h,ytemp,ak2);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	derivs(x+a3*h,ytemp,ak3);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	derivs(x+a4*h,ytemp,ak4);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	derivs(x+a5*h,ytemp,ak5);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	derivs(x+a6*h,ytemp,ak6);
	for (i=0;i<n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=0;i<n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}


/*


*/

//----- rkqs_h
void rkqs_h(float *y, float *dydx, float &x, const float &htry,
    const float &eps, float *yscal, float &hdid, float &hnext, int n,
    void (*derivs)(const float, float*, float*))
{
    const float SAFETY = 0.9f, PGROW = -0.2f, PSHRNK = -0.25f, ERRCON = 1.89e-4f;

    float errmax, h, htemp, xnew;   

    // int n = y.size();

    h = htry;

    float yerr[n], ytemp[n];

    int i;

    for (;;)
    {
        rkck_h(y, dydx, x, h, ytemp, yerr, n, derivs);

        errmax = 0.0f;

        for (i = 0; i < n; i++) errmax = max(errmax, fabs(yerr[i] / yscal[i]));

        errmax /= eps;

        if (errmax <= 1.0f) break;

        htemp = SAFETY * h * pow(errmax, PSHRNK); // Since here errmax > 1 then htemp < h !

        h = (h >= 0.0f ? max(htemp, 0.1f*h) : min(htemp, 0.1f*h));

        xnew = x + h;

        if (xnew == x) cout <<"stepsize underflow in rkqs" << endl; // xnew is defined just to check if h = 0 or not to raise an error !!!!
    }
    // Step Size Adjustment for Next Call --Conservative growth if error was not very small.--Aggressive growth (5×) if error was small enough.
    if (errmax > ERRCON) hnext = SAFETY * h * pow(errmax, PGROW); // Note: since errmax < 1 then this line will increase the value of h.... so hnext > h.
    else hnext = 5.0f * h;
    x += (hdid = h);
    for (i = 0; i < n; i++) y[i] = ytemp[i];
}


//----- odeint_h
void odeint_h(float *ystart, const float x1, const float x2, const float eps,
	const float h1, const float hmin, int &nok, int &nbad,
	void derivs(const float, float*, float*),
	void rkqs(float*, float*, float&, const float &, const float&, float*, float&, float&,
	int, void (*)(const float, float*, float*)))
{
	//const int MAXSTP = 10000;  // defined globally !
	//const float TINY = 1.0e-30; // defined globally !
	
	float xsav, x, hnext, hdid, h;

	// int nvar=ystart.size(); // Will be set as a global variable !
	float yscal[nvar], y[nvar], dydx[nvar];

	x = x1;
	h = sign_h(h1, x2 - x1);
	nok = nbad = kount = 0;

	int i, nstp;

	for (i= 0; i < nvar; i++) y[i] = ystart[i];

	// Note that xsav is the x of the last save !
	if (kmax > 0) xsav = x - dxsav * 2.0f; // This is a trick to allow us to save the first initial starting point.

	for (nstp = 0;nstp < MAXSTP; nstp++) 
	{
		derivs(x, y, dydx);
		for (i = 0;i < nvar; i++)
			yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;

		if (kmax > 0 && kount < kmax - 1 && fabs(x - xsav) > fabs(dxsav)) 
		{
			for (i = 0; i < nvar; i++) y_p[i * (kmax + 1) + kount] = y[i];
			xx_p[kount++] = x;
			xsav = x;
		}

		if ((x+h-x2)*(x+h-x1) > 0.0f) h=x2-x;

		rkqs(y, dydx, x, h, eps, yscal, hdid, hnext, nvar, derivs);

		if (hdid == h) ++nok; else ++nbad;

		if ((x - x2) * (x2 - x1) >= 0.0f)
		{
			for (i = 0; i < nvar; i++) ystart[i] = y[i];

			if (kmax != 0)
			{
				for (i = 0; i < nvar; i++) y_p[i * (kmax + 1) + kount] = y[i];
				xx_p[kount++] = x;
			}
			return;
		}
		if (fabs(hnext) <= hmin) cout <<"Step size too small in odeint" <<endl;
		h=hnext;
	}
	cout <<"Too many steps in routine odeint" << endl;
}




int main()
{

	float y0[3] = {-8.0f, 8.0f, 27.0f};
	// Note that here x1 and x2 represent time !
	const float x1 = 0.0f; // Initial time !!!!!!!!!!!!!!!!!!!!!!!
	const float x2 = 100.0f; // Final time !!!!!!!!!!!!!!!!!!!!!!!!!
	float h1 = 0.02f; // Step size !!!!!!!!!!!!!!!!!!!!!!!!

	float eps = 1e-7f;   // THIS affects the execution time !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	float hmin = 1e-5f;  //!!!!!!!!!!!!!!!!!!!!!!!!!! To be adjusted for each problem !!!!!!!!!!!!!!!!!!!!!!!!!!

	int nok = 0;
	int nbad = 0;

	auto start = chrono::high_resolution_clock::now();

	odeint_h(y0, x1, x2, eps, h1, hmin, nok, nbad, derivs, rkqs_h);

	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = end - start;
	cout << "Elapsed time: " << elapsed.count() << " seconds" << endl;
	
	cout << "kount = " << kount << endl;
	
	ofstream outfile("outX.csv"); // use pplot.py to plot the result !
  
	outfile << "x,y" << endl;

	for (int j = 0; j < kount; j++)
	{
		int n_row_1 = 0;
		int n_row_2 = 2;
		outfile << y_p[n_row_1 * (kmax+1) + j] << "," << y_p[n_row_2 * (kmax+1) + j] << endl;
	}

	outfile.close();

}
