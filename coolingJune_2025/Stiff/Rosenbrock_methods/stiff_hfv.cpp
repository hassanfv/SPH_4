#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

using namespace std;


const int nvar = 3; // see derivs_s !
const int kmax = 4000;
const float dxsav = 0.0001f;

const int MAXSTP = 10000; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
const float TINY = 1.0e-30;

int kount;

float *x_p; // Declare globally
float *y_p;


//----- lubksb
void lubksb(float *a, int *indx, float *b, int n) // void lubksb(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b)
{
	int i,ii=0,ip,j;
	float sum;

	//int n=a.nrows();
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= a[i * n + j]*b[j];
		else if (sum != 0.0)
			ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i * n + j]*b[j];
		b[i]=sum/a[i * n + i];
	}
}


//----- ludcmp
void ludcmp(float *a, int *indx, float &d, int n) // void NR::ludcmp(Mat_IO_DP &a, Vec_O_INT &indx, DP &d)
{
	const float TINY=1.0e-20;
	int i,imax,j,k;
	float big,dum,sum,temp;

	//int n=a.nrows();
	float vv[n]; // Vec_DP vv(n);
	d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i * n + j])) > big) big=temp;
		if (big == 0.0) cout << "Singular matrix in routine ludcmp" << endl;
		vv[i]=1.0f/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i * n + j];
			for (k=0;k<i;k++) sum -= a[i * n + k]*a[k * n + j];
			a[i * n + j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i * n + j];
			for (k=0;k<j;k++) sum -= a[i * n + k]*a[k * n + j];
			a[i * n + j]=sum;
			if ((dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax * n + k];
				a[imax * n + k]=a[j * n + k];
				a[j * n + k]=dum;
			}
			d = -d;
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j * n + j] == 0.0) a[j * n + j]=TINY;
		if (j != n-1) {
			dum=1.0f/(a[j * n + j]);
			for (i=j+1;i<n;i++) a[i * n + j] *= dum;
		}
	}
}


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


//----- jacobn_s
void jacobn_s(const float x, float *y, float *dfdx, float *dfdy, int n) // Note: dfdy is a Matrix !
{
	int i;

	//int n=y.size();
	for (i=0;i<n;i++) dfdx[i]=0.0f;
	dfdy[0 * n + 0] = -0.013f-1000.0f*y[2];
	dfdy[0 * n + 1] = 0.0f;
	dfdy[0 * n + 2] = -1000.0f*y[0];
	dfdy[1 * n + 0] = 0.0f;
	dfdy[1 * n + 1] = -2500.0f*y[2];
	dfdy[1 * n + 2] = -2500.0f*y[1];
	dfdy[2 * n + 0] = -0.013f-1000.0f*y[2];
	dfdy[2 * n + 1] = -2500.0f*y[2];
	dfdy[2 * n + 2] = -1000.0f*y[0]-2500.0f*y[1];
}

//----- derivs
void derivs(const float x, float *y, float *dydx)
{
	dydx[0] = -0.013f*y[0]-1000.0f*y[0]*y[2];
	dydx[1] = -2500.0f*y[1]*y[2];
	dydx[2] = -0.013f*y[0]-1000.0f*y[0]*y[2]-2500.0f*y[1]*y[2];
}


//----- stiff
void stiff(float *y, float *dydx, float &x, const float htry, // htry Not by reference !!! Is that OK ????
	const float eps, float *yscal, float &hdid, float &hnext, int n, // Note that I added n (n is actually nvar)!!
	void derivs(const float, float *, float *))
{
	const float SAFETY=0.9,GROW=1.5,PGROW= -0.25,SHRNK=0.5;
	const float PSHRNK=(-1.0/3.0),ERRCON=0.1296;
	const int MAXTRY=40;
	const float GAM=1.0/2.0,A21=2.0,A31=48.0/25.0,A32=6.0/25.0,C21= -8.0,
		C31=372.0/25.0,C32=12.0/5.0,C41=(-112.0/125.0),
		C42=(-54.0/125.0),C43=(-2.0/5.0),B1=19.0/9.0,B2=1.0/2.0,
		B3=25.0/108.0,B4=125.0/108.0,E1=17.0/54.0,E2=7.0/36.0,E3=0.0,
		E4=125.0/108.0,C1X=1.0/2.0,C2X=(-3.0/2.0),C3X=(121.0/50.0),
		C4X=(29.0/250.0),A2X=1.0,A3X=3.0/5.0;
	int i,j,jtry;
	float d,errmax,h,xsav;

	//int n=y.size();
	float *a = new float[n * n]; // Mat_DP a(n,n);
	float *dfdy = new float[n * n]; // Mat_DP dfdy(n,n);
	
	int *indx = new int[n]; // Vec_INT indx(n);
	
  // Vec_DP dfdx(n),dysav(n),err(n),ysav(n),g1(n),g2(n),g3(n),g4(n);
	float *dfdx = new float[n]; // Vec_DP dfdx(n)
	float *dysav = new float[n]; // Vec_DP dysav(n)
	float *err = new float[n]; // Vec_DP err(n)
	float *ysav = new float[n]; // Vec_DP ysav(n)
	float *g1 = new float[n]; // Vec_DP g1(n)
	float *g2 = new float[n]; // Vec_DP g2(n)
	float *g3 = new float[n]; // Vec_DP g3(n)
	float *g4 = new float[n]; // Vec_DP g4(n)

	xsav=x;
	for (i=0;i<n;i++)
	{
		ysav[i]=y[i];
		dysav[i]=dydx[i];
	}
	
	jacobn_s(xsav,ysav,dfdx,dfdy, n); // n is nvar
	
	h=htry;
	for (jtry=0;jtry<MAXTRY;jtry++)
	{
		for (i=0;i<n;i++)
		{
			for (j=0;j<n;j++) a[i * n + j] = -dfdy[i * n + j];
			a[i * n + i] += 1.0/(GAM*h);
		}
		ludcmp(a,indx,d, n);
		for (i=0;i<n;i++)
			g1[i]=dysav[i]+h*C1X*dfdx[i];
		lubksb(a,indx,g1,n);
		for (i=0;i<n;i++)
			y[i]=ysav[i]+A21*g1[i];
		x=xsav+A2X*h;
		derivs(x,y,dydx);
		for (i=0;i<n;i++)
			g2[i]=dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h;
		lubksb(a,indx,g2,n);
		for (i=0;i<n;i++)
			y[i]=ysav[i]+A31*g1[i]+A32*g2[i];
		x=xsav+A3X*h;
		derivs(x,y,dydx);
		for (i=0;i<n;i++)
			g3[i]=dydx[i]+h*C3X*dfdx[i]+(C31*g1[i]+C32*g2[i])/h;
		lubksb(a,indx,g3,n);
		for (i=0;i<n;i++)
			g4[i]=dydx[i]+h*C4X*dfdx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;
		lubksb(a,indx,g4,n);
		for (i=0;i<n;i++) {
			y[i]=ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i];
			err[i]=E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i];
		}
		
		x=xsav+h;
		if (x == xsav) cout << "stepsize not significant in stiff" << endl;
		
		errmax=0.0f;
		for (i=0;i<n;i++) errmax=max(errmax,fabs(err[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0f)
		{
			hdid=h;
			hnext=(errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);
			delete[] a;
      delete[] dfdy;
      delete[] indx;
      delete[] dfdx;
      delete[] dysav;
      delete[] err;
      delete[] ysav;
      delete[] g1;
      delete[] g2;
      delete[] g3;
      delete[] g4;
			return;
		} else
		{
			hnext=SAFETY*h*pow(errmax,PSHRNK);
			h=(h >= 0.0f ? max(hnext,SHRNK*h) : min(hnext,SHRNK*h));
		}
	}
	cout << "exceeded MAXTRY in stiff" << endl;
}


//----- odeint_h
void odeint_h(float *ystart, const float x1, const float x2, const float eps,
	const float htry, const float hmin, int &nok, int &nbad, // htry is h1 in my other script if you are confused one day !!!!!!
	void derivs(const float, float*, float*))
{
	//const int MAXSTP = 10000;  // defined globally !
	//const float TINY = 1.0e-30; // defined globally !
	
	float xsav, x, hnext, hdid, h;

	// int nvar=ystart.size(); // Will be set as a global variable !
	float yscal[nvar], y[nvar], dydx[nvar];

	x = x1;
	h = sign_h(htry, x2 - x1);
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
			x_p[kount++] = x;
			xsav = x;
		}

		if ((x+h-x2)*(x+h-x1) > 0.0f) h=x2-x;

		//rkqs(y, dydx, x, h, eps, yscal, hdid, hnext, nvar, derivs);
		stiff(y, dydx, x, h, eps, yscal, hdid, hnext, nvar, derivs);

		if (hdid == h) ++nok; else ++nbad;

		if ((x - x2) * (x2 - x1) >= 0.0f)
		{
			for (i = 0; i < nvar; i++) ystart[i] = y[i];

			if (kmax != 0)
			{
				for (i = 0; i < nvar; i++) y_p[i * (kmax + 1) + kount] = y[i];
				x_p[kount++] = x;
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

  x_p = new float[kmax + 1]; // Declare globally - see above.
  y_p = new float[nvar * (kmax + 1)]; // Declare globally - see above.

	float y0[3] = {1.0f, 1.0f, 0.0f};
	// Note that here x1 and x2 represent time !
	const float x1 = 0.0f; // Initial time !!!!!!!!!!!!!!!!!!!!!!!
	const float x2 = 50.0f; // Final time !!!!!!!!!!!!!!!!!!!!!!!!!
	float htry = 2.9e-4; // Step size !!!!!!!!!!!!!!!!!!!!!!!!

	float eps = 1e-4f;   // THIS affects the execution time !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	float hmin = 1e-6f;  //!!!!!!!!!!!!!!!!!!!!!!!!!! To be adjusted for each problem !!!!!!!!!!!!!!!!!!!!!!!!!!

	int nok = 0;
	int nbad = 0;

	auto start = chrono::high_resolution_clock::now();

	odeint_h(y0, x1, x2, eps, htry, hmin, nok, nbad, derivs);

	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = end - start;
	cout << "Elapsed time: " << elapsed.count() << " seconds" << endl;
	
	cout << "kount = " << kount << endl;
	
	ofstream outfile("outX.csv"); // use pplot.py to plot the result !
  
	outfile << "x,y" << endl;

	for (int j = 0; j < kount; j++)
	{
		int n_row_1 = 1;
		int n_row_2 = 2;
		outfile << y_p[n_row_1 * (kmax+1) + j] << "," << y_p[n_row_2 * (kmax+1) + j] << endl;
	}

	outfile.close();
	
	cout << "y = " << y0[0] << ", " << y0[1] << ", " << y0[2] << endl; // The final evolved y values!
	
	delete[] x_p;
	delete[] y_p;
	
	
}



