
// The difference with __4 is that I removed ystart from odeint_h because it was needed only once in odeint_h and also stiff would always keep the original backed up !
// The difference with __3.cpp is that here I removed functions like derivs from the argument of another function like ode_int_h!

// The difference with stiff_hfv_cuda_2.cpp is that here we remove the x dependence, i.e. now it is dy/dx = f(y) but in _2.cpp it was dy/dx = f(x, y). We exclude
// the x dependence because in cooling derivative we do not have direct dependence on time !!

#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

using namespace std;

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
		else if (sum != 0.0f)
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
	const float TINY=1.0e-20f;
	int i,imax,j,k;
	float big,dum,sum,temp;

	//int n=a.nrows();
	float vv[n]; // Vec_DP vv(n);
	d=1.0f;
	for (i=0;i<n;i++) {
		big=0.0f;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i * n + j])) > big) big=temp;
		if (big == 0.0f) cout << "Singular matrix in routine ludcmp" << endl;
		vv[i]=1.0f/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i * n + j];
			for (k=0;k<i;k++) sum -= a[i * n + k]*a[k * n + j];
			a[i * n + j]=sum;
		}
		big=0.0f;
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
		if (a[j * n + j] == 0.0f) a[j * n + j]=TINY;
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
void jacobn_s(float *y, float *dfdx, float *dfdy, int n) // Note: dfdy is a Matrix !
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
void derivs(float *y, float *dydx)
{
	dydx[0] = -0.013f*y[0]-1000.0f*y[0]*y[2];
	dydx[1] = -2500.0f*y[1]*y[2];
	dydx[2] = -0.013f*y[0]-1000.0f*y[0]*y[2]-2500.0f*y[1]*y[2];
}


//----- stiff
void stiff(float *y, float *dydx, float &x, const float htry,
	const float eps, float *yscal, float &hdid, float &hnext, int n)
{
	const float SAFETY=0.9f,GROW=1.5f,PGROW= -0.25f,SHRNK=0.5f;
	const float PSHRNK=(-1.0f/3.0f),ERRCON=0.1296f;
	const int MAXTRY=40;
	const float GAM=1.0f/2.0f,A21=2.0f,A31=48.0f/25.0f,A32=6.0f/25.0f,C21= -8.0f,
		C31=372.0f/25.0f,C32=12.0f/5.0f,C41=(-112.0f/125.0f),
		C42=(-54.0f/125.0f),C43=(-2.0f/5.0f),B1=19.0f/9.0f,B2=1.0f/2.0f,
		B3=25.0f/108.0f,B4=125.0f/108.0f,E1=17.0f/54.0f,E2=7.0f/36.0f,E3=0.0f,
		E4=125.0f/108.0f,C1X=1.0f/2.0f,C2X=(-3.0f/2.0f),C3X=(121.0f/50.0f),
		C4X=(29.0f/250.0f),A2X=1.0f,A3X=3.0f/5.0f;
	int i,j,jtry;
	float d,errmax,h,xsav;

	float *a = new float[n * n];
	float *dfdy = new float[n * n];
	
	int *indx = new int[n];
	
	float *dfdx = new float[n];
	float *dysav = new float[n];
	float *err = new float[n];
	float *ysav = new float[n];
	float *g1 = new float[n];
	float *g2 = new float[n];
	float *g3 = new float[n];
	float *g4 = new float[n];

	xsav=x;
	for (i=0;i<n;i++)
	{
		ysav[i]=y[i]; // ysav will always contain the original y, so it is a kind of back up!
		dysav[i]=dydx[i]; // dysav is also like ysav -- a back up !!
	}
	
	jacobn_s(ysav,dfdx,dfdy, n); // n is nvar
	
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
		derivs(y,dydx);
		for (i=0;i<n;i++)
			g2[i]=dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h;
		lubksb(a,indx,g2,n);
		for (i=0;i<n;i++)
			y[i]=ysav[i]+A31*g1[i]+A32*g2[i];
		x=xsav+A3X*h;
		derivs(y,dydx);
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
		if (x == xsav) cout << "stepsize not significant in stiff function!" << endl;
		
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
	cout << "exceeded MAXTRY in stiff function !" << endl;
}


//----- odeint_h
void odeint_h(float *y, const float x1, const float x2, const float eps, // y was ystart!
	const float htry, const float hmin, int nvar) // htry is h1 in my other script if you are confused one day !!!!!!
{
	const int MAXSTP = 10000;
	const float TINY = 1.0e-30f;
	
	float xsav, x, hnext, hdid, h;

	float yscal[nvar], dydx[nvar]; // , y[nvar]

	x = x1;
	h = sign_h(htry, x2 - x1);

	int i, nstp;

	//for (i= 0; i < nvar; i++) y[i] = ystart[i];

	for (nstp = 0;nstp < MAXSTP; nstp++) 
	{
		derivs(y, dydx);
		for (i = 0;i < nvar; i++)
			yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;

		if ((x+h-x2)*(x+h-x1) > 0.0f) h=x2-x;

		stiff(y, dydx, x, h, eps, yscal, hdid, hnext, nvar); // y will be backed up in stiff so that its original value is always available in stiff !

		if ((x - x2) * (x2 - x1) >= 0.0f)
		{
			//for (i = 0; i < nvar; i++) ystart[i] = y[i];

			return;
		}
		if (fabs(hnext) <= hmin) cout <<"Step size too small in odeint" <<endl;
		h=hnext;
	}
	cout <<"Too many steps in routine odeint" << endl;
}



int main()
{

  const int nvar = 3; // see derivs_s !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	float y0[3] = {1.0f, 1.0f, 0.0f}; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	cout << "y0 = " << y0[0] << ", " << y0[1] << ", " << y0[2] << endl;
	
	// Note that here x1 and x2 represent time !
	const float x1 = 0.0f; // Initial time !!!!!!!!!!!!!!!!!!!!!!!
	const float x2 = 50.0f; // Final time !!!!!!!!!!!!!!!!!!!!!!!!!
	float htry = 2.9e-4f; // Step size !!!!!!!!!!!!!!!!!!!!!!!!

	float eps = 1e-4f;   // THIS affects the execution time !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	float hmin = 1e-6f;  //!!!!!!!!!!!!!!!!!!!!!!!!!! To be adjusted for each problem !!!!!!!!!!!!!!!!!!!!!!!!!!

	auto start = chrono::high_resolution_clock::now();

	odeint_h(y0, x1, x2, eps, htry, hmin, nvar); // the final evolved y will be overwritten in y0 !

	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = end - start;
	cout << "Elapsed time: " << elapsed.count() << " seconds" << endl;
	
	cout << "y = " << y0[0] << ", " << y0[1] << ", " << y0[2] << endl; // final evolved y !

	
}



