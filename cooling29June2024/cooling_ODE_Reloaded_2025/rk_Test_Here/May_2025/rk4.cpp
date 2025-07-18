#include "nr.h"

void NR::rk4(Vec_I_DP &y, Vec_I_DP &dydx, const DP x, const DP h,
	Vec_O_DP &yout, void derivs(const DP, Vec_I_DP &, Vec_O_DP &))
{
	int i;
	DP xh,hh,h6;

	int n=y.size();
	Vec_DP dym(n),dyt(n),yt(n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i]; // dydx = k1
	derivs(xh,yt,dyt); // dyt = k2
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i]; // dyt = k2
	derivs(xh,yt,dym); // dym = k3
	for (i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i]; // dym = k3
		dym[i] += dyt[i]; // k2 + k3. Note that dyt is k2 from above and we are adding it to dym which is k3. Then dym contains both k2 and k3 !
	}
	derivs(x+h,yt,dyt); // Now dyt is overwritten as k4 because it was already added to dym as k2 so it is not needed any more and we can asigned it to be k4!
	for (i=0;i<n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]); // in paranthesis we have, in order (k1 + k4 + 2.0*(k2+k3))
}
