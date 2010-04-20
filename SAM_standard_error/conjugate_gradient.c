#include <math.h>
#define NRANSI
#include "nrutil.h"
#include <stdio.h>
#define ITMAX 10000
#define EPS 1.0e-10
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);

void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []))
{
	void linmin(double p[], double xi[], int n, double *fret,
		double (*func)(double []));
	int j,its;
	double gg,gam,fp,dgg;
	double *g,*h,*xi;

	g=dvector(1,n);
	h=dvector(1,n);
	xi=dvector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		linmin(p,xi,n,fret,func);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp= *fret;
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"
#define TOL 2.0e-4

int ncom;
double *pcom,*xicom,(*nrfunc)(double []);

void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
	double brent(double ax, double bx, double cx,
		double (*f)(double), double tol, double *xmin);
	double f1dim(double x);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
		double *fc, double (*func)(double));
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;
	extern float LINMIN_TOL;

	ncom=n;
	pcom=dvector(1,n);
	xicom=dvector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,LINMIN_TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_dvector(xicom,1,n);
	free_dvector(pcom,1,n);
}
#undef TOL
#undef NRANSI
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define MAXRESCALE 1000 /* maximum number of rescalings of scale parameter */

#include <errno.h>

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
	double (*func)(double))
{
	double ulim,u,r,q,fu,dum;
	double utemp, scale; /* scale is a rescaling factor */
	int nrescale ; /* the number of times scale has been rescaled */

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	
	scale = 1.; 
	nrescale = 0; 
	while ( errno == EDOM ){
	  scale /= 3.;
	  errno = 0;
	  (*bx) *= scale;
	  *fb = (*func)(*bx);
	  nrescale++ ; 
	  if ( nrescale == MAXRESCALE ){
	    printf("Error: Model parameter values out of range (in mnbrak)\n"); exit(1);
	  }
	}  
	errno = 0;
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	scale = GOLD;
	nrescale = 0;
	while ( errno == EDOM ){ 
	  scale /= GOLD;
	  errno = 0;
	  *cx=(*bx)+scale*(*bx-*ax);
	  *fc=(*func)(*cx);
	  nrescale++ ; 
	  if ( nrescale == MAXRESCALE ){
	    printf("Error: Model parameter values out of range (in mnbrak)\n"); exit(1);
	  }
	}
	errno = 0;
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
			scale = GOLD;
			nrescale = 0;
			while ( errno == EDOM ){ 
			  scale /= GOLD;
			  errno = 0;
			  utemp=(*cx)+scale*(*cx-*bx);
			  fu=(*func)(utemp);
			  nrescale++ ; 
			  if ( nrescale == MAXRESCALE ){
			    printf("Error: Model parameter values out of range (in mnbrak)\n"); exit(1);
			  }
			}
			errno = 0; 

		} else if ((*cx-u)*(u-ulim) > 0.0) { 
		  fu=(*func)(u);
		  if ( errno == EDOM ){
		    u=(*cx)+GOLD*(*cx-*bx); 
		    errno=0;
		    fu=(*func)(u);
		    scale = GOLD;
		    nrescale = 0;
		    while ( errno == EDOM ){
		      scale /= GOLD;
		      errno = 0;
		      utemp=(*cx)+scale*(*cx-*bx);
		      fu=(*func)(utemp);
		      nrescale++;
		      if ( nrescale == MAXRESCALE ){
			printf("Error: Model parameter values out of range (in mnbrak)\n"); exit(1);
		      }
		    }
		    errno = 0; 
		  }
		  if (fu < *fc) {
		    utemp =  *cx+GOLD*(*cx-*bx); 
		    fu=(*func)(u);
		    scale = GOLD;
		    nrescale = 0;
		    while ( errno == EDOM ){ 
		      scale /= GOLD;
		      errno = 0;
		      utemp=(*cx)+scale*(*cx-*bx);
		      fu=(*func)(utemp);
		      nrescale++; 
		      if ( nrescale == MAXRESCALE ){
			printf("Error: Model parameter values out of range (in mnbrak)\n"); exit(1);
		      }
		    }
		    errno = 0; 
		    SHFT(*bx,*cx,u,utemp)
		    SHFT(*fb,*fc,fu,(*func)(utemp))
		  }

		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);

			scale = 1.;
			nrescale = 0; 
			while ( errno == EDOM ){ 
			  scale /= 3.;
			  errno = 0;
			  u =(*bx)+ scale * GLIMIT*(*cx-*bx);
			  if ( scale* GLIMIT <= 1 ) printf ("trouble: u is less than cx\n"); 
			  fu = (*func)(u);
			  nrescale++;
			  if ( nrescale == MAXRESCALE ){
			    printf("Error: Model parameter values out of range (in mnbrak)\n"); exit(1);
			  }
			}  
			errno = 0;

		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
			scale = GOLD;
			nrescale = 0;
			while ( errno == EDOM ){
			  scale /= GOLD;
			  errno = 0;
			  u=(*cx)+scale*(*cx-*bx);
			  fu=(*func)(u);
			  nrescale++;
			  if ( nrescale == MAXRESCALE ){
			    printf("Error: Model parameter values out of range (in mnbrak)\n"); exit(1);
			  }	
			}
			errno = 0; 

		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent(double ax, double bx, double cx, double (*f)(double), double tol,
	double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brent");
	*xmin=x;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI
/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"

extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double []);

double f1dim(double x)
{
	int j;
	double f,*xt;

	xt=dvector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_dvector(xt,1,ncom);
	return f;
}
#undef NRANSI
