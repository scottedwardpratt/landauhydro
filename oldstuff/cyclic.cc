#include <iostream>
#include <new>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>

using namespace std;
 
void solve(double* a, double* b, double* c, double* d, double *u,int n,double u0,double uminus1) 
{
	int i;
	u[0]=u0;
	for(i=1;i<=n;i++){
		u[i]=(d[i-1]-b[i-1]*u[i-1]-a[i-1]*uminus1)/c[i-1];
		uminus1=u[i-1];
	}
}

 
int main ()
{
	const double Pi=4.0*atan(1.0);
	const int n=100;
	double dx=0.01,dx2=pow(dx,2); 
	double u0a,uminus1a,u0b,uminus1b,u0d,uminus1d;
	double x[n+1],a[n+1],b[n+1],c[n+1],d[n+1],d0[n+1],AA[n+1],BB[n+1],rb[n+1];
	double ua[n+1],ub[n+1],ud[n+1],u[n+1],th_u[n+1]; 
     
	int i;
	for(i=0; i<n; i++)
	{
		x[i] = (i*dx);
		th_u[i] = sin(x[i]); 
	};
     
	for(i=0; i<=n; i++)
	{
		rb[i] = 7.0*sin(2.0*i*dx);   // driving force
		AA[i] = 0.5;                 // drag force
		BB[i] = 1.0;                 // restoring force
		
		a[i] = 1.0-0.5*dx*AA[i];
		b[i] = dx2*BB[i]-2.0;
		c[i] = 1.0+0.5*dx*AA[i];
		d[i] = dx2*rb[i];
		d0[i]=0.0;
	};
    
	// Find homogenous solutions
	u0a=0.0;  //arbitrary
	uminus1a=-dx;   //arbitrary
	solve(a,b,c,d0,ua,n,u0a,uminus1a);  
	
	u0b=1.0;   //arbitrary
	uminus1b=1.0-0.5*dx*dx;  //arbitrary
	solve(a,b,c,d0,ub,n,u0b,uminus1b);
	
	// Find particular solution
	u0d=0.0;  // arbitrary
	uminus1d=0.0;  // arbitrary
	solve(a,b,c,d,ud,n,u0d,uminus1d);
	
	double fa,fb; // coefficients for homogenous solutions
	double alpha1,beta1,alpha2,beta2,det,slope;
	alpha1=uminus1a-ua[n-1];
	beta1=uminus1b-ub[n-1];
	alpha2=u0a-ua[n];
	beta2=u0b-ub[n];
	det=alpha1*beta2-alpha2*beta1;
	fa=(beta2*ud[n-1]-beta1*ud[n])/det;
	fb=(alpha1*ud[n]-alpha2*ud[n-1])/det;
	printf("--- fa=%g, fb=%g ----\n",fa,fb);
	for(i=0;i<=n;i++){
		u[i]=ud[i]+fa*ua[i]+fb*ub[i];
		if(i<100){
			slope=u[i+1]-u[i];
			printf("u[%3d]=%9.6f %9.3f\n",i,u[i],slope);
		}
		else
			printf("u[%3d]=%9.6f\n",i,u[i]);
	}	
}
