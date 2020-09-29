
     #include <iostream>
     #include <new>
     #include <cstdlib>
     #include <cmath>
     #include <cstdio>
     #include <complex>

     using namespace std;

     void solve(float* a, float* b, float* c, float* d, int n) 
     {
     n--; 
     c[0] /=b[0];
     d[0] /=b[0];

     for (int i = 1; i < n; i++) 
     {
     c[i] /= b[i] - a[i]*c[i-1];
     d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
     }

     d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

     for (int i = n; i-- > 0;) 
     {
     d[i] -= c[i]*d[i+1];
     }
     }
//------------------------------------------------------------------

     int main ()
     { 
     const double gama=3.0, fac=0.02;
     int i,nx=100;
     float dx=0.01, dx2=pow(dx,2), dx3=pow(dx,3), x0=0.0, x1=1.0, 
           dg, dg2, xmid, dt=0.003, tend=0.65, ak=1.0; 
     float x[102],xi[103],u[102],q1[102],q2[102],q2p[102],q3[102],f1[101],f2[101],f3[101],
           Pr[101],q1new[102],q2new[102],q2pnew[102],q3new[102],
           term1[100],term2[100],term3[100],term4[100],
           termt1[100],termt2[100],termt3[100],termt4[100],
           a[100],b[100],c[100],d[100],AA[100],BB[100],rb[100]; 

//-------------------------------------------------------------      
//-----------SETTING  THE  GRID--------------------------------
     for(i=0; i<=(nx+1); i++)
     {
     x[i]=(i-1.)*dx;
     };
     for(i=1; i<=(nx+1); i++)
     {
     xi[i]=0.5*(x[i]+x[i-1]);
     };
     xi[0]=x[0]-0.5*(x[1]- x[0]);
     xi[nx+2]=x[nx+1]+0.5*(x[nx+1]-x[nx]);
//-------------------------------------------------------------
//----------ASSIGNING INITIAL CONDITION------------------------
 
     dg=0.1*(x1-x0);
     dg2=pow(dg,2);
     xmid= 0.5*(x0+x1);

     for(i=0; i<nx; i++)
     {
     q1[i]=exp(-(pow((x[i]-xmid),2))/(2.*pow(dg,2)));
     q2[i]=0.5*q1[i];
     q3[i]=0.5*q1[i]*fac+0.5*pow(q2[i],2)/q1[i];
     };
//--------------------------------------------------------------
//--------------ADVECTION---------------------------------------
 
     for(float t=0.0; t<=tend; t=t+dt)
     {
//--------------------------------------------------------------
//----------PERIODIC BOUNDART CONDITION2-------------------------    
     q1[-1]  =q1[nx-1]; 
     q1[nx]  =q1[0]; 
     q1[nx+1]=q1[1];
     q1[-2]  =q1[nx-2];
     
     q2[-1]  =q2[nx-1];
     q2[nx]  =q2[0];
     q2[nx+1]=q2[1];
     q2[-2]  =q2[nx-2];
     
     q3[-1]  =q3[nx-1];
     q3[nx]  =q3[0];
     q3[nx+1]=q3[1];
     q3[-2]  =q3[nx-2];
//--------------------------------------------------------------
//---------------SETTING  VELOCITY & PRESSURE-------------------
      
     for(i=-1; i<=nx; i++)
     {
//   u[i]=fac*t*(x[i]-xmid)/(fac*pow(t,2)+dg2);     
     u[i]=(q2[i]/q1[i]+q2[i-1]/q1[i-1])/2.;
     Pr[i]=(gama-1.)*(q3[i]-0.5*pow(q2[i],2)/q1[i]);
     };
//--------------------------------------------------------------         
     
     for(i=0; i<=nx; i++)
     {
       if (u[i] >= 0)
       {
       f1[i]=q1[i-1]*u[i]+(u[i]/4.)*(1.-u[i]*dt/dx)*(q1[i]-q1[i-2]);
       f2[i]=q2[i-1]*u[i]+(u[i]/4.)*(1.-u[i]*dt/dx)*(q2[i]-q2[i-2]);
       f3[i]=q3[i-1]*u[i]+(u[i]/4.)*(1.-u[i]*dt/dx)*(q3[i]-q3[i-2]);
       } 
       else
       {
       f1[i]=q1[i]*u[i]-(u[i]/4.)*(1.+u[i]*dt/dx)*(q1[i+1]-q1[i-1]);
       f2[i]=q2[i]*u[i]-(u[i]/4.)*(1.+u[i]*dt/dx)*(q2[i+1]-q2[i-1]);
       f3[i]=q3[i]*u[i]-(u[i]/4.)*(1.+u[i]*dt/dx)*(q3[i+1]-q3[i-1]);
       };
     };
   
 
    for(i=0; i<nx; i++)
    {
     
     term1[i]=u[i]*q1[i]*(q1[i+2]-2.*q1[i+1]+2.*q1[i-1]-q1[i-2])/(2.*dx3);
     term2[i]=1.5*q1[i]*(u[i+1]-u[i-1])*(q1[i+1]+q1[i-1]-2.*q1[i])/(2.*dx3);
     term3[i]=1.5*q1[i]*(q1[i+1]-q1[i-1])*(u[i+1]+u[i-1]-2.*u[i])/(2.*dx3);
     term4[i]=0.5*pow(q1[i],2)*(u[i+2]-2.*u[i+1]+2.*u[i-1]-u[i-2])/(2.*dx3);

     termt1[i]=0.5*u[i]*pow((q1[i+1]-q1[i-1]),2)/(4.*dx2);
     termt2[i]=u[i]*q1[i]*(q1[i+1]+q1[i-1]-2.*q1[i])/dx2;
     termt3[i]=0.5*q1[i]*(q1[i+1]-q1[i-1])*(u[i+1]-u[i-1])/(4.*dx2);
     termt4[i]=0.5*pow(q1[i],2)*(u[i+1]+u[i-1]-2.*u[i])/dx2;

     q2p[i]=q2[i]+ak*(termt1[i]-termt2[i]-termt3[i]-termt4[i]);

     q1new[i]=q1[i]-dt*(f1[i+1]-f1[i])/(xi[i+1]-xi[i]);

     q2pnew[i]=q2p[i]-dt*(f2[i+1]-f2[i])/(xi[i+1]-xi[i])-(dt/(2.*dx))*(Pr[i+1]-Pr[i-1])
                   +ak*dt*q1[i]*(q1[i+2]-2.*q1[i+1]+2.*q1[i-1]-q1[i-2])/(2.*dx3);

     q3new[i]=q3[i]-dt*(f3[i+1]-f3[i])/(xi[i+1]-xi[i])-(dt/(2.*dx))*
                                    (Pr[i+1]*u[i+1]-Pr[i-1]*u[i-1])
                   +ak*dt*(term1[i]+term2[i]+term3[i]+term4[i]);
         
      q1[i]=q1new[i];
      q2p[i]=q2pnew[i];
      q3[i]=q3new[i];
        
    };

          for(i=0; i<nx; i++)
          {
          rb[i]=-q2p[i]/((ak/2.)*pow(q1[i],2));
          AA[i]=(q1[i+1]-q1[i-1])/(2.*dx*q1[i]);
          BB[i]=-pow((q1[i+1]-q1[i-1]),2)/(4.*dx2*pow(q1[i],2))
                +2.*(q1[i+1]+q1[i-1]-2.*q1[i])/(dx2*q1[i])
                -2./(ak*q1[i]);

          a[i]=1.-(dx/2.)*AA[i];
          b[i]=dx2*BB[i]-2.;
          c[i]=1.+(dx/2.)*AA[i];
          };

          for(i=1; i<nx-1; i++)
          {
          d[i]=dx2*rb[i];
          };

          d[0]=dx2*rb[0]-a[0]*u[0];
          d[nx-1]=dx2*rb[nx-1]-c[nx-1]*u[nx-1];

          solve(a,b,c,d,nx);

          for(i=0; i<nx; i++)
          {
          u[i]=d[i];
          q2[i]=q1[i]*u[i];
          };
     
    }; 

    for(i=0; i<nx; i++)
    {
    printf("%g %g %g\n",x[i],q1[i],q3[i]);
    }; 

    }

//-----------------------------------------------------------------------------------------

