#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <iostream>
#include "../software/include/randy.h"

const int D=100;  // dimension
const double DXYZ=1.0;
using namespace std;

int main(int argc, char** argv){
	CRandy randy(-1234);
	int ix,iy,iz;
	double x,y,z;
	double m2=1.0,k2=12.0*PI*PI/(D*D*DXYZ*DXYZ);
	double denom,E2=(k2+m2);
	printf("E2=%g\n",E2);
	double ***phi,***phi0,***phidummy,***b;
	phi=new double **[D];
	phi0=new double **[D];
	b=new double **[D];
	for(ix=0;ix<D;ix++){
		x=DXYZ*ix;
		phi[ix]=new double *[D];
		phi0[ix]=new double *[D];
		b[ix]=new double *[D];
		for(iy=0;iy<D;iy++){
			y=DXYZ*iy;
			phi[ix][iy]=new double[D];
			phi0[ix][iy]=new double[D];
			b[ix][iy]=new double [D];
			for(iz=0;iz<D;iz++){
				z=DXYZ*iz;
				b[ix][iy][iz]=E2*sin(2.0*PI*double(x)/double(D*DXYZ))
									*sin(2.0*PI*double(y)/double(D*DXYZ))
										*sin(2.0*PI*double(z)/double(D*DXYZ));
				phi0[ix][iy][iz]=b[ix][iy][iz]/E2;
				phi0[ix][iy][iz]+=(0.00002/E2)*randy.ran();
				phi[ix][iy][iz]=phi0[ix][iy][iz];
			}
		}
	}
	printf("arrays created\n");
	//
	int iter,niters;
	printf("How many iterations?\n");
	scanf("%d",&niters);
	int ixm,iym,izm,ixp,iyp,izp;
	double A;
	for(iter=0;iter<niters;iter++){
		printf("iter=%d\n",iter);
		for(ix=0;ix<D;ix++){
			ixm=ix-1;
			if(ixm==-1)
				ixm=D-1;
			ixp=ix+1;
			if(ixp==D)
				ixp=0;
			for(iy=0;iy<D;iy++){
				iym=iy-1;
				if(iym==-1)
					iym=D-1;
				iyp=iy+1;
				if(iyp==D)
					iyp=0;
				for(iz=0;iz<D;iz++){
					izm=iz-1;
					if(izm==-1)
						izm=D-1;
					izp=iz+1;
					if(izp==D)
						izp=0;
					denom=m2+6.0/(DXYZ*DXYZ);
					phi[ix][iy][iz]=b[ix][iy][iz]/denom;
					phi[ix][iy][iz]+=(phi0[ixm][iy][iz]+phi0[ixp][iy][iz]
						+phi0[ix][iym][iz]+phi0[ix][iyp][iz]+phi0[ix][iy][izm]+phi0[ix][iy][izp])/(denom*DXYZ*DXYZ);
				}
			}
		}
		phidummy=phi;
		phi=phi0;
		phi0=phidummy;
	}
	
	double ratio;
	for(ix=0;ix<D;ix++){
		for(iy=0;iy<D;iy++){
			for(iz=0;iz<D;iz++){
				ratio=E2*phi0[ix][iy][iz]/b[ix][iy][iz];
				printf("%3d %3d %3d %10.6f =? %10.6f  1.0=?%g\n",ix,iy,iz,phi0[ix][iy][iz]*E2,b[ix][iy][iz],ratio);
			}
		}
	}
	printf("E2=%g\n",E2);
	
	return 0;
}
