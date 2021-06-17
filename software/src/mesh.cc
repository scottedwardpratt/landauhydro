#include "landau.h"
CLandau *CLandauMesh::landau=NULL;
CEoS *CLandauMesh::eos=NULL;
double CLandauMesh::DXYZ=0.0;
int CLandauMesh::NX=1;
int CLandauMesh::NY=1;
int CLandauMesh::NZ=1;
int CLandauMesh::NDIM=3;

CLandauMesh::CLandauMesh(CLandau *landauset,double tset){
	landau=landauset;
	t=tset;
	int ix,iy,iz; 
	NX=landau->NX; NY=landau->NY; NZ=landau->NZ;
	DXYZ=landau->DXYZ;
	cell.resize(NX);
	for(ix=0;ix<NX;ix++){
		cell[ix].resize(NY);
		for(iy=0;iy<NY;iy++){
			cell[ix][iy].resize(NZ);
		}
	}
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				
				if(ix>0)
					cell[ix][iy][iz].neighborMinus[1]=&cell[ix-1][iy][iz];
				else
					cell[ix][iy][iz].neighborMinus[1]=&cell[NX-1][iy][iz];
				
				if(ix<NX-1)
					cell[ix][iy][iz].neighborPlus[1]=&cell[ix+1][iy][iz];
				else
					cell[ix][iy][iz].neighborPlus[1]=&cell[0][iy][iz];
				 
				if(iy>0)
					cell[ix][iy][iz].neighborMinus[2]=&cell[ix][iy-1][iz];
				else
					cell[ix][iy][iz].neighborMinus[2]=&cell[ix][NY-1][iz];
				
				if(iy<NY-1) 
					cell[ix][iy][iz].neighborPlus[2]=&cell[ix][iy+1][iz];
				else
					cell[ix][iy][iz].neighborPlus[2]=&cell[ix][0][iz];
				
				if(iz>0)
					cell[ix][iy][iz].neighborMinus[3]=&cell[ix][iy][iz-1];
				else
					cell[ix][iy][iz].neighborMinus[3]=&cell[ix][iy][NZ-1];
				
				if(iz<NZ-1)
					cell[ix][iy][iz].neighborPlus[3]=&cell[ix][iy][iz+1];
				else 
					cell[ix][iy][iz].neighborPlus[3]=&cell[ix][iy][0];
								
				cell[ix][iy][iz].Zero();
			}
		}
	}
}

void CLandauMesh::Initialize(double tset){
	int ix,iy,iz;
	int nx=1,ny=0,nz=0;
	CLandauCell *c;
	double x,y,z,Lx,Ly,Lz,kx,ky,kz,jB0=0.5,T0=2.0/27.0,cs2,cs20,omega0,SoverB0,epsilon0=1.5,P0,Arho=0.01;
	t=tset;
	Lx=NX*DXYZ;
	Ly=NY*DXYZ;
	Lz=NZ*DXYZ;
	kx=2.0*PI*nx/Lx;
	ky=2.0*PI*ny/Ly;
	kz=2.0*PI*nz/Lz;
	eos->eos(epsilon0,jB0,T0,P0,SoverB0,cs20);
	omega0=sqrt((kx*kx+ky*ky+kz*kz)*cs20);
	
	for(ix=0;ix<NX;ix++){
		x=DXYZ*ix;
		for(iy=0;iy<NY;iy++){
			y=DXYZ*iy;
			for(iz=0;iz<NZ;iz++){
				z=DXYZ*iz;
				c=&cell[ix][iy][iz];
				c->Zero();
				c->jB[0]=jB0+Arho*cos(kx*x)*cos(ky*y)*cos(kz*z)*cos(omega0*t);
				c->T=T0*pow(c->jB[0]/jB0,2.0/3.0);
				c->Pdens[0]=1.5*c->jB[0]*c->T;
				eos->eos(c->Pdens[0],c->jB[0],c->T,c->Pr,c->SoverB,cs2);
				c->Pdens[1]=((eos->mass*cs20*kx*Arho)/(omega0))*sin(kx*x)*cos(ky*y)*cos(kz*z)*sin(omega0*t);
				c->Pdens[2]=((eos->mass*cs20*ky*Arho)/(omega0))*cos(kx*x)*sin(ky*y)*cos(kz*z)*sin(omega0*t);
				c->Pdens[3]=((eos->mass*cs20*kz*Arho)/(omega0))*cos(kx*x)*cos(ky*y)*sin(kz*z)*sin(omega0*t);
			}
		}
	}
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&cell[ix][iy][iz];
				CalculateUJMEpsilonSE();
			}
		}
	}
}

void CLandauMesh::WriteInfo(){
	int ix,iy,iz;
	char filename[140];
	sprintf(filename,"rhoB_%g.dat",t);
	FILE *fptr=fopen(filename,"w");
	CLandauCell *c;
	//	printf("----------TIME=%g -------------\n",currentmesh->t);
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(cell[ix][iy][iz]);
				fprintf(fptr,"%3d %3d %3d %12.5e",ix,iy,iz,c->jB[0]);
				fprintf(fptr," %12.5e",c->jB[1]);
				if(NY>1)
					fprintf(fptr," %12.5e",c->jB[2]);
				if(NZ>1)
					fprintf(fptr," %12.5e",c->jB[3]);
				fprintf(fptr,"\n");
			}
		} 
		fclose(fptr);
	}
}

void CLandauMesh::WriteXSliceInfo(int iy,int iz){
	int ix;
	double maxdens=0.0,mindens=1000000.0;
	char filename[140]; 
	sprintf(filename,"output/xslice_t%g.dat",t);
	FILE *fptr=fopen(filename,"w");
	CLandauCell *c;
	//	printf("----------TIME=%g -------------\n",currentmesh->t);
	for(ix=0;ix<NX;ix++){
		c=&(cell[ix][iy][iz]);
		fprintf(fptr,"%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",ix*DXYZ,c->jB[0],c->jB[1],c->epsilon,c->u[1],c->Pr,c->T);
		if(c->jB[0]>maxdens)
			maxdens=c->jB[0];
		if(c->jB[0]<mindens)
			mindens=c->jB[0];
	} 
	printf("t=%g: %g < dens < %g\n",t,mindens,maxdens);
	fclose(fptr);
}

void CLandauMesh::CalculateBtotEtot(){
	double Btot=0.0,Stot=0.0,Etot=0.0;
	vector<double> Ptot(NDIM+1);
	double dOmega=pow(DXYZ,NDIM);
	int ix,iy,iz,i;
	CLandauCell *c;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(cell[ix][iy][iz]);
				Btot+=c->jB[0]*dOmega;
				Stot+=c->SoverB*c->jB[0]*dOmega;
				Etot+=c->Pdens[0]*dOmega;
				for(i=0;i<=NDIM;i++)
					Ptot[i]+=c->Pdens[i]*dOmega;
			}
		}
	}
	printf("t=%g, Btot=%g, Etot=%g, Stot=%g\n",t,Btot,Etot,Stot);
	for(i=1;i<=NDIM;i++)
		printf(", Ptot[%d]=%g",i,Ptot[i]);
	printf("\n");
}

void CLandauMesh::PrintInfo(){
	int ix,iy,iz;
	CLandauCell *c;
	printf("----------TIME=%g -------------\n",t);
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(cell[ix][iy][iz]);
				printf("%3d %3d %3d %12.5e",ix,iy,iz,c->jB[0]);
				printf(" %12.5e",c->jB[1]);
				if(NY>1)
					printf(" %12.5e",c->jB[2]);
				if(NZ>1)
					printf(" %12.5e",c->jB[3]);
				printf("\n");
			}
		} 
	}
}

void CLandauMesh::CalculateUJMEpsilonSE(){
	int ix,iy,iz,i;
	CLandauCell *c;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				for(i=1;i<=NDIM;i++){
					c=&cell[ix][iy][iz];
					c->u[i]=c->Pdens[i]/(eos->mass*c->jB[0]);
					c->jB[i]=c->u[i]*c->jB[0];
				}
			}
		}
	}
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&cell[ix][iy][iz];
				c->CalcM();
				c->CalcEpsilonSE();
			}
		}
	}
}

