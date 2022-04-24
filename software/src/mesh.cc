#include "landau.h"
#include "eos.h"
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
	int ix,iy,iz,i,j;
	int nx=1,ny=0,nz=0;
	CLandauCell *c;
	double x,y,z,Lx,Ly,Lz,kx,ky,kz,jB0=0.2,Drho,T0,grad2rhoB,Dtemp;
	T0=landau->parmap->getD("LANDAU_INIT_TEMP",0.15);
	jB0=landau->parmap->getD("LANDAU_INIT_RHO",0.25);
	Drho=landau->parmap->getD("LANDAU_INIT_DRHO",0.1);
	Dtemp=landau->parmap->getD("LANDAU_INIT_DTEMP",0.1);
	t=tset;
	Lx=NX*DXYZ;
	Ly=NY*DXYZ;
	Lz=NZ*DXYZ;
	kx=2.0*PI*nx/Lx;
	ky=2.0*PI*ny/Ly;
	kz=2.0*PI*nz/Lz;
	for(ix=0;ix<NX;ix++){
		x=DXYZ*(ix+0.5);
		for(iy=0;iy<NY;iy++){
			y=DXYZ*(iy+0.5);
			for(iz=0;iz<NZ;iz++){
				z=DXYZ*(iz+0.5);
				c=&cell[ix][iy][iz];
				c->Zero();
				c->jB[0]=jB0*(1.0+Drho*cos(kx*x)*cos(ky*y)*cos(kz*z));
				grad2rhoB=-(kx*kx+ky*ky+kz*kz)*Drho*jB0*cos(kx*x)*cos(ky*y)*cos(kz*z);
				c->T=T0*(1.0-Dtemp*cos(kx*x)*cos(ky*y)*cos(kz*z));
				eos->CalcEoS_of_rho_T(c);
				
				c->Pdens[0]=c->epsilonk-0.5*eos->kappa*c->jB[0]*grad2rhoB;
				c->Pdens[1]=c->Pdens[2]=c->Pdens[3]=0.0;
				//printf("ix=%3d, Pdens[0]=%g, epsilonk=%g, rhoB=%g, grad2rhoB=%g\n",ix,c->Pdens[0],c->epsilonk,c->jB[0],grad2rhoB);
			}
		}
	}
	UpdateQuantities();
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&cell[ix][iy][iz];
				for(i=1;i<=NDIM;i++){
					c->kflow[i]=c->kflow_target[i];
					for(j=1;j<=NDIM;j++){
						c->pivisc[i][j]=c->pitarget[i][j];
					}
				}
			}
		}
	}
	//CalculateBtotEtot();
	printf("Initialization finished for t=%g\n",t);
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
	double maxdens=0.0,mindens=1000000.0,x;
	char filename[140]; 
	sprintf(filename,"%s/xslice_t%06.1f.dat",landau->output_dirname.c_str(),t);
	FILE *fptr=fopen(filename,"w");
	CLandauCell *c;
	//	printf("----------TIME=%g -------------\n",currentmesh->t);
	for(ix=0;ix<NX;ix++){
		c=&(cell[ix][iy][iz]);
		x=(ix+0.5)*DXYZ;
		fprintf(fptr,"%8.4f %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5f\n",
		x,c->jB[0],c->jB[1],c->epsilonk,c->u[1],c->Pr,c->SE[1][1],c->T,c->kflow[1],c->cs2);
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
	//printf("_________________ Calculating Btot, Etot, Stot for t=%g.  _______________\n",t);
	printf("$$$$$$ Btot=%g, Etot=%g, Stot=%g $$$$$$$",Btot,Etot,Stot);
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

void CLandauMesh::UpdateQuantities(){
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
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&cell[ix][iy][iz];
				c->Calckflow_target();
				c->Calcpi_target();
			}
		}
	}
}

