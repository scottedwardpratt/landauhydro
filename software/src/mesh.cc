#include "landau.h"
CLandau *CLandauMesh::landau=NULL;
CEoS *CLandauMesh::eos=NULL;

CLandauMesh::CLandauMesh(CLandau *landau,double tset){
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

void CLandauMesh::InitializeDensities(){
	int ix,iy,iz,i,j,NDIM=landau->NDIM;
	double T,Pr;
	CLandauCell *c;
	double T0=0.01;
	double x,y,z,Lx,Ly,Lz;
	Lx=NX*DXYZ;
	Ly=NY*DXYZ;
	Lz=NZ*DXYZ;

	for(ix=0;ix<NX;ix++){
		x=DXYZ*ix;
		for(iy=0;iy<NY;iy++){
			y=DXYZ*iy;
			for(iz=0;iz<NZ;iz++){
				z=DXYZ*iz;
				c=&cell[ix][iy][iz];
				c->Zero();
				c->jB[0]=1.0+0.01*cos(2.0*PI*x/Lx)*cos(2.0*PI*y/Ly)*cos(2.0*PI*z/Lz);
				c->epsilon=eos->mass*c->jB[0]+c->jB[0]*1.5*T0;
				//printf("--------------\nBefore: mass=%g, epsilon=%g, rhoB=%g, T0=%g\n",
				//eos->mass,c->epsilon,c->jB[0],T0);	
				eos->eos(c->epsilon,c->jB[0],T,Pr);
				//printf("after: T=%g =? %g, Pr=%g\n",T,T0,Pr);
				for(i=0;i<=NDIM;i++){
					c->u[i]=0.0;
					if(i>0)
						c->jB[i]=0.0;
					for(j=0;j<=NDIM;j++){
						c->SE[i][j]=0.0;
						if(i==j)
							c->SE[i][j]=Pr;
					}
				}
			}
		}
	}
}

void CLandauMesh::WriteInfo(){
	int ix,iy,iz,i;
	double x,y,z;
	char filename[140];
	sprintf(filename,"rhoB_%g.dat",t);
	FILE *fptr=fopen(filename,"w");
	CLandauCell *c;
	//	printf("----------TIME=%g -------------\n",currentmesh->t);
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				x=ix*DXYZ;
				y=iy*DXYZ;
				z=iy*DXYZ;
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

void CLandauMesh::PrintInfo(){
	int ix,iy,iz,i;
	double x,y,z;
	CLandauCell *c;
	printf("----------TIME=%g -------------\n",t);
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				x=ix*DXYZ;
				y=iy*DXYZ;
				z=iy*DXYZ;
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


