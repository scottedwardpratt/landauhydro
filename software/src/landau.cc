#include "landau.h"

using namespace std;

CLandau::CLandau(CparameterMap *parmapset){
	parmap=parmapset;
	NDIM=parmap->getI("LANDAU_NDIM",3);
	DXYZ=parmap->getD("LANDAU_DXYZ",1);
	NX=parmap->getI("LANDAU_NX",100);
	NY=parmap->getI("LANDAU_NY",100);
	NZ=parmap->getI("LANDAU_NZ",100);
	DELT=parmap->getD("LANDAU_DELT",0.01);
	NRungeKutta=parmap->getI("LANDAU_NRUNGEKUTTA",2);
	eos=new CEoS(parmap);
	if(NDIM<3){
		NZ=1;
		if(NDIM<2)
			NY=1;
	}
	CLandauCell::DXYZ=DXYZ;
	CLandauCell::NDIM=NDIM;
	CLandauCell::eos=eos;
	CLandauMesh::eos=eos;
	CLandauMesh::landau=this;
	CreateMeshes(parmap->getD("LANDAU_T0",0.0));
}

//-------------------------------------------------------------
void CLandau::solve(double* a, double* b, double *c, double *d, int n){
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
//-------------------------------------------------------------
void CLandau::CreateMeshes(double tset){
	oldmesh=new CLandauMesh(this,tset-DELT);
	currentmesh=new CLandauMesh(this,tset);
	newmesh=new CLandauMesh(this,tset+DELT);
	currentmesh->InitializeDensities();
	oldmesh->InitializeDensities();
}

void CLandau::CycleMeshes(){
	CLandauMesh *crap;
	crap=oldmesh;
	oldmesh=currentmesh;
	currentmesh=newmesh;
	newmesh=crap;
}

void CLandau::Propagate(){
	PropagateRhoB();
	PropagateT0i();
	PropagateT00();
	CalcEpsilonU();
	newmesh->t=currentmesh->t+DELT;
}

void CLandau::PropagateFirst(){
	DELT=0.5*DELT;
	PropagateRhoB();
	PropagateT0i();
	PropagateT00();
	CalcEpsilonU();
	InterpolateOldMesh();
	DELT=2.0*DELT;
	PropagateRhoB();
	PropagateT0i();
	PropagateT00();
	CalcEpsilonU();
}

void CLandau::InterpolateOldMesh(){
	int ix,iy,iz,j,k;
	CLandauCell *oldc,*newc,*c;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(currentmesh->cell[ix][iy][iz]);
				oldc=&(oldmesh->cell[ix][iy][iz]);
				newc=&(newmesh->cell[ix][iy][iz]);
				oldc->epsilon=2.0*c->epsilon-newc->epsilon;
				oldc->Pr=2.0*c->Pr-newc->Pr;
				oldc->T=2.0*c->T-newc->T;
				for(k=0;k<=NDIM;k++){
					oldc->jB[k]=2.0*c->jB[k]-newc->jB[k];
					oldc->u[k]=2.0*c->u[k]-newc->u[k];
					oldc->M[k]=2.0*c->M[k]-newc->M[k];
					for(j=0;j<=NDIM;j++){
						oldc->SE[j][k]=2.0*c->SE[j][k]-newc->SE[j][k];
					}
				}
			}
		} 	
	}
}

//-----------------------------------------------------------------
void CLandau::PropagateRhoB(){ 
	int ix,iy,iz;
	CLandauCell *c,*oldc,*newc;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(currentmesh->cell[ix][iy][iz]);
				oldc=&(oldmesh->cell[ix][iy][iz]);
				newc=&(newmesh->cell[ix][iy][iz]);
				newc->jB[0]=oldc->jB[0]-2.0*DELT*c->DelDotJB();
				//printf("check: ix=%d, DeltDotJB=%g\n",ix,c->DelDotJB());
			}
		}
	} 

}

//----------------------------------------------------------------
void CLandau::PropagateT0i(){
	int ix,iy,iz,k;
	vector<double> DeliTij(NDIM+1);
	CLandauCell *c,*oldc,*newc;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(currentmesh->cell[ix][iy][iz]);
				oldc=&(oldmesh->cell[ix][iy][iz]);
				newc=&(newmesh->cell[ix][iy][iz]);
				c->CalcDeliTij(DeliTij);
				for(k=1;k<=NDIM;k++){
					newc->SE[0][k]=oldc->SE[0][k]-2.0*DELT*DeliTij[k];
					newc->SE[k][0]=newc->SE[0][k];
				}
				printf("%d %g\n",ix,DeliTij[1]);
			}
		}
	}
}

//---------------------------------------------------------------
void CLandau::PropagateT00(){
	int ix,iy,iz;
	CLandauCell *c,*oldc,*newc;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(currentmesh->cell[ix][iy][iz]);
				oldc=&(oldmesh->cell[ix][iy][iz]);
				newc=&(newmesh->cell[ix][iy][iz]);
				newc->SE[0][0]=oldc->SE[0][0]-2.0*DELT*c->DelDotT0i(); 
			}
		}
	} 
}

//----------------------------------------------------------------
void CLandau::PrintInfo(){
	int ix,iy,iz;
	CLandauCell *cnew,*cold,*ccurrent;
	printf("----------TIME=%g -------------\n",newmesh->t);
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				ccurrent=&(currentmesh->cell[ix][iy][iz]);
				cold=&(oldmesh->cell[ix][iy][iz]);
				cnew=&(newmesh->cell[ix][iy][iz]);
				printf("%3d %3d %3d ",ix,iy,iz);
				
				printf(" %12.5e ",cold->jB[0]);
				printf(" %12.5e ",ccurrent->jB[0]);
				printf(" %12.5e ",cnew->jB[0]);
				
				printf(" %12.5e ",cold->jB[1]);
				printf(" %12.5e ",ccurrent->jB[1]);
				printf(" %12.5e ",cnew->jB[1]);
				
				if(NY>1){
					printf(" %12.5e ",cold->jB[2]);
					printf(" %12.5e ",ccurrent->jB[2]);
					printf(" %12.5e ",cnew->jB[2]);
				}
				if(NZ>1){
					printf(" %12.5e ",cold->jB[3]);
					printf(" %12.5e ",ccurrent->jB[3]);
					printf(" %12.5e ",cnew->jB[3]);
				}
				printf("\n");
			}
		} 
	}
}

void CLandau::WriteInfo(){
	int ix,iy,iz;
	double time=newmesh->t;
	double maxdens=0.0;
	char filename[100];
	sprintf(filename,"output/current_t%g.dat",time);
	CLandauCell *cnew,*cold,*ccurrent;
	FILE *fptr=fopen(filename,"w");
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				ccurrent=&(currentmesh->cell[ix][iy][iz]);
				cold=&(oldmesh->cell[ix][iy][iz]);
				cnew=&(newmesh->cell[ix][iy][iz]);
				fprintf(fptr,"%3d %3d %3d ",ix,iy,iz);
				
				fprintf(fptr," %12.5e ",cold->jB[0]);
				fprintf(fptr," %12.5e ",ccurrent->jB[0]);
				fprintf(fptr," %12.5e ",cnew->jB[0]);
				if(cnew->jB[0]>maxdens)
					maxdens=cnew->jB[0];
				
				fprintf(fptr," %12.5e ",cold->jB[1]);
				fprintf(fptr," %12.5e ",ccurrent->jB[1]);
				fprintf(fptr," %12.5e ",cnew->jB[1]);
				
				if(NY>1){
					fprintf(fptr," %12.5e ",cold->jB[2]);
					fprintf(fptr," %12.5e ",ccurrent->jB[2]);
					fprintf(fptr," %12.5e ",cnew->jB[2]);
				}
				if(NZ>1){
					fprintf(fptr," %12.5e ",cold->jB[3]);
					fprintf(fptr," %12.5e ",ccurrent->jB[3]);
					fprintf(fptr," %12.5e ",cnew->jB[3]);
				}
				fprintf(fptr,"\n");
			}
		} 
	}
	printf("max density=%g\n",maxdens);
}

void CLandau::CalcEpsilonU(){
	int ix,iy,iz,irk;
	CLandauCell *c,*oldc,*newc;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				newc=&(newmesh->cell[ix][iy][iz]);
				c=&(currentmesh->cell[ix][iy][iz]);
				oldc=&(oldmesh->cell[ix][iy][iz]);
				newc->CalcM();
				newc->Pr=2.0*c->Pr-oldc->Pr;
				printf("Pr[%d]=%g\n",ix,newc->Pr);
			}
		}
	} 
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				newc=&(newmesh->cell[ix][iy][iz]);
				for(irk=0;irk<NRungeKutta;irk++){
					newc->CalcEpsilonU();
				}
			}
		} 
		
	}
}

void CLandau::WriteData1D(){
	char filename[100];
	sprintf(filename,"output/data1d_t%g.dat",currentmesh->t);
	FILE *fptr=fopen(filename,"w");
	int ix,iy=0,iz=0;
	for(ix=0;ix<NX;ix++){
		fprintf(fptr,"%8.3e %10.4e\n",ix*DXYZ,currentmesh->cell[ix][iy][iz].jB[0]);
	}
	fclose(fptr);
}
