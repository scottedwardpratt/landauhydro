#include "landau.h"
#include "eos.h"

using namespace std;

CLandau::CLandau(CparameterMap *parmapset){
	parmap=parmapset;
	string EOSDEF=parmap->getS("LANDAU_EOS","FreeGas");
	NDIM=parmap->getI("LANDAU_NDIM",3);
	DXYZ=parmap->getD("LANDAU_DXYZ",1);
	NX=parmap->getI("LANDAU_NX",100);
	NY=parmap->getI("LANDAU_NY",100);
	NZ=parmap->getI("LANDAU_NZ",100);
	//NT=parmap->getI("LANDAU_NT",1000);
	DELT=parmap->getD("LANDAU_DELT",0.01);
	printf("DELT=%g\n",DELT);
	TMAX=parmap->getD("LANDAU_TIMEMAX",1000.0);
	NT=lrint(TMAX/DELT);
	printf("NT=%d\n",NT);
	NRungeKutta=parmap->getI("LANDAU_NRUNGEKUTTA",2);
	output_dirname=parmap->getS("LANDAU_OUTPUT_DIRNAME","output");
	
	if(EOSDEF=="FreeGas"){
		eos=new CEoS_FreeGas(parmap);
	}
	else if(EOSDEF=="VdW")
		eos=new CEoS_VdW(parmap);
	else if(EOSDEF=="SCOTT1")
		eos=new CEoS_Scott(parmap);
	else{
		printf("LANDAU_EOS=%s, not recognized\n",EOSDEF.c_str());
		exit(1);
	}
		
	if(NDIM<3){
		NZ=1;
		if(NDIM<2)
			NY=1;
	}
	CLandauCell::DXYZ=DXYZ;
	CLandauCell::NDIM=NDIM;
	CLandauCell::eos=eos;
	CLandauMesh::eos=eos;
	CLandauMesh::NX=NX;
	CLandauMesh::NY=NY;
	CLandauMesh::NZ=NY;
	CLandauMesh::NDIM=NDIM;
	CLandauMesh::landau=this;
	CreateMeshes(parmap->getD("LANDAU_TIME0",0.0));
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
	oldmesh->Initialize(tset-DELT);
	currentmesh->Initialize(tset);
}

void CLandau::CycleMeshes(){
	CLandauMesh *crap;
	crap=oldmesh;
	oldmesh=currentmesh;
	currentmesh=newmesh;
	newmesh=crap;
}

void CLandau::Propagate(){
	PropagateRhoBPdens();
	newmesh->UpdateQuantities();
	newmesh->t=currentmesh->t+DELT;
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
				oldc->epsilonk=2.0*c->epsilonk-newc->epsilonk;
				if(oldc->epsilonk<0){
					printf("In InterpolateOldMesh, epsilonk=%g\n",oldc->epsilonk);
				}
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

void CLandau::PropagateRhoBPdens(){
	int ix,iy,iz,i,j,k;
	CLandauCell *c,*oldc,*newc;
	double DivKFlow,vdotgradki,vdotgradpiij;
	vector<double> DeliTij(4);
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				
				c=&(currentmesh->cell[ix][iy][iz]);
				oldc=&(oldmesh->cell[ix][iy][iz]);
				newc=&(newmesh->cell[ix][iy][iz]);
			//if(oldc->Pdens[0]<0.0){
				printf("BEFORE: old Pdens[0]=%g, new Pdens[0]=%g\n",oldc->Pdens[0],newc->Pdens[0]);
				//exit(1);
				//}
				
				newc->jB[0]=oldc->jB[0]-2.0*DELT*c->DelDotJB();
				if(newc->jB[0]<0.0 || newc->jB[0]!=newc->jB[0]){
					printf("negative or ill-defined density in CLandau::PropagateRhoBPdens(), rho=%g, t=%g\n",
					newc->jB[0],newmesh->t);
					printf("A: oldc->jB[0]=%g, c->jB[0]=%g\n",oldc->jB[0],c->jB[0]);
					newc->PrintInfo();
					exit(1);
				}
				c->CalcDeliTij(DeliTij);
				for(i=0;i<=NDIM;i++){
					newc->Pdens[i]=oldc->Pdens[i]-2.0*DELT*DeliTij[i];
					if(DeliTij[i]!=DeliTij[i]){
						printf("DeliTij=Nan\n");
						exit(1);
					}
				}
				if(newc->Pdens[1]!=newc->Pdens[1]){
					printf("ill-defined momentum density in CLandau::PropagateRhoBPdens(), %g, t=%g\n",
					newc->Pdens[1],newmesh->t);
					printf("B: oldc->Pdens[1]=%g, c->Pdens[1]=%g\n",oldc->Pdens[1],c->Pdens[1]);
					newc->PrintInfo();
					exit(1);
				}
				DivKFlow=c->CalcDivKFlow();
				newc->Pdens[0]-=2.0*DELT*DivKFlow;
				printf("DivKFlow=%g\n",DivKFlow);
				for(i=1;i<=NDIM;i++){
					newc->Kflow[i]=newc->alpha_K*(
						oldc->Kflow[i]/oldc->alpha_K
						-(2.0*DELT/c->tau_K)*(c->Kflow[i]/c->alpha_K-c->Kflow_target[i]/c->alpha_K)
						);

					vdotgradki=0.0;
					for(j=1;j<=NDIM;j++){
						
						vdotgradki+=c->u[j]*(c->neighborPlus[j]->Kflow[i]-c->neighborMinus[j]->Kflow[i])/(2.0*DXYZ);

						newc->pivisc[i][j]=newc->alpha_eta*(
							oldc->pivisc[i][j]/oldc->alpha_eta-(2.0*DELT/c->tau_eta)*(c->pivisc[i][j]/c->alpha_eta-c->pi_target[i][j]/c->alpha_eta));
						vdotgradpiij=0.0;
						for(k=1;k<=NDIM;k++){
							vdotgradpiij+=c->u[k]*(c->neighborPlus[k]->pivisc[i][j]-c->neighborMinus[k]->pivisc[i][j])/(2.0*DXYZ);
						}
						newc->pivisc[i][j]-=vdotgradpiij*DELT;
					}
					newc->Kflow[i]-=vdotgradki*DELT;
				}
				printf("tau_K=%g\n",c->tau_K);
				printf("oldc->alpha_K=%g, newc->alpha_K=%g\n",oldc->alpha_K,newc->alpha_K);
				printf("oldc->Kflow[1]=%g, c->Kflow[1]=%g, newc->Kflow[1]=%g, Kflow_target=%g\n",oldc->Kflow[1],c->Kflow[1],newc->Kflow_target[1]);
				if(newc->Pdens[0]<0.0){
					printf("after: Pdens[0]<0.0,=%g\n",newc->Pdens[0]);
					exit(1);
				}
			}
		}
	}
}

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

void CLandau::AverageMeshes_EvenQuantities(double weight){
	int ix,iy,iz,i,j;
	CLandauCell *c,*newc,*oldc;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(currentmesh->cell[ix][iy][iz]);
				oldc=&(oldmesh->cell[ix][iy][iz]);
				newc=&(newmesh->cell[ix][iy][iz]);
				i=0;
				c->Pdens[i]=(1.0-weight)*c->Pdens[i]+0.5*weight*(oldc->Pdens[i]+newc->Pdens[i]);
				c->jB[i]=(1.0-weight)*c->jB[i]+0.5*weight*(oldc->jB[i]+newc->jB[i]);
				for(i=1;i<=NDIM;i++){
					for(j=1;j<=NDIM;j++)
						c->pivisc[i][j]=(1.0-weight)*c->pivisc[i][j]+0.5*weight*(oldc->pivisc[i][j]+newc->pivisc[i][j]);
				}
			}
		}
	}	
	currentmesh->UpdateQuantities();	
}

void CLandau::AverageMeshes_OddQuantities(double weight){
	int ix,iy,iz,i;
	CLandauCell *c,*newc,*oldc;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(currentmesh->cell[ix][iy][iz]);
				oldc=&(oldmesh->cell[ix][iy][iz]);
				newc=&(newmesh->cell[ix][iy][iz]);
				for(i=1;i<=NDIM;i++){
					c->Pdens[i]=(1.0-weight)*c->Pdens[i]+0.5*weight*(oldc->Pdens[i]+newc->Pdens[i]);
					c->jB[i]=(1.0-weight)*c->jB[i]+0.5*weight*(oldc->jB[i]+newc->jB[i]);
				}
				for(i=1;i<=NDIM;i++){
					c->Kflow[i]=(1.0-weight)*c->Kflow[i]+0.5*weight*(oldc->Kflow[i]+newc->Kflow[i]);
				}
			}
		}
	}	
	currentmesh->UpdateQuantities();	
}
