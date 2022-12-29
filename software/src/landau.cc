#include "landau.h"
#include "eos.h"

using namespace std;

CLandau::CLandau(CparameterMap *parmapset){
	parmap=parmapset;
	CMeshParameters::DX0=parmap->getD("LANDAU_DX",1);
	CMeshParameters::NX=parmap->getI("LANDAU_NX",100);
	CMeshParameters::DT=parmap->getD("LANDAU_DT",1);
	CMeshParameters::TMAX=parmap->getD("LANDAU_TMAX",1);

	NRungeKutta=parmap->getI("LANDAU_NRUNGEKUTTA",2);
	output_dirname=parmap->getS("LANDAU_OUTPUT_DIRNAME","output");
	
	string EOSDEF=parmap->getS("LANDAU_EOS","FreeGas");
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
	CIntegralCell::eos=eos;
	
	CIntegralMesh::landau=CHalfIntegralMesh::landau=this;
	CreateMeshes(parmap->getD("LANDAU_TIME0",0.0));
}

void CLandau::Evolve(){
	double DT=CMeshParameters::DT,TMAX=CMeshParameters::TMAX;
	while(newIntegralMesh->t+0.1*DT<TMAX){
		CycleIntegralMeshes();
		newIntegralMesh->t+=DT;
		PropagateRhoSdensPI();
		CycleHalfIntegralMeshes();
		newHalfIntegralMesh->t+=DT;
		PropagateVxKx();
		newIntegralMesh->t=oldIntegralMesh->t+CMeshParameters::DT;
	}
}

void CLandau::CreateMeshes(double tset){
	double DT=CMeshParameters::DT;
	oldIntegralMesh=new CIntegralMesh(tset-DT);
	oldHalfIntegralMesh=new CHalfIntegralMesh(tset-0.5*DT);
	newHalfIntegralMesh=new CHalfIntegralMesh(tset);
	newIntegralMesh=new CIntegralMesh(tset+0.5*DT);
}

void CLandau::CycleIntegralMeshes(){
	CIntegralMesh *imesh_tmp;
	imesh_tmp=newIntegralMesh;
	newIntegralMesh=oldIntegralMesh;
	oldIntegralMesh=imesh_tmp;
}
void CLandau::CycleHalfIntegralMeshes(){
	CHalfIntegralMesh *himesh_tmp;
	himesh_tmp=newHalfIntegralMesh;
	newHalfIntegralMesh=oldHalfIntegralMesh;
	oldHalfIntegralMesh=himesh_tmp;
}

//-----------------------------------------------------------------

void CLandau::PropagateRhoSdensPI(){
	int ix,jx,NX=CMeshParameters::NX,DT=CMeshParameters::DT;
	for(ix=0;ix<NX;ix++){
		jx=ix+1;
		if(jx==NX)
			jx=0;
		newIntegralMesh->cell[ix]->Delx=oldIntegralMesh->cell[ix]->Delx
			+DT*(newHalfIntegralMesh->cell[jx]->vx-newHalfIntegralMesh->cell[ix]->vx);
		newIntegralMesh->cell[ix]->rho=oldIntegralMesh->cell[ix]->Q/newIntegralMesh->cell[ix]->Delx;
	}
	for(ix=0;ix<NX;ix++){
		newIntegralMesh->cell[ix]->CalcGrad2Rho();
	}
}

void CLandau::PropagateSdens(){
	int ix,jx,NX=CMeshParameters::NX,DT=CMeshParameters::DT;
	double oldS,newPi,oldPi,newDelx,oldDelx,newEta,oldEta;
	CIntegralCell *oldcell,*newcell;
	for(ix=0;ix<NX;ix++){
		jx=ix+1;
		if(jx==NX)
			jx=0;
		oldcell=oldIntegralMesh->cell[ix];
		newcell=newIntegralMesh->cell[ix];
		
		oldS=oldcell->S;
		newcell->S=oldS
			+DT*(newHalfIntegralMesh->cell[jx]->Kx-newHalfIntegralMesh->cell[ix]->Kx)/(0.5*(newcell->T+oldcell->T));
		
		newEta=newIntegralMesh->cell[ix]->eta;
		newDelx=newIntegralMesh->cell[ix]->Delx;
		newPi=newIntegralMesh->cell[ix]->Pi;
		oldEta=oldIntegralMesh->cell[ix]->eta;
		oldDelx=oldIntegralMesh->cell[ix]->Delx;
		oldPi=oldIntegralMesh->cell[ix]->Pi;
		newIntegralMesh->cell[ix]->S=oldIntegralMesh->cell[ix]->S
			+0.5*DT*((newDelx*newPi*newPi/newEta)+(oldDelx*oldPi*oldPi/oldEta));
	}
}

void CLandau::PropagatePi(){
	int ix,jx,NX=CMeshParameters::NX;
	//double DT=CMeshParameters::DT;
	//double oldPi,xnewEta,oldEta;
	//double tau_eta,alpha_eta;
	//CIntegralCell *oldcell,*newcell;
	for(ix=0;ix<NX;ix++){
		jx=ix+1;
		if(jx==NX)
			jx=0;
		//oldcell=oldIntegralMesh->cell[ix];
		//newcell=newIntegralMesh->cell[ix];
		//tau_eta=0.5*(oldcell->tau_eta+newcell->tau_eta);
		//alpha_eta=0.5*(oldcell->alpha_eta+newcell->alpha_eta);
		newHalfIntegralMesh->cell[jx]->GetOmega();		
	}
}

void CLandau::PropagateVxKx(){
	int ix,jx,NX=CMeshParameters::NX;
	CIntegralCell *bulkcella,*bulkcellb;
	CHalfIntegralCell *oldcell,*newcell;
	double gradP,gradT,tau_gmma,alpha_gmma,tmpterm;
	// First calc Vx
	for(ix=0;ix<NX;ix++){
		jx=ix-1;
		if(jx==-1)
			jx=NX-1;
		oldcell=oldHalfIntegralMesh->cell[ix];
		newcell=newHalfIntegralMesh->cell[ix];
		bulkcellb=newIntegralMesh->cell[jx];
		bulkcella=newIntegralMesh->cell[ix];
		
		gradP=bulkcellb->Pr-bulkcella->Pr;
		gradP+=0.5*CEos::kappa*(bulkcellb->grad2Rho-bulkcella->grad2Rho);
		gradP=gradP/(0.5*(bulkcella->Delx+bulkcellb->Delx));
		newcell->vx=oldcell->vx
			-gradP/(CEos::mass*0.5*(bulkcella->rho+bulkcellb->rho));
	}
	// Now calc Kx
	for(ix=0;ix<NX;ix++){
		jx=ix-1;
		if(jx==-1)
			jx=NX-1;
		oldcell=oldHalfIntegralMesh->cell[ix];
		newcell=newHalfIntegralMesh->cell[ix];
		bulkcellb=newIntegralMesh->cell[jx];
		bulkcella=newIntegralMesh->cell[ix];
		
		alpha_K_factor=0.5*(bulkcella->alpha_K_factor+bulkcellb->alpha_K_factor);
		gradT=2.0*(bulkcellb->T-bulkcella->T)/(bulkcella->Delx+bulkcellb->Delx);		
		
		
		tau_gmma=0.5*(bulkcella->tau_gmma+bulkcellb->tau_gmma);
		gmma=0.5*(bulkcella->gmma+bulkcellb->gmma);
		alpha_gmma=0.5*(bulkcella->alpha_gmma+bulkcellb->alpha_gmma);
		tmpterm=0.5*((1.0/tau_alpha)+alpha_gmma));
	
		Denom=(1.0/CMeshParameters::DT)+tmpterm;
		newcell->Kx=oldcell->Kx*((1.0/tau_alpha)-tmpterm)+gmma*gradT/tau_gmma;
		newcell->Kx=newcell->Kx/Denom;		
	}
	
}

/*
void CLandau::PrintInfo(){
	int ix,iy,iz;
	CIntegalCell *cnew,*cold;
	CHalfIntegralCell *ccurrent;
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

void CLandau::WriteData(){
	char filename[100];
	sprintf(filename,"output/data1d_t%g.dat",currentmesh->t);
	FILE *fptr=fopen(filename,"w");
	int ix,iy=0,iz=0;
	for(ix=0;ix<NX;ix++){
		fprintf(fptr,"%8.3e %10.4e\n",ix*DXYZ,currentmesh->cell[ix][iy][iz].jB[0]);
	}
	fclose(fptr);
}
*/

