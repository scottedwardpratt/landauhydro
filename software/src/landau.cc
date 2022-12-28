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
	double t,DT=CMeshParameters::DT,TMAX=CMeshParameters::TMAX;
	while(newIntegralMesh->t+0.1*DT<TMAX){
		CycleIntegralMeshes();
		newIntegralMesh->t+=DT;
		PropagateRhoSdensPI();
		CycleHalfIntegralMeshes();
		newHalfIntegralMesh->t+=DT;
		PropagateVxKx();
		newmesh->t=currentmesh->t+DELT;
	}
}

void CLandau::CreateMeshes(double tset){
	oldIntegralMesh=new CLandauMesh(tset-DELT);
	oldHalfIntegralMesh=new CLandauMesh(tset-0.5*DELT);
	newHalfIntegralMesh=new CLandauMesh(tset);
	newIntegralMesh=new CLandauMesh(tset+0.5*DELT);
}

void CLandau::CycleIntegralMeshes(){
	CIntegralMesh *imesh_tmp;
	imesh_temp=newIntegralMesh;
	newIntegralMesh=oldIntegralMesh;
	oldIntegralMesh=imesh_tmp;
}
void CLandau::CycleHalfIntegralMeshes(){
	CHalfIntegralMesh *himesh_tmp;
	himesh_temp=newHalfIntegralMesh;
	newHalfIntegralMesh=oldHalfIntegralMesh;
	oldHalIntegralMesh=himesh_tmp;
}

//-----------------------------------------------------------------

void CLandau::PropagateRho(){
	int ix,jx,NX=CMeshParameters::NX,DT=CMeshParameters::DT;
	for(ix=0;ix<NX;ix++){
		jx=ix+1;
		if(jx==NX)
			jx=0;
		newIntegralMesh->cell[ix].DelX=oldIntegralMesh->cell[ix].DelX
			+DT*(newHalfIntegralMesh->cell[jx].Vx-newHalIntegralMesh->cell[ix].Vx);
		newIntegralMesh->cell[ix].rho=oldIntegralMesh->cell[ix].Q/newIntegralMesh->cell[ix].Delx;
	}
	for(ix=0;ix<NX;ix++){
		newIntegralMesh->cell[ix].CalcGrad2Rho();
	}
}

void CLandau::PropagateSDens(){
	int ix,jx,NX=CMeshParameters::NX,DT=CMeshParameters::DT;
	double oldS,newPi,oldPi,newDelX,oldDelX,newEta,oldEta;
	for(ix=0;ix<NX;ix++){
		jx=ix+1;
		if(jx==NX)
			jx=0;
		oldcell=oldIntegralMesh->cell[ix];
		newcell=newIntegralMesh->cell[ix];
		
		oldS=oldcell->S;
		newcell->S=oldS
			+DT*(newHalfIntegralMesh->cell[jx].Kx-newHalIntegralMesh->cell[ix].Kx)/(0.5*(newcell->T+oldcell->T));
		
		newEta=newIntegralMesh->cell[ix].eta;
		newDelX=newIntegralMesh->cell[ix].DelX;
		newPi=newIntegralMesh->cell[ix].Pi;
		oldEta=oldIntegralMesh->cell[ix].eta;
		oldDelX=oldIntegralMesh->cell[ix].DelX;
		oldPi=oldIntegralMesh->cell[ix].Pi
		newIntegralMesh->cell[ix].S=oldIntegralMesh->cell[ix].S
			+0.5*DT*((newDelX*newPi*newPi/newEta)+(oldDelX*oldPi*oldPi/oldEta));
	}
}

void CLandau::PropagatePi(){
	int ix,jx,NX=CMeshParameters::NX,DT=CMeshParameters::DT;
	double oldPi,newDelX,oldDelX,newEta,oldEta,tmpterm,omega;
	for(ix=0;ix<NX;ix++){
		jx=ix+1;
		if(jx==NX)
			jx=0;
		oldcell=oldIntegralMesh->cell[ix];
		newcell=newIntegralMesh->cell[ix];
		tau_eta=0.5*(oldcell->tau_eta+newcell->tau_eta);
		alpha_eta=0.5*(oldcell->alpha_eta+newcell->alpha_eta);
		halfIntegralCell->GetOmega(omega);
		
		
		
	}
}

void CLandau::PropagateRhoSdensPI(){
	int ix,NX=CMeshParameters::NX;
	double DX=CMeshParameters::DX;
	CIntegralCell *newcell,oldcell;
	PropagateRho();
	for(ix=0;ix<NX;ix++){
		oldcell=oldIntegralMesh->cell[ix];
		newcell=newIntegralMesh->cell[ix];
		
		newcell->T=oldcell->T;
		newcell->zeta=oldcell->zeta;
		newcell->Pi=oldcell->Pi;
		newcell->gamma=oldcell->gamma;
		//newcell->eta=oldcell->eta;
	}
	for(int irk=0;irk<NRungeKuta;irk++){
		PropagateSdens();
		newIntegralMesh->UpdateQuantities();
		PropagatePi();
		newIntegralMesh->UpdateQuantities();
	}
}

void CLandau::PropagateVxKx(){
	int ix,jx,NX=CMeshParameters::NX;
	CIntegralCell *bulkcella,*bulkcellb;
	double gradP,gradT,tau_gamma,alpha_gamma,tmpterm;
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
		gradP+=0.5*eos->kappa*(bulkcellb->grad2rho-bulkcella->grad2rho);
		gradP=gradP/(0.5*(bulkcella->DelX+bulkcellb->DelX));
		newcell->Vx=oldcell->Vx
			-gradP/(eos->mass*0.5*(bulkcella->rho+bulkcellb->rho));
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
		
		oldKx=oldcell
		
		alpha_K_factor=0.5*(bulkcella->alpha_K_factor+bulkcellb->alpha_K_factor);
		gradT=2.0*(bulkcellb->T-bulkcella->T)/(bulkcella->DelX+bulkcellb->Delx);		
		
		
		tau_gamma=0.5*(bulkcella->tau_gamma+bulkcellb->tau_gamma);
		gamma=0.5*(bulkcella->gamma+bulkcellb->gamma);
		alpha_gamma=0.5*(bulkcella->alpha_gamma+bulkcellb->alpha_gamma);
		tmpterm=0.5*((1.0/tau_alpha)+alpha_gamma));
	
		Denom=(1.0/DT)+tmpterm;
		newcell->Kx=oldcell->Kx*((1.0/tau_alpha)-tmpterm)+gamma*gradT/tau_gamma;
		newcell->Kx=newcell->Kx/Denom;		
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
