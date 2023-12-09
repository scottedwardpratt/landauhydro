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
	CIntegralCell::eos=CIntegralMesh::eos=eos;
	
	CIntegralMesh::landau=CHalfIntegralMesh::landau=this;
	CreateMeshes(parmap->getD("LANDAU_TIME0",0.0));
}

void CLandau::Evolve(){
	int irk,nrk=parmap->getD("LANDAU_NRUNGEKUTTA",2);
	double DT=CMeshParameters::DT,TMAX=CMeshParameters::TMAX;
	while(newIntegralMesh->t+0.1*DT<TMAX){

		CycleIntegralMeshes();
		newIntegralMesh->t+=DT;	
		// SDens is altered by entropy production, which requires knowing Pi and temperature at half-integral times
		for(irk=0;irk<nrk;irk++){
			PropagateRho();
			if(irk==0){
				EstimatePiS();
				CalcBulkQuantities();
			}
			PropagateSdens();
			CalcBulkQuantities();
			CalcPiTarget();
			CalcKxTarget();
			PropagatePi();
		}

		CycleHalfIntegralMeshes();
		newHalfIntegralMesh->t+=DT;	
		PropagateVxKx();

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

void CLandau::EstimatePiS(){
	int ix,NX=CMeshParameters::NX,jx,kx;
	CIntegralCell *oldcell,*newcell;
	for(ix=0;ix<NX;ix++){
		oldcell=oldIntegralMesh->cell[ix];
		newcell=newIntegralMesh->cell[ix];
		newcell->sdens=2.0*oldcell->sdens-newcell->sdens;
		for(jx=1;jx<4;jx++){
			newcell->pi_shear[jx][jx]=2.0*oldcell->pi_shear[jx][jx]-newcell->pi_shear[jx][jx];
			for(kx=jx+1;kx<4;kx++){
				newcell->pi_shear[ix][jx]=2.0*oldcell->pi_shear[ix][jx]-newcell->pi_shear[ix][jx];
				newcell->pi_shear[jx][ix]=newcell->pi_shear[ix][jx];
			}
		}
	}
}

void CLandau::PropagateRho(){
	int ix,jx,NX=CMeshParameters::NX,DT=CMeshParameters::DT;
	for(ix=0;ix<NX;ix++){
		jx=ix+1;
		if(jx==NX)
			jx=0;
		newIntegralMesh->cell[ix]->Delx=oldIntegralMesh->cell[ix]->Delx
			+DT*(newHalfIntegralMesh->cell[jx]->Vx-newHalfIntegralMesh->cell[ix]->Vx);
		newIntegralMesh->cell[ix]->rho=oldIntegralMesh->cell[ix]->Q/newIntegralMesh->cell[ix]->Delx;
	}
	for(ix=0;ix<NX;ix++){
		newIntegralMesh->cell[ix]->CalcGrad2Rho();
	}
}

void CLandau::PropagateSdens(){
	int ix,jx,j,k,NX=CMeshParameters::NX,DT=CMeshParameters::DT;
	double eta,T,delx,pijk;
	CIntegralCell *oldcell,*newcell;
	CHalfIntegralCell *halfcell,*halfcellplus;
	for(ix=0;ix<NX;ix++){
		jx=ix+1;
		if(jx==NX)
			jx=0;
		oldcell=oldIntegralMesh->cell[ix];
		newcell=newIntegralMesh->cell[ix];
		halfcell=newHalfIntegralMesh->cell[ix];
		halfcellplus=halfcell->neighborPlus;
		
		eta=0.5*(newcell->eta+oldcell->eta);
		T=0.5*(newcell->T+oldcell->T);
		delx=0.5*(newcell->Delx+oldcell->Delx);
	
		newcell->S=oldcell->S+DT*(halfcellplus->Kx-halfcell->Kx)/T;
		for(j=1;j<4;j++){
			for(k=1;k<4;k++){
				pijk=0.5*(newcell->pi_shear[j][k]+oldcell->pi_shear[j][k]);
				newcell->S+=DT*delx*pijk*pijk/eta;
			}
		}
	}
}

void CLandau::PropagatePi(){
	int ix,NX=CMeshParameters::NX,j,k;
	double DT=CMeshParameters::DT;
	double eta,tau_eta,alpha_eta,deldotv,delx,pipi,pitarget;
	CIntegralCell *oldcell,*newcell;
	CHalfIntegralCell *halfcell;
	
	for(ix=0;ix<NX;ix++){
		oldcell=oldIntegralMesh->cell[ix];
		newcell=newIntegralMesh->cell[ix];
		halfcell=newHalfIntegralMesh->cell[ix];
		
		eta=0.5*(newcell->eta+oldcell->eta);
		tau_eta=0.5*(oldcell->tau_eta+newcell->tau_eta);
		alpha_eta=0.5*(oldcell->alpha_eta+newcell->alpha_eta);
		delx=0.5*(newcell->Delx+oldcell->Delx);
		deldotv=0.0;
		for(j=1;j<4;j++){
			for(k=1;k<4;k++){
				halfcell->omega[j][k]=0.0;
			}
		}
		
		deldotv=(halfcell->neighborPlus->Vx-halfcell->Vx)/delx;
		halfcell->omega[1][1]=(4.0/3.0)*deldotv;
		halfcell->omega[2][2]=halfcell->omega[3][3]=-0.5*halfcell->omega[1][2];
		
		for(j=1;j<4;j++){
			for(k=1;k<4;k++){
				pitarget=-eta*halfcell->omega[j][k];
				pipi=0.5*(oldcell->pi_shear[j][k]+newcell->pi_shear[j][k]);
				newcell->pi_shear[j][k]-=(DT/(tau_eta*alpha_eta))*(pipi-pitarget);
			}
		}
		
		newHalfIntegralMesh->cell[ix]->GetOmega();
	}
}

void CLandau::PropagateVxKx(){
	int ix,jx,NX=CMeshParameters::NX;
	CIntegralCell *bulkcella,*bulkcellb;
	CHalfIntegralCell *oldcell,*newcell;
	double gradP,gradT,tau_gmma,alpha_gmma,tmpterm,alpha_Kx_factor,gmma,Denom;
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
		gradP+=0.5*eos->kappa*(bulkcellb->grad2Rho-bulkcella->grad2Rho);
		gradP=gradP/(0.5*(bulkcella->Delx+bulkcellb->Delx));
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
		
		alpha_Kx_factor=0.5*(bulkcella->alpha_Kx_factor+bulkcellb->alpha_Kx_factor);
		gradT=2.0*(bulkcellb->T-bulkcella->T)/(bulkcella->Delx+bulkcellb->Delx);		
		
		
		tau_gmma=0.5*(bulkcella->tau_gmma+bulkcellb->tau_gmma);
		gmma=0.5*(bulkcella->gmma+bulkcellb->gmma);
		alpha_gmma=alpha_Kx_factor*0.5*(bulkcella->alpha_gmma+bulkcellb->alpha_gmma);
		tmpterm=0.5*((1.0/tau_gmma)+alpha_gmma);
	
		Denom=(1.0/CMeshParameters::DT)+tmpterm;
		tau_gmma=eos->Kfactor/alpha_gmma;
		newcell->Kx=oldcell->Kx*((1.0/tau_gmma)-tmpterm)+gmma*gradT/tau_gmma;
		newcell->Kx=newcell->Kx/Denom;
	}
	
}

void CLandau::CalcPiTarget(){
	
}

