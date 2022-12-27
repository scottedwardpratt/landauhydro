#include "landau.h"
#include "eos.h"
CLandau *CLandauMesh::landau=NULL;
CEoS *CLandauMesh::eos=NULL;
int CLandauMesh::NX=1;
int CLandauMesh::NY=1;
int CLandauMesh::NZ=1;
int CLandauMesh::NDIM=3;

void CIntegralMesh(){
	int ix; 
	int NX=CMeshParameters::NX;
	double DX=CMeshParameters::DX;

	cell.resize(NX);

	for(ix=0;ix<NX;ix++){
		if(ix>0)
			cell[ix].neighborMinus=&cell[ix-1];
		else
			cell[ix].neighborMinus=&cell[NX-1];
		if(ix<NX-1)
			cell[ix].neighborPlus=&cell[ix+1];
		else
			cell[ix].neighborPlus=&cell[0];
		cell[ix].Zero();		
	}
}

void CIntegralMesh::Initialize(double tset){
	int ix,i,j;
	int nx=1;
	CIntegralCell *c;
	double x,kx,rho0=0.2,Drho,T0,grad2rhoB,Dtemp;
	T0=landau->parmap->getD("LANDAU_INIT_TEMP",0.15);
	rho0=landau->parmap->getD("LANDAU_INIT_RHO",0.25);
	Drho=landau->parmap->getD("LANDAU_INIT_DRHO",0.1);
	Dtemp=landau->parmap->getD("LANDAU_INIT_DTEMP",0.1);
	t=tset;
	Lx=NX*DX;
	kx=2.0*PI*nx/Lx;
	for(ix=0;ix<NX;ix++){
		x=DX*(ix+0.5);
		c=&cell[ix];
		c->Zero();
		c->rho=rho00*(1.0+Drho*cos(kx*x));
		grad2rhoB=-kx*kx*Drho*rho0*cos(kx*x);
		c->T=T0*(1.0-Dtemp*cos(kx*x));
		eos->CalcEoS_of_rho_T(c);
		//printf("ix=%3d, Pdens[0]=%g, epsilonk=%g, rhoB=%g, grad2rhoB=%g\n",ix,c->Pdens[0],c->epsilonk,c->jB[0],grad2rhoB);
	}
}

void CIntegralMesh::UpdateQuantities(){
	int ix,iy,iz,i;
	CIntegralCell *c;
	for(ix=0;ix<NX;ix++){
		cell[ix]->UpdateBulkQuantities();
	}
	for(ix=0;ix<NX;ix++){
		c=&cell[ix];
		c->Calc_Kflow_target();
		c->Calc_pi_target();
	}
}

void CIntegralMesh::Zero(){
	
}

void CIntegralMesh::CalculateBtotStot(){
	double Btot=0.0,Stot=0.0;
	int ix;
	CIntegralCell *c;
	for(ix=0;ix<NX;ix++){
		c=&(cell[ix]);
		Btot+=c->rho*cell->Delx;
		Stot+=c->SoverB*c->rho*cell->Delx;
	}
	//printf("_________________ Calculating Btot, Etot, Stot for t=%g.  _______________\n",t);
	printf("## Btot=%g, Stot=%g\n",Btot,Stot);
}

///////////////////////

CHalfIntegralMesh::CHalfIntegralMesh(){
	int ix; 
	int NX=CMeshParameters::NX;
	double DX=CMeshParameters::DX;

	cell.resize(NX);

	for(ix=0;ix<NX;ix++){
		if(ix>0)
			cell[ix].neighborMinus=&cell[ix-1];
		else
			cell[ix].neighborMinus=&cell[NX-1];
		if(ix<NX-1)
			cell[ix].neighborPlus=&cell[ix+1];
		else
			cell[ix].neighborPlus=&cell[0];
		cell[ix].Zero();		
	}
}

void CHalfIntegralMesh::Initialize(double tset){
	CHalfIntegralCell *c;
	double Dvel,Lx,kx;
	int ix,nx=1;
	t=tset;
	Lx=NX*DX;
	Dvel=landau->parmap->getD("LANDAU_INIT_DVEL",0.0);
	kx=2.0*PI*nx/Lx;
	for(ix=0;ix<NX;ix++){
		x=DX*(ix+0.5);
		c=&cell[ix];
		c->vx=Dvel*sin(kx*x);
		c->Kx=0.0;
	}		
}