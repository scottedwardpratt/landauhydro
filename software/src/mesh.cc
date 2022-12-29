#include "landau.h"
#include "eos.h"
CLandau *CIntegralMesh::landau=NULL;
CLandau *CHalfIntegralMesh::landau=NULL;
CEoS *CIntegralMesh::eos=NULL;
int CMeshParameters::NX=1;
double CMeshParameters::DX0=0.0;
double CMeshParameters::DT=0.0;
double CMeshParameters::TMAX=0.0;

CIntegralMesh::CIntegralMesh(double tset){
	t=tset;
	int ix; 
	int NX=CMeshParameters::NX;
	cell.resize(NX);
	for(ix=0;ix<NX;ix++){
		if(ix>0)
			cell[ix]->neighborMinus=cell[ix-1];
		else
			cell[ix]->neighborMinus=cell[NX-1];
		if(ix<NX-1)
			cell[ix]->neighborPlus=cell[ix+1];
		else
			cell[ix]->neighborPlus=cell[0];
		cell[ix]->Zero();		
	}
}

void CIntegralMesh::Initialize(double tset){
	int ix;
	int nx=1;
	CIntegralCell *c;
	double x,Lx,kx,rho0=0.2,Drho,T0,Dtemp;
	T0=landau->parmap->getD("LANDAU_INIT_TEMP",0.15);
	rho0=landau->parmap->getD("LANDAU_INIT_RHO",0.25);
	Drho=landau->parmap->getD("LANDAU_INIT_DRHO",0.1);
	Dtemp=landau->parmap->getD("LANDAU_INIT_DTEMP",0.1);
	t=tset;
	Lx=CMeshParameters::NX*CMeshParameters::DX0;
	kx=2.0*PI*nx/Lx;
	for(ix=0;ix<CMeshParameters::NX;ix++){
		x=CMeshParameters::DX0*(ix+0.5);
		c=cell[ix];
		c->Zero();
		c->rho=rho0*(1.0+Drho*cos(kx*x));
		//grad2rho=-kx*kx*Drho*rho0*cos(kx*x);
		c->T=T0*(1.0-Dtemp*cos(kx*x));
		eos->CalcEoS_of_rho_T(c);
	}
}

void CIntegralMesh::UpdateQuantities(){
	int ix;
	CIntegralCell *c;
	for(ix=0;ix<CMeshParameters::NX;ix++){
		cell[ix]->UpdateQuantities();
	}
	for(ix=0;ix<CMeshParameters::NX;ix++){
		c=cell[ix];
		c->Calc_Kx_target();
	}
}

void CIntegralMesh::CalculateBtotStot(){
	double Btot=0.0,Stot=0.0;
	int ix;
	CIntegralCell *c;
	for(ix=0;ix<CMeshParameters::NX;ix++){
		c=cell[ix];
		Btot+=c->rho*cell[ix]->Delx;
		Stot+=c->SoverB*c->rho*cell[ix]->Delx;
	}
	//printf("_________________ Calculating Btot, Etot, Stot for t=%g.  _______________\n",t);
	printf("## Btot=%g, Stot=%g\n",Btot,Stot);
}

///////////////////////

CHalfIntegralMesh::CHalfIntegralMesh(double tset){
	t=tset;
	int ix; 
	int NX=CMeshParameters::NX;
	cell.resize(NX);
	for(ix=0;ix<NX;ix++){
		if(ix<NX-1)
			cell[ix]->neighborPlus=cell[ix+1];
		else
			cell[ix]->neighborPlus=cell[0];
		cell[ix]->Zero();	
	}
}

void CHalfIntegralMesh::Initialize(double tset){
	CHalfIntegralCell *c;
	double Dvel,Lx,kx,x;
	int ix,nx=1;
	t=tset;
	Lx=CMeshParameters::NX*CMeshParameters::DX0;
	Dvel=landau->parmap->getD("LANDAU_INIT_DVEL",0.0);
	kx=2.0*PI*nx/Lx;
	for(ix=0;ix<CMeshParameters::NX;ix++){
		c=cell[ix];
		x=ix*CMeshParameters::DX0;
		c->vx=Dvel*sin(kx*x);
		c->Kx=0.0;
	}		
}

void CHalfIntegralMesh::UpdateQuantities(){
	int ix;
	for(ix=0;ix<CMeshParameters::NX;ix++){
		cell[ix]->Calc_pi_target();
	}
}