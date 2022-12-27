#include "landau.h"
#include "eos.h"
double CLandauCell::DXYZ=0.0;
int CLandauCell::NDIM=3;
CEoS *CLandauCell::eos=NULL;
double CLandauCell::Tlowest=999999999999.0;
double CLandauCell::Thighest=0.00000;

CLandauCell::CLandauCell(){
	int i,j;
	u.resize(NDIM+1);
	Pdens.resize(NDIM+1);
	M.resize(NDIM+1);
	jB.resize(NDIM+1);
	SE.resize(NDIM+1);
	pivisc.resize(NDIM+1);
	pi_target.resize(NDIM+1);
	Kflow.resize(NDIM+1);
	Kflow_target.resize(NDIM+1);
	neighborMinus.resize(NDIM+1);
	neighborPlus.resize(NDIM+1);
	for(i=0;i<=NDIM;i++){
		Kflow[i]=0.0;
		Pdens[i]=0.0;
		SE[i].resize(NDIM+1);
		pi_target[i].resize(NDIM+1);
		pivisc[i].resize(NDIM+1);
		for(j=0;j<=NDIM;j++){
			SE[i][j]=0.0;
			pivisc[i][j]=0.0;
		}
	}
	Zero();
}

void CIntegralCell::Zero(){
	S=Q=rho=sigma=alphaZeta=alphaK=tauZta=tauK=gammaPr=Pi=epsilon=grad2Rho=0.0;
}

void CHalfIntegralCell::Zero(){
	vx=Kx=0.0;
}

void CHalfIntegralCell::PrintInfo(){
	printf("vx=%g, Kx=%g\n",vx,Kx);
}

void CIntegralCell::CalcGrad2Rho(){
	double dxplus,dxminus,gradrhoplus,gradrhominus;
	dxplus=0.5*neighorPlus->DelX+0.5*DelX;
	dxminus=0.5*neighborMinus->DelX+0.5*DelX;
	gradrhoplus=(neighborPlus->rho-rho)/dxplus;
	gradrhominus=(rho-neighborMinus->rho)/dxminus;
	gradrhoPlus=(gradrhoplus-gradrhominus)/DelX;
}

void CIntegralCell::UpdateBulkQuantities(){
	eos->CalcEoS_of_rho_sdens(this);
	eos->CalcEtaZetaK(this);
}

void CHalfIntegralCell::GetOmega(double &omega){
	double DX=CMeshParameters::DX;
	omega=0.5*(neighborPlus->Vx-neighborMinus->Vx)/DX;
	omega=4.0*omega/3.0;
}
void CHalfIntegralCell::GetDelDotV(double &deldotv){
	double DX=CMeshParameters::DX;
	deldotv=0.5*(neighborPlus->Vx-neighborMinus->Vx)/DX;
}


void CHalfIntegralCell::CalcDxTxx(vector<double> &DxTxx){
	DxTxx=0.5*(neighborPlus->SE[1][1]-neighborMinus->SE[1][1])/DXYZ;
}

double CLandauCell::DelDotU(CIntegralCell *newcell,CIntegralCell *oldcell){
	double answer=newcell->neighborPlus->Vx-newcellneighborMinus->Vx;
	answer+=oldcell->neighborPlus->Vx-oldcell->neighborMinus->Vx;
	return answer/(2.0*DXYZ);
}

void CHalfIntegralCell::Calc_Kflow_target(){
		K_target=-K*(neighborPlus->T-neighborMinus->T)/(2.0*DXYZ);
	}

double CLandauCell::CalcDivKFlow(){
	int i;
	double DivKFlow=0.0;
	for(i=1;i<=NDIM;i++){
		DivKFlow+=(neighborPlus[i]->Kflow[i]-neighborMinus[i]->Kflow[i]);
	}
	return DivKFlow/(2.0*DXYZ);
}

void CIntegralCell::CalcEpsilonSE(){
	int i,j;
	double grad2rhoB,epsilon,gradrhoB2=0.0,mass=eos->mass;
	vector<double> gradrhoB(NDIM+1);
	epsilon=Pdens[0];
	printf("-------\nix=%d: before, epsilon=%g\n",ix,epsilon);
	for(i=1;i<=NDIM;i++)
		epsilon-=(M[i]*u[i]+0.5*eos->mass*jB[0]*u[i]*u[i]);
	printf("now, epsilon=%g\n",epsilon);
	grad2rhoB=Grad2RhoB();
	CalcGradRhoB(gradrhoB);
	for(i=1;i<=NDIM;i++){
		gradrhoB2+=gradrhoB[i]*gradrhoB[i];
	}
	epsilonk=epsilon+0.5*eos->kappa*jB[0]*grad2rhoB;
	if(epsilonk<0.0){
		printf("Yikes!!! epsilon=%g, epsilonk=%g, jB[0]=%g\n",epsilon,epsilonk,jB[0]);
		exit(1);
	}
	eos->CalcEoS_of_rho_epsilon(this);
	if(T<Tlowest){
		Tlowest=T;
		if(Tlowest<0.0){
			printf("Tlowest=%g\n",T);
			exit(1);
		}
	}
	if(T>Thighest){
		Thighest=T;
	}
	for(i=1;i<=NDIM;i++){
		for(j=1;j<=NDIM;j++){
			SE[i][j]=pivisc[i][j];
			SE[i][j]+=eos->kappa*gradrhoB[i]*gradrhoB[j];
			if(i==j)
				SE[i][j]+=Pr-eos->kappa*(jB[0]*grad2rhoB+0.5*gradrhoB2);
			SE[i][j]+=mass*jB[0]*u[i]*u[j];
		}
	}
	for(i=1;i<=NDIM;i++){
		SE[0][i]=epsilon*u[i];
		for(j=1;j<=NDIM;j++){
			SE[0][i]+=SE[i][j]*u[j];
		}
		SE[i][0]=SE[0][i];
	}
}

void CHalfIntegralCell::Calc_pi_target(){
	int i,j;
	double deldotv=DelDotU();
	for(i=1;i<NDIM;i++){
		for(j=i;j<NDIM+1;j++){
			if(j!=i){
				pi_target[i][j]=-eta*(neighborPlus[j]->u[i]-neighborMinus[j]->u[i])/(2.0*DXYZ);
				pi_target[i][j]-=eta*(neighborPlus[i]->u[j]-neighborMinus[i]->u[j])/(2.0*DXYZ);
				pi_target[j][i]=pi_target[i][j];
			}
		}
		pi_target[i][i]+=(2.0*eta/3.0)*deldotv;
	}
}
