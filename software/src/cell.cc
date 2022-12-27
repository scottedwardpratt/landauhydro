#include "landau.h"
#include "eos.h"
CEoS *CIntegralCell::eos=NULL;

CIntegralCell::CIntegralCell(){
	int i,j;
	pivisc.resize(4);
	pi_target.resize(4);
	alpha_eta=alpha_zeta=alpha_gamma=tau_zeta=tau_gamma=0..0;
	for(i=0;i<4;i++){
		SE[i].resize(4);
		pi_target[i].resize(4);
		pivisc[i].resize(4);
		for(j=0;j<4;j++){
			SE[i][j]=0.0;
			pivisc[i][j]=0.0;
		}
	}
	Zero();
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

void CIntegralCell::Zero(){
	S=Q=rho=sigma=alphaZeta=alphaK=tauZta=tauK=gammaPr=Pi=epsilon=grad2Rho=Kx_target0.0;
}

void CIntegralCell::CalcGrad2Rho(){
	double dxplus,dxminus,gradrhoplus,gradrhominus;
	dxplus=0.5*neighorPlus->Delx+0.5*Delx;
	dxminus=0.5*neighborMinus->Delx+0.5*Delx;
	gradrhoplus=(neighborPlus->rho-rho)/dxplus;
	gradrhominus=(rho-neighborMinus->rho)/dxminus;
	gradrhoPlus=(gradrhoplus-gradrhominus)/Delx;
}

void CIntegralCell::UpdateBulkQuantities(){
	eos->CalcEoS_of_rho_sdens(this);
	eos->CalcEtaZetaK(this);
}

// Half-Integral Cell


void CHalfIntegralCell::Zero(){
	vx=Kx=pi_bulk_target=0.0;
}

void CHalfIntegralCell::GetOmega(double &omega){
	omega=0.5*(neighborPlus->Vx-neighborMinus->Vx)/Delx;
	omega=4.0*omega/3.0;
}

void CHalfIntegralCell::GetDelDotV(double &deldotv){
	double DX=CMeshParameters::DX;
	deldotv=0.5*(neighborPlus->Vx-neighborMinus->Vx)/Delx;
}


void CHalfIntegralCell::CalcDxTxx(vector<double> &DxTxx){
	DxTxx=0.5*(neighborPlus->SE[1][1]-neighborMinus->SE[1][1])/Delx;
}

double CHalfIntegralCell::DelDotU(CIntegralCell *newcell,CIntegralCell *oldcell){
	double answer=newcell->neighborPlus->Vx-newcellneighborMinus->Vx;
	answer+=oldcell->neighborPlus->Vx-oldcell->neighborMinus->Vx;
	return answer/(2.0*Delx);
}

void CHalfIntegralCell::Calc_K_target(){
		K_target=-K*(neighborPlus->T-neighborMinus->T)/(2.0*Delx);
	}

double CHalfLandauCell::CalcDivK(){
	int i;
	double DivF=(neighborPlus->Kx-neighborMinus->Kx);
	return DivK/(2.0*Delx);
}

void CHalfIntegralCell::Calc_pi_target(){
	int i,j;
	double deldotv=DelDotU();
	for(i=1;i<NDIM;i++){
		for(j=i;j<NDIM+1;j++){
			if(j!=i){
				pi_target[i][j]=-eta*(neighborPlus[j]->u[i]-neighborMinus[j]->u[i])/(2.0*DX);
				pi_target[i][j]-=eta*(neighborPlus[i]->u[j]-neighborMinus[i]->u[j])/(2.0*DX);
				pi_target[j][i]=pi_target[i][j];
			}
		}
		pi_target[i][i]+=(2.0*eta/3.0)*deldotv;
	}
}


void CHalfIntegralCell::PrintInfo(){
	printf("vx=%g, Kx=%g\n",vx,Kx);
}

