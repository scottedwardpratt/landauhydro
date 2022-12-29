#include "landau.h"
#include "eos.h"
CEoS *CIntegralCell::eos=NULL;

CIntegralCell::CIntegralCell(){
	int i;
	pi_shear.resize(4);
	SE.resize(4);
	for(i=0;i<4;i++){
		SE[i].resize(4);
		pi_shear[i].resize(4);
	}
	Zero();
}


void CIntegralCell::Zero(){
	int i,j;
	S=Q=rho=sigma=alpha_Zeta=alpha_Kx=tau_Zeta=tau_K=gmma_Pr=Pi=epsilon=grad2Rho=Kx_target=0.0;
	alpha_eta=alpha_zeta=alpha_gmma=tau_zeta=tau_gmma=0.0;
	for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			SE[i][j]=0.0;
			pi_shear[i][j]=0.0;
		}
	}
}



void CIntegralCell::CalcEpsilonSE(){
	int i,j;
	double grad2rho,epsilon,gradrho2=0.0,mass=CEos::mass;
	vector<double> gradrho(NDIM+1);
	epsilon=Pdens[0];
	printf("-------\nix=%d: before, epsilon=%g\n",ix,epsilon);
	for(i=1;i<=NDIM;i++)
		epsilon-=(M[i]*u[i]+0.5*CEos::mass*rho*u[i]*u[i]);
	printf("now, epsilon=%g\n",epsilon);
	grad2rho=Grad2RhoB();
	CalcGradRhoB(gradrho);
	for(i=1;i<=NDIM;i++){
		gradrho2+=gradrho[i]*gradrho[i];
	}
	epsilonk=epsilon+0.5*CEos::kappa*rho*grad2rho;
	if(epsilonk<0.0){
		printf("Yikes!!! epsilon=%g, epsilonk=%g, rho=%g\n",epsilon,epsilonk,rho);
		exit(1);
	}
	CEos::CalcEoS_of_rho_epsilon(this);
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
			SE[i][j]=pi_shear[i][j];
			SE[i][j]+=CEos::kappa*gradrho[i]*gradrho[j];
			if(i==j)
				SE[i][j]+=Pr-CEos::kappa*(rho*grad2rho+0.5*gradrho2);
			SE[i][j]+=mass*rho*u[i]*u[j];
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

void CIntegralCell::CalcGrad2Rho(){
	double dxplus,dxminus,gradrhoplus,gradrhominus;
	dxplus=0.5*neighorPlus->Delx+0.5*Delx;
	dxminus=0.5*neighborMinus->Delx+0.5*Delx;
	gradrhoplus=(neighborPlus->rho-rho)/dxplus;
	gradrhominus=(rho-neighborMinus->rho)/dxminus;
	gradrhoPlus=(gradrhoplus-gradrhominus)/Delx;
}

void CIntegralCell::UpdateBulkQuantities(){
	CEos::CalcEoS_of_rho_sdens(this);
	CEos::CalcEtaZetaK(this);
}

// Half-Integral Cell

CHalfIntegralCell::CHalfIntegralCell(){
	int i;
	pi_shear_target.resize(4);
	for(i=0;i<4;i++){
		pi_shear_target[i].resize(4);
	}
}

void CHalfIntegralCell::Zero(){
	vx=Kx=pi_bulk_target=0.0;
	for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			pi_shear_target[i][j]=0.0;
		}
	}
}

void CHalfIntegralCell::GetOmega(){
	omega=0.5*(neighborPlus->Vx-neighborMinus->Vx)/Delx;
	omega=4.0*omega/3.0;
}

void CHalfIntegralCell::GetDelDotV(){
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

void CHalfIntegralCell::Calc_pi_shear_target(){
	int i,j;
	double deldotv=DelDotU();
	for(i=1;i<NDIM;i++){
		for(j=i;j<NDIM+1;j++){
			if(j!=i){
				pi_shear_target[i][j]=-eta*(neighborPlus[j]->u[i]-neighborMinus[j]->u[i])/(2.0*Delx);
				pi_shear_target[i][j]-=eta*(neighborPlus[i]->u[j]-neighborMinus[i]->u[j])/(2.0*Delx);
				pi_shear_target[j][i]=pi_shear_target[i][j];
			}
		}
		pi_shear_target[i][i]+=(2.0*eta/3.0)*deldotv;
	}
}

void CHalfIntegralCell::Zero(){
	vx=Kx=0.0;
}


void CHalfIntegralCell::PrintInfo(){
	printf("vx=%g, Kx=%g\n",vx,Kx);
}

