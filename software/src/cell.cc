#include "landau.h"
#include "eos.h"
double CLandauCell::DXYZ=0.0;
int CLandauCell::NDIM=3;
CEoS *CLandauCell::eos=NULL;

CLandauCell::CLandauCell(){
	int i;
	u.resize(NDIM+1);
	M.resize(NDIM+1);
	jB.resize(NDIM+1);
	SE.resize(NDIM+1);
	neighborMinus.resize(NDIM+1);
	neighborPlus.resize(NDIM+1);
	for(i=0;i<=NDIM;i++){
		SE[i].resize(NDIM+1);
	}
	Zero();
}

void CLandauCell::Zero(){
	int i,j;
	epsilon=Pr=0.0;
	for(i=0;i<=NDIM;i++){
		u[i]=M[i]=jB[i]=0.0;
		for(j=0;j<=NDIM;j++){
			SE[i][j]=0.0;
		}
	}
}

void CLandauCell::PrintInfo(){
	printf("rhoB=%g, epsilon=%g, u=(%g,%g,%g)\n",jB[0],epsilon,u[1],u[2],u[3]);
}

double CLandauCell::Grad2RhoB(){
	double answer=0.0;
	int i;
	for(i=1;i<=NDIM;i++){
		answer+=neighborPlus[i]->jB[0]+neighborMinus[i]->jB[0]-2.0*jB[0];
	}
	return answer/(2.0*DXYZ*DXYZ);
}

void CLandauCell::CalcGradRhoB(vector<double> &GradRhoB){
	int i;
	for(i=1;i<=NDIM;i++){
		GradRhoB[i]=0.5*(neighborPlus[i]->jB[0]-neighborMinus[i]->jB[0])/DXYZ;
	}
}

void CLandauCell::CalcDeliTij(vector<double> &DeliTij){
	int i,j;
	for(i=1;i<=NDIM;i++){
		DeliTij[i]=0.0;
		for(j=1;j<=NDIM;j++){
			DeliTij[i]+=0.5*(neighborPlus[j]->SE[i][j]-neighborMinus[j]->SE[i][j])/DXYZ;
		}
	}
}

double CLandauCell::DelDotT0i(){
	int i;
	double answer=0.0;
	for(i=1;i<=NDIM;i++)
		answer+=neighborPlus[i]->SE[0][i]-neighborMinus[i]->SE[0][i];
	answer=answer/(2.0*DXYZ);
	return answer;
}

double CLandauCell::DelDotU(){
	int i;
	double answer=0.0;
	for(i=1;i<=NDIM;i++){
		answer+=neighborPlus[i]->u[i]-neighborMinus[i]->u[i];
	}
	return answer/(2.0*DXYZ);
}

double CLandauCell::DelDotJB(){
	int i;
	double answer=0.0;
	for(i=1;i<=NDIM;i++){
		answer+=neighborPlus[i]->jB[i]-neighborMinus[i]->jB[i];
	}
	return answer/(2.0*DXYZ);
}

void CLandauCell::CalcM(){
	int i,j;
	vector<double> GradRhoB(NDIM+1);
	CalcGradRhoB(GradRhoB);
	for(i=1;i<=NDIM;i++){
		M[i]=0.0;
		M[i]=0.5*eos->kappa*jB[0]*jB[0]*(neighborPlus[i]->DelDotU()-neighborMinus[i]->DelDotU())/(2.0*DXYZ);
		M[i]+=0.5*eos->kappa*jB[0]*DelDotU()*GradRhoB[i];
		for(j=1;j<=NDIM;j++){
			M[i]-=0.5*eos->kappa*jB[0]*GradRhoB[j]
				*(neighborPlus[i]->u[j]-neighborMinus[i]->u[j]+neighborPlus[j]->u[i]-neighborMinus[j]->u[i])
					/(2.0*DXYZ);
		}
	}
}

void CLandauCell::CalcEpsilonU(){
	int i;
	double mass=eos->mass,vsquared=0.0;
	CalcM();
	for(i=1;i<=NDIM;i++){
		u[i]=(SE[0][i]-M[i])/(jB[0]*mass);
		jB[i]=u[i]*jB[0];
		vsquared+=u[i]*u[i];
	}
	epsilon=SE[0][0]-jB[0]*mass*(1.0+0.5*vsquared);
	for(i=1;i<=NDIM;i++)
		epsilon-=2.0*M[i]*u[i];
}

void CLandauCell::CalcTij(){
	int i,j;
	double vsquared=0.0,epsilonk,grad2rhoB,gradrhoB2=0.0;
	vector<double> gradrhoB(NDIM+1);
	grad2rhoB=Grad2RhoB();
	CalcGradRhoB(gradrhoB);
	for(i=1;i<=NDIM;i++){
		gradrhoB2+=gradrhoB[i]*gradrhoB[i];
	}
	epsilonk=epsilon+0.5*jB[0]*grad2rhoB;
	eos->eos(epsilonk,jB[0],T,Pr);
	for(i=1;i<=NDIM;i++)
		vsquared+=u[i]*u[i];
	
	for(i=1;i<=NDIM;i++){
		for(j=i;j<=NDIM;j++){
			SE[i][j]=eos->mass*u[i]*u[j]+eos->kappa*gradrhoB[i]*gradrhoB[j];
			if(i==j)
				SE[i][j]+=Pr-eos->kappa*(jB[0]*grad2rhoB+0.5*gradrhoB2);
			else
				SE[j][i]=SE[i][j];
		}
	}	
}