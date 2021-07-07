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
	kflow.resize(NDIM+1);
	kflow_target.resize(NDIM+1);
	neighborMinus.resize(NDIM+1);
	neighborPlus.resize(NDIM+1);
	for(i=0;i<=NDIM;i++){
		SE[i].resize(NDIM+1);
		for(j=0;j<=NDIM;j++){
			SE[i][j]=0.0;
		}
	}
	Zero();
}

void CLandauCell::Zero(){
	int i,j;
	epsilonk=Pr=tau_K=0.0;
	for(i=0;i<=NDIM;i++){
		u[i]=M[i]=jB[i]=Pdens[i]=0.0;
		for(j=0;j<=NDIM;j++){
			SE[i][j]=0.0;
		}
	}
}

void CLandauCell::PrintInfo(){
	printf("jB=(%g,%g,%g,%g), epsilonk=%g\n",jB[0],jB[1],jB[2],jB[3],epsilonk);
	printf("u=(%g,%g,%g), T=%g, Pr=%g, cs2=%g, SoverB=%g, kflow=(%g,%g,%g)\n",
	u[1],u[2],u[3],T,Pr,cs2,SoverB,kflow[1],kflow[2],kflow[3]);
}

double CLandauCell::Grad2RhoB(){
	double answer=0.0;
	int i;
	for(i=1;i<=NDIM;i++){
		answer+=neighborPlus[i]->jB[0]+neighborMinus[i]->jB[0]-2.0*jB[0];
	}
	return answer/(DXYZ*DXYZ);
}

void CLandauCell::CalcGradRhoB(vector<double> &GradRhoB){
	int i;
	for(i=1;i<=NDIM;i++){
		GradRhoB[i]=0.5*(neighborPlus[i]->jB[0]-neighborMinus[i]->jB[0])/DXYZ;
	}
}

void CLandauCell::CalcDeliTij(vector<double> &DeliTij){
	int i,j;
	for(i=0;i<=NDIM;i++){
		DeliTij[i]=0.0;
		for(j=1;j<=NDIM;j++){			
			DeliTij[i]+=0.5*(neighborPlus[j]->SE[i][j]-neighborMinus[j]->SE[i][j])/DXYZ;
		}
	}
}

double CLandauCell::DelDotU(){
	int i;
	double answer=0.0;
	for(i=1;i<=NDIM;i++){
		answer+=neighborPlus[i]->u[i]-neighborMinus[i]->u[i];
	}
	return answer/(2.0*DXYZ);
}

void CLandauCell::Calckflow_target(){
	int i;
	double Kflow_target;
	for(i=1;i<=NDIM;i++){
		Kflow_target=CEoS::Kfactor*(35.0/4.0)*(T/CEoS::mass)*(neighborPlus[i]->T-neighborMinus[i]->T)/(2.0*DXYZ);	
		kflow_target[i]=Kflow_target/(sigma_K*jB[0]);
	}
}

double CLandauCell::CalcDivKFlow(){
	int i;
	double DivKFlow=0.0,Kplus,Kminus;
	CLandauCell *neighbor;
	for(i=1;i<=NDIM;i++){
		neighbor=neighborPlus[i];
		Kplus=neighbor->kflow[i]*neighbor->sigma_K*neighbor->jB[0];
		neighbor=neighborMinus[i];
		Kminus=neighbor->kflow[i]*neighbor->sigma_K*neighbor->jB[0];
		DivKFlow+=(Kplus-Kminus);
	}
	return DivKFlow/(2.0*DXYZ);
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
		M[i]=-0.5*eos->kappa*jB[0]*jB[0]*(neighborPlus[i]->DelDotU()-neighborMinus[i]->DelDotU())/(2.0*DXYZ);
		M[i]+=0.5*eos->kappa*jB[0]*DelDotU()*GradRhoB[i];
		for(j=1;j<=NDIM;j++){
			M[i]-=0.5*eos->kappa*jB[0]*GradRhoB[j]
				*(neighborPlus[i]->u[j]-neighborMinus[i]->u[j]+neighborPlus[j]->u[i]-neighborMinus[j]->u[i])
					/(2.0*DXYZ);
		}
	}
}

void CLandauCell::CalcEpsilonSE(){
	int i,j;
	double grad2rhoB,epsilon,gradrhoB2=0.0,mass=eos->mass;
	vector<double> gradrhoB(NDIM+1);
	epsilon=Pdens[0];
	for(i=1;i<=NDIM;i++)
		epsilon-=(M[i]*u[i]+0.5*eos->mass*jB[0]*u[i]*u[i]);
	grad2rhoB=Grad2RhoB();
	CalcGradRhoB(gradrhoB);
	for(i=1;i<=NDIM;i++){
		gradrhoB2+=gradrhoB[i]*gradrhoB[i];
	}
	epsilonk=epsilon+0.5*eos->kappa*jB[0]*grad2rhoB;
	eos->CalcEoS(this);
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
			SE[i][j]=eos->kappa*gradrhoB[i]*gradrhoB[j];
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
