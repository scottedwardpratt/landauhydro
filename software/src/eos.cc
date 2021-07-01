#include "eos.h"
double CEoS::kappa=0.0;
double CEoS::mass=0.0;
double CEoS::Kfactor=0.0;

CEoS::CEoS(CparameterMap *parmapset){
	parmap=parmapset;  
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	Kfactor=parmap->getD("EOS_KFACTOR",0.01);
}

CEoS_FreeGas::CEoS_FreeGas(CparameterMap *parmapset){
	parmap=parmapset; 
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	Kfactor=parmap->getD("EOS_KFACTOR",1.0);
}

CEoS_VdW::CEoS_VdW(CparameterMap *parmapset){
	parmap=parmapset;  
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	Kfactor=parmap->getD("EOS_KFACTOR",1.0);
	a=parmap->getD("EOS_VDW_A",1.0);
	rho0=parmap->getD("EOS_VDW_RHO0",1.0);
}

void CEoS_FreeGas::CalcEoS(CLandauCell *cell){
	double epsilonk=cell->epsilonk,rhoB=cell->jB[0];
	double T,Pr,SoverB,cs2,K;
	Pr=epsilonk/1.5;
	T=Pr/rhoB;
	if(T<0.0){
		printf("in Free Gas CalcEoS, T<0!! =%g\n",T);
		exit(1);
	}
	if(rhoB<0.0){
		printf("in Free Gas CalcEoS, rhoB<0!! =%g\n",rhoB);
		exit(1);
	}
	SoverB=log(pow(T,1.5)/rhoB);
	cs2=(5.0/3.0)*T/mass;
	K=Kfactor*sqrt(T/mass);
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;
	cell->K=K;
}

void CEoS_VdW::CalcEoS(CLandauCell *cell){
	double epsilonk=cell->epsilonk,rhoB=cell->jB[0];
	double T,Pr,SoverB,cs2,K;
	T=2.0*(epsilonk+a*rhoB*rhoB)/(3.0*rhoB);
	if(T<0.0){
		printf("in VdW CalcEoS, T<0!! =%g\n",T);
		exit(1);
	}
	if(rhoB<0.0){
		printf("in VdW CalcEoS, rhoB<0!! =%g\n",rhoB);
		exit(1);
	}
	Pr=rhoB*T/(1.0-rhoB/rho0)-a*rhoB*rhoB;
	SoverB=1.5*log(T)+log((rho0/rhoB)-1.0);
	cs2=((5.0*T/3.0)/((1.0-rhoB/rho0)*(1.0-rhoB/rho0))-2.0*a*rhoB)/mass;
	K=Kfactor*sqrt(T/mass);
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;
	cell->K=K;
}
