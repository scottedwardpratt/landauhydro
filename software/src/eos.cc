#include "eos.h"

CEoS::CEoS(CparameterMap *parmapset){
	parmap=parmapset;  
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
}

CEoS_FreeGas::CEoS_FreeGas(CparameterMap *parmapset){
	parmap=parmapset;  
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
}

CEoS_VdW::CEoS_VdW(CparameterMap *parmapset){
	parmap=parmapset;  
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	a=parmap->getD("EOS_VDW_A",1.0);
	rho0=parmap->getD("EOS_VDW_RHO0",1.0);
}

void CEoS_FreeGas::eos(double epsilon,double rhoB,double &T,double &Pr,double &SoverB,double &cs2){
	Pr=epsilon/1.5;
	T=Pr/rhoB;
	SoverB=log(pow(T,1.5)/rhoB);
	cs2=(5.0/3.0)*T/mass;
}

void CEoS_VdW::eos(double epsilon,double rhoB,double &T,double &Pr,double &SoverB,double &cs2){
	T=2.0*(epsilon+a*rhoB*rhoB)/(3.0*rhoB);
	Pr=rhoB*T/(1.0-rhoB/rho0)-a*rhoB*rhoB;
	SoverB=1.5*log(T)+log((rho0/rhoB)-1.0);
	cs2=((5.0*T/3.0)/((1.0-rhoB/rho0)*(1.0-rhoB/rho0))-2.0*a*rhoB)/mass;
}
