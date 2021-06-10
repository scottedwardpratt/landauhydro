#include "eos.h"

CEoS::CEoS(CparameterMap *parmapset){
	parmap=parmapset;  
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
}

void CEoS::eos(double epsilon,double rhoB,double &T,double &Pr,double &SoverB,double &cs2){
	Pr=epsilon/1.5;
	T=Pr/rhoB;
	SoverB=log(pow(T,1.5)/rhoB);
	cs2=(5.0/3.0)*T/mass;
}
