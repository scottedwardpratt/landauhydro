#include "eos.h"

CEoS::CEoS(CparameterMap *parmapset){
	parmap=parmapset;  
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
}

void CEoS::eos(double epsilon,double rhoB,double &T,double &Pr){
	Pr=(epsilon-rhoB*mass)/1.5;
	T=Pr/rhoB;
}
