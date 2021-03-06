#include "eos.h"
double CEoS::kappa=0.0;
double CEoS::mass=0.0;
double CEoS::Kfactor=0.0;
double CEoS::Etafactor=0.0;
double CEoS::Zetafactor=0.0;

CEoS::CEoS(CparameterMap *parmapset){
	parmap=parmapset;  
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	Kfactor=parmap->getD("EOS_KFACTOR",0.2);
	Etafactor=parmap->getD("EOS_ETAFACTOR",1.0);
	Zetafactor=parmap->getD("EOS_ZETAFACTOR",0.0);
}

CEoS_FreeGas::CEoS_FreeGas(CparameterMap *parmapset){
	parmap=parmapset;
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	Kfactor=parmap->getD("EOS_KFACTOR",1.0);
	Etafactor=parmap->getD("EOS_ETAFACTOR",1.0);
	Zetafactor=parmap->getD("EOS_ZETAFACTOR",0.0);
}

CEoS_VdW::CEoS_VdW(CparameterMap *parmapset){
	parmap=parmapset;
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	Kfactor=parmap->getD("EOS_KFACTOR",1.0);
	a=parmap->getD("EOS_VDW_A",1.0);
	rho0=parmap->getD("EOS_VDW_RHO0",1.0);
	Etafactor=parmap->getD("EOS_ETAFACTOR",1.0);
	Zetafactor=parmap->getD("EOS_ZETAFACTOR",0.0);
}

void CEoS_FreeGas::CalcEoS_of_rho_epsilon(CLandauCell *cell){
	double epsilonk=cell->epsilonk,rhoB=cell->jB[0];
	double T,Pr,SoverB,cs2,K,tau_K,sigma_K,eta,tau_eta,sigma_eta,zeta,tau_zeta,sigma_zeta;
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
	
	sigma_K=sqrt(35.0*rhoB/(4.0*mass*T));
	tau_K=Kfactor*sqrt(mass/T)/rhoB;
	K=(sigma_K*sigma_K/(T*T))*tau_K; // usual conductivity
	K=K/(sigma_K*rhoB);   // scaled viscosity
	
	sigma_eta=sqrt(rhoB*T*T);
	tau_eta=Etafactor*sqrt(mass/T)/rhoB;
	eta=(sigma_eta*sigma_eta/T)*tau_eta; // usual viscosity
	eta=eta/(sigma_eta*rhoB);  // scaled viscosity
	
	sigma_zeta=1.0;
	tau_zeta=0.0;
	zeta=0.0;
	
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;
	
	cell->K=K;
	cell->tau_K=tau_K;
	cell->sigma_K=sigma_K;
	
	cell->eta=eta;
	cell->tau_eta=tau_eta;
	cell->sigma_eta=sigma_eta;
	
	cell->zeta=zeta;
	cell->sigma_zeta=sigma_zeta;
	cell->tau_zeta=tau_zeta;
}

void CEoS_FreeGas::CalcEoS_of_rho_T(CLandauCell *cell){
	double rhoB=cell->jB[0];
	cell->epsilonk=1.5*rhoB*cell->T;
	CalcEoS_of_rho_epsilon(cell);
}

void CEoS_VdW::CalcEoS_of_rho_epsilon(CLandauCell *cell){
	double epsilonk=cell->epsilonk,rhoB=cell->jB[0];
	double T,Pr,SoverB,cs2,K,tau_K,sigma_K,eta,tau_eta,sigma_eta,zeta,tau_zeta,sigma_zeta;
	T=2.0*(epsilonk+a*rhoB*rhoB)/(3.0*rhoB);
	if(T<0.0){
		printf("in VdW CalcEoS, T<0!! =%g\n",T);
		printf("epsilonk=%g, rhoB=%g\n",epsilonk,rhoB);
		exit(1);
	}
	if(rhoB<0.0){
		printf("in VdW CalcEoS, rhoB<0!! =%g\n",rhoB);
		exit(1);
	}
	Pr=rhoB*T/(1.0-rhoB/rho0)-a*rhoB*rhoB;
	SoverB=1.5*log(T)+log((rho0/rhoB)-1.0);
	cs2=((5.0*T/3.0)/((1.0-rhoB/rho0)*(1.0-rhoB/rho0))-2.0*a*rhoB)/mass;
	
	sigma_K=sqrt(35.0*rhoB*T*T*T/(4.0*mass));
	tau_K=sqrt(mass/T)/pow(rhoB,1.0/3.0);
	K=Kfactor*(sigma_K*sigma_K/(T*T))*tau_K; // usual conductivity
	K=K/(sigma_K*rhoB);   // scaled viscosity
	
	sigma_eta=sqrt(rhoB*T*T);
	tau_eta=sqrt(mass/T)/pow(rhoB,1.0/3.0);
	eta=(sigma_eta*sigma_eta/T)*tau_eta; // usual viscosity
	eta=Etafactor*eta/(sigma_eta*rhoB);  // scaled viscosity
	
	sigma_zeta=1.0;
	tau_zeta=0.0;
	zeta=0.0;
	
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;
	
	cell->K=K;
	cell->tau_K=tau_K;
	cell->sigma_K=sigma_K;
	
	cell->eta=eta;
	cell->tau_eta=tau_eta;
	cell->sigma_eta=sigma_eta;
	
	cell->zeta=zeta;
	cell->sigma_zeta=sigma_zeta;
	cell->tau_zeta=tau_zeta;
}

void CEoS_VdW::CalcEoS_of_rho_T(CLandauCell *cell){
	double rhoB=cell->jB[0];
	cell->epsilonk=1.5*rhoB*cell->T-a*rhoB*rhoB;
	CalcEoS_of_rho_epsilon(cell);
}
