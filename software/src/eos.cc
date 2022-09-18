#include "eos.h"
#include "msu_commonutils/misc.h"

CEoS::CEoS(CparameterMap *parmapset){
	parmap=parmapset;  
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	Kfactor=parmap->getD("EOS_KFACTOR",0.2);
	etafactor=parmap->getD("EOS_ETAFACTOR",1.0);
	zetafactor=parmap->getD("EOS_ZETAFACTOR",0.0);
}

CEoS_FreeGas::CEoS_FreeGas(CparameterMap *parmapset){
	parmap=parmapset;
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	Kfactor=parmap->getD("EOS_KFACTOR",1.0);
	etafactor=parmap->getD("EOS_ETAFACTOR",1.0);
	zetafactor=parmap->getD("EOS_ZETAFACTOR",0.0);
}

CEoS_VdW::CEoS_VdW(CparameterMap *parmapset){
	parmap=parmapset;
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	Kfactor=parmap->getD("EOS_KFACTOR",1.0);
	a=parmap->getD("EOS_VDW_A",1.0);
	rho0=parmap->getD("EOS_VDW_RHO0",1.0);
	etafactor=parmap->getD("EOS_ETAFACTOR",1.0);
	zetafactor=parmap->getD("EOS_ZETAFACTOR",0.0);
}

CEoS_Scott::CEoS_Scott(CparameterMap *parmapset){
	parmap=parmapset;
	kappa=parmap->getD("EOS_KAPPA",0.0);
	mass=parmap->getD("EOS_MASS",1.0);
	Kfactor=parmap->getD("EOS_KFACTOR",1.0);
	a=parmap->getD("EOS_VDW_A",1.0);
	rho0=parmap->getD("EOS_VDW_RHO0",1.0);
	etafactor=parmap->getD("EOS_ETAFACTOR",1.0);
	zetafactor=parmap->getD("EOS_ZETAFACTOR",0.0);
}

void CEoS_FreeGas::CalcEoS_of_rho_epsilon(CLandauCell *cell){
	double epsilonk=cell->epsilonk,rhoB=cell->jB[0];
	double T,Pr,SoverB,cs2;
	Pr=epsilonk/1.5;
	T=Pr/rhoB;
	if(T<0.0){
		printf("in Free Gas CalcEoS, T<0!! =%g\n",T);
		printf("epsilonk=%g, rhoB=%g\n",epsilonk,rhoB);
		exit(1);
	}
	if(rhoB<0.0){
		printf("in Free Gas CalcEoS, rhoB<0!! =%g\n",rhoB);
		exit(1);
	}
	SoverB=log(pow(T,1.5)/rhoB);
	cs2=(5.0/3.0)*T/mass;
	
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;
	
	CalcEtaZetaK(cell);
}

void CEoS_FreeGas::CalcEoS_of_rho_T(CLandauCell *cell){
	double rhoB=cell->jB[0];
	cell->epsilonk=1.5*rhoB*cell->T;
	CalcEoS_of_rho_epsilon(cell);
}

void CEoS_VdW::CalcEoS_of_rho_epsilon(CLandauCell *cell){
	double epsilonk=cell->epsilonk,rhoB=cell->jB[0];
	double T,Pr,SoverB,cs2;
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
	
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;

	CalcEtaZetaK(cell);
}

void CEoS_VdW::CalcEoS_of_rho_T(CLandauCell *cell){
	double rhoB=cell->jB[0];
	cell->epsilonk=1.5*rhoB*cell->T-a*rhoB*rhoB;
	CalcEoS_of_rho_epsilon(cell);
}

void CEoS_Scott::CalcEoS_of_rho_epsilon(CLandauCell *cell){
	double epsilonk=cell->epsilonk,rhoB=cell->jB[0];
	double T,Pr,SoverB,cs2;
	double F,Fprime,G,Gprime,H,x;
	x=rhoB/rho0;
	// P=rho*T*F(x)-a*rho0*rho*G(x)
	F=(1.0+0.5*x);
	G=x*x*(1.0-x)/((1.0+x)*(1.0+x));
	//H= x*\int dx' G(x')/x'^2
	H=(2.0*x*x/(1.0+x))-x*log(1.0+x);
	T=(2.0/(3.0*rhoB))*(epsilonk+a*rho0*rho0*H);
	if(T<0.0){
		printf("in Scott CalcEoS, T<0!! =%g\n",T);
		printf("epsilonk=%g, rhoB=%g\n",epsilonk,rhoB);
		exit(1);
	}
	if(rhoB<0.0){
		printf("in Scott CalcEoS, rhoB<0!! =%g\n",rhoB);
		exit(1);
	}
	Pr=rhoB*T*F-a*rho0*rho0*G;
	SoverB=1.5*log(T)-log(rhoB)-0.5*x;
	Fprime=0.5/rho0;
	Gprime=(-x*x*x-3.0*x*x+2.0*x)/pow(1.0+x,3);
	cs2=(T/mass)*(F+x*Fprime+2.0*F*F/3.0)-a*rho0*Gprime/mass;
	
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;

	CalcEtaZetaK(cell);
}

void CEoS_Scott::CalcEoS_of_rho_T(CLandauCell *cell){
	double x,H,rhoB=cell->jB[0];
	x=rhoB/rho0;
	//H= x*\int dx' G(x')/x'^2
	H=(2.0*x*x/(1.0+x))-x*log(1.0+x);
	cell->epsilonk=1.5*rhoB*cell->T-a*rho0*rho0*H;
	CalcEoS_of_rho_epsilon(cell);
}

void CEoS::CalcEtaZetaK(CLandauCell *cell){
	double rhoB=cell->jB[0],T=cell->T;

	cell->tau_K=Kfactor*sqrt(mass/T)/rhoB;
	cell->alpha_K=sqrt(rhoB*T*T*T);
	cell->K=(21.0/4.0)*(T/mass)*cell->tau_K;

	cell->tau_eta=etafactor*sqrt(mass/T)/rhoB;
	cell->alpha_eta=(4.0/15.0)*sqrt(rhoB*rhoB*T*T);
	cell->eta=cell->alpha_eta*cell->tau_eta;//// DANGER -- CHECK!!!!

	cell->tau_zeta=1.0;
	cell->alpha_zeta=1.0;
	cell->zeta=0.0;

}