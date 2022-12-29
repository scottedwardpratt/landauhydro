#include "eos.h"
#include "msu_commonutils/misc.h"

double CEoS::kappa=0.0;
double CEoS::mass=0.0;
double CEoS::Kfactor=0.0;
double CEoS::etafactor=0.0;
double CEoS::zetafactor=0.0;
double CEoS::gmma=0.0;

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

void CEoS_FreeGas::CalcEoS_of_rho_epsilon(CIntegralCell *cell){
	double epsilonk=cell->epsilonk,rho=cell->rho;
	double T,Pr,SoverB,cs2;
	Pr=epsilonk/1.5;
	T=Pr/rho;
	if(T<0.0){
		printf("in Free Gas CalcEoS, T<0!! =%g\n",T);
		printf("epsilonk=%g, rho=%g\n",epsilonk,rho);
		exit(1);
	}
	if(rho<0.0){
		printf("in Free Gas CalcEoS, rho<0!! =%g\n",rho);
		exit(1);
	}
	SoverB=log(pow(T,1.5)/rho);
	cs2=(5.0/3.0)*T/mass;
	
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;
	
	CalcEtaZetaK(cell);
}

void CEoS_FreeGas::CalcEoS_of_rho_T(CIntegralCell *cell){
	double rho=cell->rho;
	cell->epsilonk=1.5*rho*cell->T;
	CalcEoS_of_rho_epsilon(cell);
}

void CEoS_VdW::CalcEoS_of_rho_epsilon(CIntegralCell *cell){
	double epsilonk=cell->epsilonk,rho=cell->rho;
	double T,Pr,SoverB,cs2;
	T=2.0*(epsilonk+a*rho*rho)/(3.0*rho);
	if(T<0.0){
		printf("in VdW CalcEoS, T<0!! =%g\n",T);
		printf("epsilonk=%g, rho=%g\n",epsilonk,rho);
		exit(1);
	}
	if(rho<0.0){
		printf("in VdW CalcEoS, rho<0!! =%g\n",rho);
		exit(1);
	}
	Pr=rho*T/(1.0-rho/rho0)-a*rho*rho;
	SoverB=1.5*log(T)+log((rho0/rho)-1.0);
	cs2=((5.0*T/3.0)/((1.0-rho/rho0)*(1.0-rho/rho0))-2.0*a*rho)/mass;
	
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;

	CalcEtaZetaK(cell);
}

void CEoS_VdW::CalcEoS_of_rho_T(CIntegralCell *cell){
	double rho=cell->rho;
	cell->epsilonk=1.5*rho*cell->T-a*rho*rho;
	CalcEoS_of_rho_epsilon(cell);
}

void CEoS_VdW::CalcEoS_of_rho_sdens(CIntegralCell *cell){
	double rho=cell->rho;
	double x=rho/rho0;
	double SoverB=cell->sdens/rho;
	double T;
	T=pow((x/(1.0-x))*exp(SoverB),2.0/3.0);
	cell->Pr=rho*T/(1.0-rho/rho0)-a*rho*rho;
	SoverB=1.5*log(T)+log((rho0/rho)-1.0);
	cell->sdens=SoverB*rho;
	cell->cs2=((5.0*T/3.0)/((1.0-x)*(1.0-x))-2.0*a*rho)/mass;
	cell->T=T;
	cell->epsilonk=1.5*rho*cell->T-a*rho*rho;
}




















void CEoS_Scott::CalcEoS_of_rho_epsilon(CIntegralCell *cell){
	double epsilonk=cell->epsilonk,rho=cell->rho;
	double T,Pr,SoverB,cs2;
	double F,Fprime,G,Gprime,H,x;
	x=rho/rho0;
	// P=rho*T*F(x)-a*rho0*rho*G(x)
	F=(1.0+0.5*x);
	G=x*x*(1.0-x)/((1.0+x)*(1.0+x));
	//H= x*\int dx' G(x')/x'^2
	H=(2.0*x*x/(1.0+x))-x*log(1.0+x);
	T=(2.0/(3.0*rho))*(epsilonk+a*rho0*rho0*H);
	if(T<0.0){
		printf("in Scott CalcEoS, T<0!! =%g\n",T);
		printf("epsilonk=%g, rho=%g\n",epsilonk,rho);
		exit(1);
	}
	if(rho<0.0){
		printf("in Scott CalcEoS, rho<0!! =%g\n",rho);
		exit(1);
	}
	Pr=rho*T*F-a*rho0*rho0*G;
	SoverB=1.5*log(T)-log(rho)-0.5*x;
	Fprime=0.5/rho0;
	Gprime=(-x*x*x-3.0*x*x+2.0*x)/pow(1.0+x,3);
	cs2=(T/mass)*(F+x*Fprime+2.0*F*F/3.0)-a*rho0*Gprime/mass;
	
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;

	CalcEtaZetaK(cell);
}

void CEoS_Scott::CalcEoS_of_rho_T(CIntegralCell *cell){
	double x,H,rho=cell->rho;
	x=rho/rho0;
	//H= x*\int dx' G(x')/x'^2
	H=(2.0*x*x/(1.0+x))-x*log(1.0+x);
	cell->epsilonk=1.5*rho*cell->T-a*rho0*rho0*H;
	CalcEoS_of_rho_epsilon(cell);
}

void CEoS::CalcEtaZetaK(CIntegralCell *cell){
	double rho=cell->rho,T=cell->T;

	cell->tau_gmma=Kfactor*sqrt(mass/T)/rho;
	cell->alpha_gmma=sqrt(rho*T*T*T);
	cell->gmma=(21.0/4.0)*(T/mass)*cell->tau_gmma;

	cell->tau_eta=etafactor*sqrt(mass/T)/rho;
	cell->alpha_eta=(4.0/15.0)*sqrt(rho*rho*T*T);
	cell->eta=cell->alpha_eta*cell->tau_eta;//// DANGER -- CHECK!!!!

	cell->tau_zeta=1.0;
	cell->alpha_zeta=1.0;
	cell->zeta=0.0;

}