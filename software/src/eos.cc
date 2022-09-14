#include "eos.h"

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

CEoS_Scott::CEoS_Scott(CparameterMap *parmapset){
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
	
	//sigma_K=sqrt(35.0*rhoB*T*T*T/(4.0*mass));
	sigma_K=1.0/rhoB;
	tau_K=sqrt(mass/T)/rhoB;
	K=(T/mass)*tau_K; // usual conductivity
	K=K/(sigma_K*rhoB);   // scaled conductivity
	
	sigma_eta=sqrt(Etafactor*rhoB*T*T);
	tau_eta=sqrt(mass/T)/rhoB;
	eta=(sigma_eta*sigma_eta/T)*tau_eta; // usual viscosity
	eta=eta/(sigma_eta*rhoB);  // scaled viscosity
	
	sigma_zeta=1.0;
	tau_zeta=0.0;
	zeta=0.0;
	
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
	double T,Pr,SoverB,cs2,K,tau_K,sigma_K,eta,tau_eta,sigma_eta,zeta,tau_zeta,sigma_zeta,cs2factor;
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
	
	//sigma_K=sqrt(35.0*rhoB*T*T*T/(4.0*mass));
	cs2factor=0.1*sqrt(35.0*rhoB*T*T*T/(4.0*mass));
	if(cs2>0.0){
		cs2factor=cs2factor+0.9*cs2;
	}
	
	CalcEtaZetaK(cell);
	
	cell->T=T;
	cell->Pr=Pr;
	cell->SoverB=SoverB;
	cell->cs2=cs2;

}

void CEoS_VdW::CalcEoS_of_rho_T(CLandauCell *cell){
	double rhoB=cell->jB[0];
	cell->epsilonk=1.5*rhoB*cell->T-a*rhoB*rhoB;
	CalcEoS_of_rho_epsilon(cell);
}

void CEoS_Scott::CalcEoS_of_rho_epsilon(CLandauCell *cell){
	double epsilonk=cell->epsilonk,rhoB=cell->jB[0];
	double T,Pr,SoverB,cs2,K,tau_K,sigma_K,eta,tau_eta,sigma_eta,zeta,tau_zeta,sigma_zeta,cs2factor;
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
	
	//sigma_K=sqrt(35.0*rhoB*T*T*T/(4.0*mass));
	cs2factor=0.1*sqrt(35.0*rhoB*T*T*T/(4.0*mass));
	if(cs2>0.0){
		cs2factor=cs2factor+0.9*cs2;
	}
	
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
	tau_K=sqrt(mass/T)/rhoB;
	K=Kfactor*(T/mass)*tau_K;
	sigma_K=1.0;  //sqrt(21.0*rhoB*T*T*cs2factor/4.0); // scaling factor
	K_scaled=K/sigma_K;   // scaled condictivity
	
	tau_eta=sqrt(mass/T)/rho_B;
	eta=(rhoB*T)*tau_eta; // usual viscosity
	sigma_eta=sqrt(rhoB*T)
	eta_scaled=Etafactor*eta/sigma_eta;  // scaled viscosity
	
	sigma_zeta=1.0;
	tau_zeta=0.0;
	zeta=0.0;
	zeta_scaled=0.0;
}