#ifndef __EOS_H__
#define __EOS_H__

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
#include <unordered_map>
#include <list>
#include <Eigen/Dense>
#include "defs.h"
#include "msu_commonutils/parametermap.h"
#include "landau.h"

class CEoS{
public:
	CparameterMap *parmap;
	double cross_section;
	double kappa,mass,Kfactor,etafactor,zetafactor;
	CEoS(){
		etafactor=Kfactor=zetafactor=kappa=0.0;
		cross_section=1.0;
		mass=1.0;
	};
	void CalcEtaZetaK(CLandauCell *cell);
	CEoS(CparameterMap *parmapin);
	void CalcEtaZetaK();
	virtual void CalcEoS_of_rho_epsilon(CLandauCell *cell){ // Calculates quantities in terms of cell->epsilon and cell->rhoB
		// gives quantities in terms of epsilon and rhoB
		cell->T=cell->Pr=cell->SoverB=cell->cs2=0.0;
	}
	virtual void CalcEoS_of_rho_T(CLandauCell *cell){
		cell->epsilonk=cell->Pr=cell->SoverB=cell->cs2=cell->eta=cell->zeta=0.0;
	}
};

class CEoS_FreeGas : public CEoS{
public:
	CEoS_FreeGas(CparameterMap *parmapin);
	void CalcEoS_of_rho_epsilon(CLandauCell *cell);
	void CalcEoS_of_rho_T(CLandauCell *cell);
};

class CEoS_VdW : public CEoS{
public:
	double a;
	double rho0;
	CEoS_VdW(CparameterMap *parmapin);
	void CalcEoS_of_rho_epsilon(CLandauCell *cell);
	void CalcEoS_of_rho_T(CLandauCell *cell);
};

class CEoS_Scott : public CEoS{
public:
	double a;
	double rho0;
	CEoS_Scott(CparameterMap *parmapin);
	void CalcEoS_of_rho_epsilon(CLandauCell *cell);
	void CalcEoS_of_rho_T(CLandauCell *cell);
};

#endif
