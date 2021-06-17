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
#include "parametermap.h"

class CEoS{
public:
	CparameterMap *parmap;
	CEoS(){};
	CEoS(CparameterMap *parmapin);
	double kappa,mass;
	virtual void eos(double epsilon,double rhoB,double &T,double &Pr,double &SoverB,double &cs2){
		// gives quantities in terms of epsilon and rhoB
		T=Pr=SoverB=cs2=0.0;
	};
};

class CEoS_FreeGas : public CEoS{
public:
	CEoS_FreeGas(CparameterMap *parmapin);
	void eos(double epsilon,double rhoB,double &T,double &Pr,double &SoverB,double &cs2);
};

class CEoS_VdW : public CEoS{
public:
	double a;
	double rho0;
	CEoS_VdW(CparameterMap *parmapin);
	void eos(double epsilon,double rhoB,double &T,double &Pr,double &SoverB,double &cs2);
};

#endif
