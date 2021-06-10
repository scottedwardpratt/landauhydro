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
	void eos(double epsilon,double rhoB,double &T,double &Pr,double &SoverB,double &cs2);
};


#endif
