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
	CParameterMap *parmap;
	CEoS(){};
	CEoS(CParameterMap *parmapin);
	double kappa,mass;
	void eos(double epsilon,double rhoB,double &T,double &Pr);
};


#endif
