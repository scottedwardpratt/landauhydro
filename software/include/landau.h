#ifndef __LANDAU_H__
#define __LANDAU_H__

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
#include <unordered_map>
#include <list>
#include <Eigen/Dense>
#include "msu_commonutils/randy.h"
#include "defs.h"
#include "msu_commonutils/parametermap.h"
class CEoS;

class CIntegralMesh;
class CHalfIntegralMesh;
class CIntegralCell;
class CHalfIntegralCell;

class CMeshParameters{
public:
	static int NX;
	static double DX0;
	static double DT;
	static double TMAX;
};

class CLandau{
public:
	CparameterMap *parmap;
	static CEoS *eos;
	string output_dirname;
	int NRungeKutta;
	
	CIntegralMesh *newIntegralMesh,*oldIntegralMesh;
	CHalfIntegralMesh *newHalfIntegralMesh,*oldHalfIntegralMesh;
	void PropagateIntegralMesh();
	void PropagateHalfIntegralMesh();
	void CycleIntegralMeshes();
	void CycleHalfIntegralMeshes();
	
	void WriteData();
	void PrintInfo();
	void WriteInfo();
	void Evolve();
	void PropagateRhoSdensPI();
	void PropagateVxKx();
	void PropagateSdens();
	void PropagatePI();
	CLandau(CparameterMap *parmapset);
	void CreateMeshes(double tset);
};

class CIntegralMesh{
public:
	static CLandau *landau;
	static CEoS *eos;
	double t;
	vector<CIntegralCell> *cell;
	CIntegralMesh();
	void Initialize(double t0);
	void UpdateQuantities();
	void CalculateBtotStot();
	void Zero();
};

class CHalfIntegralMesh{
public:
	static CLandau *landau;
	double t;
	vector<CIntegralCell> *cell;
	CHalfIntegralMesh();
	void Initialize(double test);
	void UpdateQuantities();
	void Zero();
};

class CIntegralCell{
public:
	CIntegralCell();
	static CLandau *landau;
	static CEoS *eos;
	CIntegralCell *neighborMinus,*neighborPlus;
	// These quantities refer to lower boundary
	double x,Kx_target;
	// These quantities refer to volume
	double Delx;
	double S,Q,rho,alpha_eta,alpha_zeta,alpha_gamma,tau_zeta,tau_gamma,tau_eta;
	double Pi,epsilon,epsilonk,grad2Rho;
	double T,Pr,SoverB,cs2,eta,zeta;
	vector<vector<double>> pi_shear;
	vector<vector<double>> SE;
	double pi_bulk;
	//
	void Zero();
	void CalcGrad2Rho();
	void UpdateBulkQuantities();
	void PrintInfo();
};

class CHalfIntegralCell{
public:
	static CLandau *landau;
	static CEoS *eos;
	// These quantities refer to lower boundaries
	double vx,Kx;
	// These quantities refer to volume
	vector<double> pi_shear_target;
	double pi_bulk_target;
	void GetOmega(double &omega); // omega=partial_iv_j+partial_jv_i-(2/3)del.v*delta_ij. o
	void GetDelDotV();
	void PrintInfo();
};

#endif
