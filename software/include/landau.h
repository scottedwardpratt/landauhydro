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
	void EstimatePiS();
	void CalcBulkQuantities();
	void PropagateRho();
	void PropagateSdens();
	void PropagatePi();
	void CalcKxTarget();

	void PropagateVxKx();
	void CalcPiTarget();
	CLandau(CparameterMap *parmapset);
	void CreateMeshes(double tset);
};

class CIntegralMesh{
public:
	static CLandau *landau;
	static CEoS *eos;
	double t;
	vector<CIntegralCell *> cell;
	CIntegralMesh(double tset);
	void Initialize(double t0);
	void UpdateQuantities();
	void CalculateBtotStot();
	void Zero();
};

class CHalfIntegralMesh{
public:
	static CLandau *landau;
	double t;
	vector<CHalfIntegralCell *> cell;
	CHalfIntegralMesh(double tset);
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
	double x; // position of lower boundary
	
	// These quantities refer to volume
	double Delx;
	// These are all functions of DelX, S and Q
	double S,Q,rho,sdens;
	// These are functions of rho and sdens
	double eta,zeta,Pi,gmma;
	double alpha_eta,alpha_zeta,alpha_gmma,tau_zeta,tau_gmma,tau_eta;
	double epsilon,epsilonk,grad2Rho,alpha_Kx_factor;
	double T,Pr,SoverB,cs2;
	double Kx_target; // refers to lower boundary
	//
	// These are evolved from Eq.s of Motion
	double pi_bulk;
	vector<vector<double>> pi_shear;
	vector<vector<double>> SE;
	//
	void Zero();
	void CalcGrad2Rho();
	void Calc_Kx_target();
	void UpdateQuantities();
	void PrintInfo();
};

class CHalfIntegralCell{
public:
	CHalfIntegralCell();
	double Vx,Kx;
	static CLandau *landau;
	static CEoS *eos;
	// These quantities refer to lower boundaries
	// These quantities refer to volume
	vector<vector<double>> pi_shear_target;
	double pi_bulk_target;
	static vector<vector<double>> omega;
	void GetOmega(); // omega=partial_iv_j+partial_jv_i-(2/3)del.v*delta_ij. o
	void PrintInfo();
	void Calc_pi_shear_target();
	CHalfIntegralCell *neighborPlus,*neighborMinus;
	void Zero();
};

#endif
