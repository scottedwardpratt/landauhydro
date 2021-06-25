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
#include "randy.h"
#include "defs.h"
#include "parametermap.h"
class CEoS;

class CLandauMesh;
class CLandauCell;

class CLandau{
public:
	CparameterMap *parmap;
	int NX,NY,NZ,NDIM,NT;
	int NRungeKutta;
	double DELT,DXYZ,TMAX;
	CLandauMesh *currentmesh,*newmesh,*oldmesh;
	CEoS *eos;
	void CycleMeshes();
	void CreateMeshes(double t);
	
	void Propagate();
	void InterpolateOldMesh(); // changes oldmesh to fit with new and current
	void PropagateRhoBPdens();
	void AverageMeshes(double weight);
	
	void WriteData1D();
	void PrintInfo();
	void WriteInfo();
	
	void solve(double *a, double *b, double *c, double *d, int n);
	void ClearMesh();
	CLandau(CparameterMap *parmapset);
	vector<vector<vector<CLandauCell>>> cell;
	double ap[100],bp[100],cp[100],dp[100],AA[100],BB[100],rb[100];
};

class CLandauMesh{
public:
	vector<vector<vector<CLandauCell>>> cell;
	static CLandau *landau;
	static int NX,NY,NZ,NDIM;
	static double DXYZ;
	double t;
	CLandauMesh(CLandau *landauset,double tset);
	void Initialize(double t);
	static CEoS *eos;
	void WriteInfo();
	void WriteXSliceInfo(int iy,int iz);
	void PrintInfo();
	void CalculateUJMEpsilonSE();
	void CalculateBtotEtot();
};

class CLandauCell{
public:
	static double Tlowest,Thighest;
	CLandauCell();
	int ix,iy,iz;  // coordinates of cell
	vector<double> u;  // velocity
	vector<double> M; // M is T_0i in restframe (only due to kappa)
	vector<double> Pdens; // Momentum density (not same as T0i in non-rel theory)
	vector<vector<double>> SE;  // in lab frame (includes KE...)
	double epsilonk,Pr,T,SoverB,cs2,K;
	vector<double> jB;
	vector<double>Kflow;
	vector<CLandauCell *> neighborPlus,neighborMinus;
	static double DXYZ;
	static int NDIM;
	static CEoS *eos;
	
	void CalcM(); // Calculates M from T0i, jB[0] and U
	void CalcEpsilonSE(); // Calculates Epsilon and SE tensor from Pdens, U, M
	
	double DelDotU();
	double DelDotJB();
	double DelDotT0i();
	double Grad2RhoB();
	
	void CalcGradPr(vector<double> &GradPr);
	void CalcGradRhoB(vector<double> &GradRhoB);
	void CalcDeliTij(vector<double> &DeliTij);
	double CalcDivKFlow();
	void CalcKFlow();

	void PrintInfo();
	void Zero();
};

#endif
