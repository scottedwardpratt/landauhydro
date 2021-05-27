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
#include "eos.h"

class CLandauMesh;
class CLandauCell;

class CLandau{
public:
	CparameterMap *parmap;
	int NX,NY,NZ,NDIM;
	int NRungeKutta;
	double DELT,DXYZ;
	CLandauMesh *currentmesh,*newmesh,*oldmesh;
	CEoS *eos;
	void CycleMeshes();
	void CreateMeshes(double t);
	
	void Propagate();
	void PropagateFirst(); // Sets up 3 meshes for first time step
	void InterpolateOldMesh(); // changes oldmesh to fit with new and current
	void PropagateRhoB();
	void PropagateT00();
	void PropagateT0i();
	void CalcEpsilonU();
	
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
	int NX,NY,NZ;
	static double DXYZ;
	double t;
	CLandauMesh(CLandau *landauset,double tset);
	void CalcU();
	void CalcEpsilon();
	void InitializeDensities();
	static CEoS *eos;
	void WriteInfo();
	void PrintInfo();
};

class CLandauCell{
public:
	CLandauCell();
	int ix,iy,iz;  // coordinates of cell
	CLandau *landau; // to 
	vector<double> u;  // velocity
	vector<double> M; // M is T_0i in restframe
	vector<vector<double>> SE;  // T_ij
	double epsilon,Pr,T;
	vector<double> jB;                
	vector<CLandauCell *> neighborPlus,neighborMinus;
	static double DXYZ;
	static int NDIM;
	static CEoS *eos;
	
	void CalcU(); // Calculates u from T0i, M and jB[0] (even though M came from estimated U)
	void CalcM(); // Calculates M from T0i, jB[0] and U
	void CalcEpsilonU(); // Calculates Epsilon and U from T00, T0i, jB[0]
	void CalcTij(); // From epsilon, jB, U
	
	double DelDotU();
	double DelDotJB();
	double DelDotT0i();
	double Grad2RhoB();
	
	void CalcGradPr(vector<double> &GradPr);
	void CalcGradRhoB(vector<double> &GradRhoB);
	void CalcDeliTij(vector<double> &DeliTij);

	void PrintInfo();
	void Zero();
};

#endif
