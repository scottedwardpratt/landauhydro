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

class CLandauMesh;
class CLandauCell;

class CLandau{
public:
	CParameterMap *parmap;
	int NX,NY,NZ;
	double DELT,DXYZ;
	CLandauMesh *currentmesh,*newmesh,*oldmesh;
	CLandauCell *neighborxl,*neighborxr,*neighboryl,*neighboryr,*neighborzl,*neighborzr;
	void CycleMeshes();
	void CreateMeshes(double t);
	void Propagate();
	void PropagaterhoB();
	void PropagateU();
	void ClearMesh();
	CLandau(CParameterMap *parmapset);
	vector<vector<vector<CLandauCell>>> cell; 
	double epsilon,P,rhoB;
       
};

class CLandauMesh{
public:
	vector<vector<vector<CLandauCell>>> cell;
	static CLandau *landau;
	int NX,NY,NZ;
	double DXYZ;
	double t;
	CLandauMesh(CLandau *landauset,double tset);
	void PrintInfo();
	void InitializeDensities();
};

class CLandauCell{
public:
	vector<vector<vector<CLandauCell>>> cell;
	int i,j,k;  // coordinates of cell
	CLandau *landau; // to 
	vector<double> u;
	vector<double> fluxB;
	vector<vector<double>> SE;
	double epsilon,P,rhoB,fluxd;
	CLandauCell *neighborxl,*neighborxr,*neighboryl,*neighboryr,*neighborzl,*neighborzr;
	static double DXYZ;
	double DrhoBDX();
	double D2rhoBDX2();
	double DrhoBDY();
	double D2rhoBDY2();
	double DrhoBDZ();
	double D2rhoBDZ2();
	double dUxdX();
	double dUxdY();
	double dUxdZ();
	double dUydX();
	double dUydY();
	double dUydZ();
	double dUzdX();
	double dUzdY();
	double dUzdZ();
	double DelDotU(){
		return dUxdX()+dUydY()+dUzdZ();
	}
	double DelDotFluxB(){
		return rhoB*DelDotU()+u[1]*DrhoBDX()+u[2]*DrhoBDY()+u[3]*DrhoBDZ();
	}
	void PrintInfo();
	void Zero();    
	CLandauCell(){
		u.resize(4);
		fluxB.resize(4);
		SE.resize(4);
		for(int i=0;i<4;i++)
			SE[i].resize(4);
	}
};

#endif
