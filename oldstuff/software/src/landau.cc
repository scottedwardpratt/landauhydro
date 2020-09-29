#include "landau.h"

using namespace std;

CLandau::CLandau(CParameterMap *parmapset){
	parmap=parmapset;
	DXYZ=parmap->getD("LANDAU_DXYZ",0.25);
	NX=parmap->getI("LANDAU_NX",50);
	NY=parmap->getI("LANDAU_NY",1);
	NZ=parmap->getI("LANDAU_NZ",1);
	CreateMeshes(parmap->getD("LANDAU_T0",0.0));
	CLandauCell::DXYZ=DXYZ;
}

void CLandau::CreateMeshes(double tset){
	oldmesh=new CLandauMesh(this,tset-DELT);
	currentmesh=new CLandauMesh(this,tset);
	newmesh=new CLandauMesh(this,tset+DELT);
	
}

void CLandau::CycleMeshes(){
	CLandauMesh *crap;
	crap=oldmesh;
	oldmesh=currentmesh;
	currentmesh=newmesh;
	newmesh=crap;
}

void CLandau::Propagate(){
	newmesh->t=currentmesh->t+DELT;
	PropagaterhoB();
	PropagateU();
}

void CLandau::PropagaterhoB(){
	// Sukanya writes this
}

void CLandau::PropagateU(){
	// Sukanya writes this
}
