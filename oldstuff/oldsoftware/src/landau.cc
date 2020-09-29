#include "landau.h"

using namespace std;

CLandau::CLandau(CParameterMap *parmapset){
	parmap=parmapset;
	DXYZ=parmap->getD("LANDAU_DXYZ",1);
	NX=parmap->getI("LANDAU_NX",10);
	NY=parmap->getI("LANDAU_NY",1);
	NZ=parmap->getI("LANDAU_NZ",1);
	DELT=parmap->getD("LANDAU_DELT",0.1);
	CreateMeshes(parmap->getD("LANDAU_T0",0.0));
	CLandauCell::DXYZ=DXYZ;

}
//-------------------------------------------------------------
//CLandau::CLandau(float a) {
//         height=a;
//}
//-------------------------------------------------------------
void CLandau::CreateMeshes(double tset){
        
	oldmesh=new CLandauMesh(this,tset-DELT);
	currentmesh=new CLandauMesh(this,tset);
	newmesh=new CLandauMesh(this,tset+DELT);
	currentmesh->InitializeDensities();
	oldmesh->InitializeDensities();
        
}

void CLandau::CycleMeshes(){
	CLandauMesh *crap;
	crap=oldmesh;  
	oldmesh=currentmesh;
	currentmesh=newmesh;
	newmesh=crap;
}

void CLandau::Propagate(){
 	PropagaterhoB();
	PropagateU();
        
	newmesh->t=currentmesh->t+DELT;
}
//-----------------------------------------------------------------
void CLandau::PropagaterhoB(){   
	int ix,iy,iz;
	CLandauCell *c,*oldc,*newc;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(currentmesh->cell[ix][iy][iz]);
				oldc=&(oldmesh->cell[ix][iy][iz]);
				newc=&(newmesh->cell[ix][iy][iz]);
				newc->rhoB=oldc->rhoB-2.0*DELT*c->DelDotFluxB();
				if(iy==4 && iz==4)
					printf("ix,iy,iz=(%d,%d,%d): new rhoB=%g, fluxB[1]=%g,\n",ix,iy,iz,newc->rhoB,c->fluxB[1]);
			}
		}
	}        
}

void CLandau::PropagateU(){
	int ix,iy,iz,k;
	CLandauCell *c,*oldc,*newc;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				c=&(currentmesh->cell[ix][iy][iz]);
				oldc=&(oldmesh->cell[ix][iy][iz]);
				newc=&(newmesh->cell[ix][iy][iz]);
				newc->rhoB=oldc->rhoB-2.0*DELT*c->DelDotFluxB();
				for(k=1;k<4;k++){
					newc->u[k]=oldc->u[k];
					newc->fluxB[k]=newc->u[k]*newc->rhoB;
				}
				if(iy==4 && iz==4)
					printf("ix,iy,iz=(%d,%d,%d): new flux=(%g,%g,%g)\n",ix,iy,iz,newc->fluxB[1],newc->fluxB[2],newc->fluxB[3]);
			}
		}
	}        
}

