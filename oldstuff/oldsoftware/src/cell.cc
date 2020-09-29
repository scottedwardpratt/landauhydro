#include "landau.h"
double CLandauCell::DXYZ=0.0;

void CLandauCell::PrintInfo(){
	printf("-----CellInfo: i=%d -----\n",i);
	printf("rhoB=%g, epsilon=%g, u=(%g,%g,%g)\n",rhoB,epsilon,u[1],u[2],u[3]);
}
                        
double CLandauCell::DrhoBDX(){
	return 0.5*(neighborxr->rhoB-neighborxl->rhoB)/DXYZ;                       
}

double CLandauCell::D2rhoBDX2(){
	return (neighborxr->rhoB-2.0*rhoB+neighborxl->rhoB)/(DXYZ*DXYZ);
}
double CLandauCell::DrhoBDY(){
	return 0.5*(neighboryr->rhoB-neighboryl->rhoB)/DXYZ;
}
double CLandauCell::D2rhoBDY2(){
	return (neighboryr->rhoB-2.0*rhoB+neighboryl->rhoB)/(DXYZ*DXYZ);
}
double CLandauCell::DrhoBDZ(){
	return 0.5*(neighborzr->rhoB-neighborzl->rhoB)/DXYZ;
}
double CLandauCell::D2rhoBDZ2(){
	return (neighborzr->rhoB-2.0*rhoB+neighborzl->rhoB)/(DXYZ*DXYZ);
}
double CLandauCell::dUxdX(){
	return 0.5*(neighborxr->u[1]-neighborxl->u[1])/DXYZ;
}
double CLandauCell::dUxdY(){
	return 0.5*(neighboryr->u[1]-neighboryl->u[1])/DXYZ;
}
double CLandauCell::dUxdZ(){
	return 0.5*(neighborzr->u[1]-neighborzl->u[1])/DXYZ;
}
double CLandauCell::dUydX(){
	return 0.5*(neighborxr->u[2]-neighborxl->u[2])/DXYZ;
}
double CLandauCell::dUydY(){
	return 0.5*(neighboryr->u[2]-neighboryl->u[2])/DXYZ;
}
double CLandauCell::dUydZ(){
	return 0.5*(neighborzr->u[2]-neighborzl->u[2])/DXYZ;
}
double CLandauCell::dUzdX(){
	return 0.5*(neighborxr->u[3]-neighborxl->u[3])/DXYZ;
}
double CLandauCell::dUzdY(){
	return 0.5*(neighboryr->u[3]-neighboryl->u[3])/DXYZ;
}
double CLandauCell::dUzdZ(){
	return 0.5*(neighborzr->u[3]-neighborzl->u[3])/DXYZ;
}
void CLandauCell::Zero(){
	u[0]=1.0; u[1]=u[2]=u[3]=0.0;
	fluxB[0]=fluxB[1]=fluxB[2]=fluxB[3]=0.0;
	rhoB=epsilon=0.0;
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			SE[i][j]=0.0;
}
