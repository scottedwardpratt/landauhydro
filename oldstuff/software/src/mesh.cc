#include "landau.h"

void CLandauMesh::PrintInfo(){
	printf("Howdy! --- t=%g ---\n",t);
}

CLandauMesh::CLandauMesh(CLandau *landau,double tset){
	t=tset;
	NX=landau->NX; NY=landau->NY; NZ=landau->NZ;
	DXYZ=landau->DXYZ;
	int ix,iy,iz;
	cell.resize(NX);
	for(ix=0;ix<NX;ix++){
		cell[ix].resize(NY);
		for(iy=0;iy<NY;iy++){
			cell[ix][iy].resize(NZ);
			for(iz=0;iz<NZ;iz++){
				if(ix>0)
					cell[ix][iy][iz].neighborxl=&cell[ix-1][iy][iz];
				else
					cell[ix][iy][iz].neighborxl=&cell[NX-1][iy][iz];
				if(ix<NX-1)
					cell[ix][iy][iz].neighborxr=&cell[ix+1][iy][iz];
				else
					cell[ix][iy][iz].neighborxr=&cell[0][iy][iz];
				if(iy>0)
					cell[ix][iy][iz].neighboryl=&cell[ix][iy-1][iz];
				else
					cell[ix][iy][iz].neighboryl=&cell[ix][NY-1][iz];
				if(iy<NY-1)
					cell[ix][iy][iz].neighboryr=&cell[ix][iy+1][iz];
				else
					cell[ix][iy][iz].neighboryr=&cell[ix][0][iz];
				if(iz>0)
					cell[ix][iy][iz].neighborzl=&cell[ix][iy][iz-1];
				else
					cell[ix][iy][iz].neighborzl=&cell[ix][iy][NZ-1];
				if(iz<NZ-1)
					cell[ix][iy][iz].neighborzr=&cell[ix][iy][iz+1];
				else
					cell[ix][iy][iz].neighborzr=&cell[ix][iy][0];
				cell[ix][iy][iz].Zero();
			}
		}
	}
}


