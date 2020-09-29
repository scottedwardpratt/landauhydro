#include "landau.h"

void CLandauMesh::PrintInfo(){
	printf("Howdy! --- t=%g ---\n",t); 
}

CLandauMesh::CLandauMesh(CLandau *landau,double tset){      
  t=tset;
  int ix,iy,iz; 
  NX=landau->NX; NY=landau->NY; NZ=landau->NZ;
  DXYZ=landau->DXYZ;
  //------------------------------------------------------------- 
  cell.resize(NX);
  for(ix=0;ix<NX;ix++){
    cell[ix].resize(NY);
    for(iy=0;iy<NY;iy++){
      cell[ix][iy].resize(NZ);
    }
  }
  for(ix=0;ix<NX;ix++){
    for(iy=0;iy<NY;iy++){
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

void CLandauMesh::InitializeDensities(){
  int ix,iy,iz,k;
  vector<double> x(NX),y(NX),z(NX);
	CLandauCell *c;
  //-------------------------------------------------------------
  //----------ASSIGNING INITIAL CONDITION------------------------
  double x0=0.0,x1=NX*DXYZ,dg,dg2,xmid;

  for(ix=0;ix<NX;ix++){
    x[ix]=(ix+0.5)*DXYZ;
    y[ix]=(ix+0.5)*DXYZ;
    z[ix]=(ix+0.5)*DXYZ;
  }

  dg=0.1*(x1-x0);
  dg2=pow(dg,2);
  xmid= 0.5*(x0+x1);  
  //-------------------------------------------------------------
  for(ix=0;ix<NX;ix++){
    for(iy=0;iy<NY;iy++){
      for(iz=0;iz<NZ;iz++){        
	cell[ix][iy][iz].rhoB=1.0+0.1*sin(2.0*PI*x[ix]/x1)*sin(2.0*PI*y[iy]/x1)*sin(2.0*PI*z[iz]/x1);    
	printf("rhoB[%d][%d][%d]=%g\n",ix,iy,iz,cell[ix][iy][iz].rhoB);  
      } 
    }
  }
  //-------------------------------------------------------------------------
  for(ix=0;ix<NX;ix++){
    for(iy=0;iy<NY;iy++){
      for(iz=0;iz<NZ;iz++){
				c=&cell[ix][iy][iz];      
				c->u[1]=1.0;
				c->u[2]=c->u[3]=0.0;
				for(k=1;k<4;k++)
					c->fluxB[k]=c->rhoB*c->u[k];
				printf("flux_x(%d,%d,%d)=%g\n",ix,iy,iz,cell[ix][iy][iz].fluxB[1]);                            
      } 
    }
  }   
  //---------------------------------------------------------------------------
}


