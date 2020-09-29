#include "landau.h"
using namespace std;

int main(int argc, char *argv[]){
	CParameterMap parmap;
	parmap.ReadParsFromFile("parameters.dat");
	CLandau landau(&parmap);
	
	int it,nt;
	
	printf("Enter number of time steps: (delt=%g): ",landau.DELT);
	scanf("%d",&nt);
	for(it=1;it<=nt;it++){
		landau.Propagate();  // Updates newmesh using currentmesh and oldmesh
		landau.CycleMeshes();
	}
	landau.currentmesh->PrintInfo();
	return 0;
}
