#include "landau.h"
using namespace std;

int main(int argc, char *argv[]){
	int nprint=1000,iprint=0;
	CParameterMap parmap;
	parmap.ReadParsFromFile("parameters.dat");
	CLandau landau(&parmap);
	
	int it,nt;
	
	printf("Enter number of time steps: (delt=%g): ",landau.DELT);
	scanf("%d",&nt);
	landau.PropagateFirst();
	landau.WriteInfo();
	for(it=1;it<=nt;it++){
		landau.Propagate(); // Updates newmesh using currentmesh and oldmesh
		iprint+=1;
		if(iprint==nprint){
			landau.WriteInfo();
			iprint=0;
		}
		landau.CycleMeshes();
	}
	//landau.currentmesh->PrintInfo();
//----------------------------------------------------------------------------------
//        CLandau rect(10.0);
//        cout << "rect area: " << rect.Propagate() << endl;
//----------------------------------------------------------------------------------
	printf("YIPPEE!!!!! I made it all the way through!\n");
	return 0;
}
