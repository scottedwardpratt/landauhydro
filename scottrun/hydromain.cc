#include "landau.h"
using namespace std;

int main(){
	int nprint,iprint=0,it;
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters.dat");
	CLandau landau(&parmap);
	nprint=landau.NT/100;
	landau.currentmesh->WriteXSliceInfo(0,0);
	for(it=1;it<=landau.NT;it++){
		landau.Propagate(); // Updates newmesh using currentmesh and oldmesh
		iprint+=1;
		if(iprint==nprint){
			landau.newmesh->WriteXSliceInfo(0,0);
			landau.newmesh->CalculateBtotEtot();
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
