#include "landau.h"
using namespace std;

int main(int argc,char *argv[]){
	string parsfilename;
	if(argc !=2){
		printf("Usage: hydro parsfilename\n");
		exit(1);
	}
	else
		parsfilename=argv[1];
	int nprint,iprint=0,it,ismooth,nsmooth;
	double tprint=1.0;
	CparameterMap parmap;
	parmap.ReadParsFromFile(parsfilename);
	CLandau landau(&parmap);
	nprint=lrint(tprint/landau.DELT);
	landau.currentmesh->WriteXSliceInfo(0,0);
	landau.oldmesh->CalculateBtotEtot();
	landau.currentmesh->CalculateBtotEtot();
	for(it=1;it<=landau.NT;it++){
		landau.Propagate(); // Updates newmesh using currentmesh and oldmesh
		iprint+=1;
		nsmooth=1;
		if(it==1)
			nsmooth=5;
		for(ismooth=0;ismooth<nsmooth;ismooth++){
			landau.AverageMeshes(0.5);
			landau.Propagate();
		}
		if(iprint==nprint){
			landau.newmesh->WriteXSliceInfo(0,0);
			landau.newmesh->CalculateBtotEtot();
			iprint=0;
		}
		landau.CycleMeshes();
	}
	//printf("Tlowest=%g, Thighest=%g\n",CLandauCell::Tlowest,CLandauCell::Thighest);
	//landau.currentmesh->PrintInfo();
//----------------------------------------------------------------------------------
//        CLandau rect(10.0);
//        cout << "rect area: " << rect.Propagate() << endl;
//----------------------------------------------------------------------------------
	printf("YIPPEE!!!!! I made it all the way through!\n");
	return 0;
}
