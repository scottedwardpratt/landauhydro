#include <cstdio>
#include <cmath>
#include <iostream>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include<Eigen/SparseLU>
#include<Eigen/SparseQR>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>


const double PI=4.0*atan(1.0);
const int D=40;  // dimension
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> CTriplet;
void GetCoefficients(SpMat& A, Eigen::VectorXd& F);

void GetXYZ(int i,int &X,int &Y,int &Z){
	X=floorl(i/(D*D));
	Y=floorl((i-D*D*X)/D);
	Z=floorl(i-D*D*X-D*Y);
}
int GetI(int X,int Y,int Z){
	return X*D*D+Y*D+Z;
}

int main(int argc, char** argv){
	int i,j,X,Y,Z,M = D*D*D;  // dimension
	vector<CTriplet> coefficients;            // list of non-zero coefficients
	Eigen::VectorXd b(M),bb(M),phi(M);        // the right hand side-vector resulting from the constraints	
	SpMat A(M,M);
	
	/* ... fill A and b ... */
	GetCoefficients(A,b);
	Eigen::BiCGSTAB<SpMat> solver;
	solver.compute(A);
	phi = solver.solve(b);
	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error()      << std::endl;
	/* ... update b ... */
	phi = solver.solve(b); // solve again
	
	double m2=1.0,k2=12.0*PI*PI/(D*D);
	double E2=(k2+m2);
	double ratio;
	bb=A*phi;
	for(i=0;i<M;i++){
		GetXYZ(i,X,Y,Z);
		ratio=E2*phi(i)/b(i);
		printf("%7d (%3d,%3d,%3d): phi=%10.6f: f=%10.6f =? %10.6f, ratio=%g\n",i,X,Y,Z,phi(i)*E2,b(i),bb(i),ratio);
	}
	
	return 0;
}

void GetCoefficients(SpMat& A, Eigen::VectorXd& b){
	int i,j,x,y,z,xx,yy,zz,M=D*D*D;
	double m2=1.0;
	double k2=12.0*PI*PI/(D*D);
	for(x=0;x<D;x++){
		for(y=0;y<D;y++){
			for(z=0;z<D;z++){
				
				i=GetI(x,y,z);

				b(i)=(k2+m2)*sin(2.0*PI*double(x)/double(D));
				b(i)*=sin(2.0*PI*double(y)/double(D));
				b(i)*=sin(2.0*PI*double(z)/double(D));
				
				A.insert(i,i)=6.0+m2;
				
				xx=x+1;
				if(xx==D)
					xx=0;
				j=GetI(xx,y,z);
				A.insert(i,j)=-1.0;
				
				xx=x-1;
				if(xx==-1)
					xx=D-1;
				j=GetI(xx,y,z);
				A.insert(i,j)=-1.0;
				
				yy=y+1;
				if(yy==D)
					yy=0;
				j=GetI(x,yy,z);
				A.insert(i,j)=-1.0;
				
				yy=y-1;
				if(yy==-1)
					yy=D-1;
				j=GetI(x,yy,z);
				A.insert(i,j)=-1.0;
				
				zz=z+1;
				if(zz==D)
					zz=0;
				j=GetI(x,y,zz);
				A.insert(i,j)=-1.0;
				
				zz=z-1;
				if(zz==-1)
					zz=D-1;
				j=GetI(x,y,zz);
				A.insert(i,j)=-1.0;
				
			}
		}
	}
}


/*
void GetCoefficients(SpMat& A, Eigen::VectorXd& F){
	int i,j,M=D*D*D;
	double element;
	for(i=0;i<M;i++){
		F(i)=sqrt(double(i*(M-i+1)))+1.0;
		element=double(i)+1.0+double(i*(M-i));
		A.insert(i,i)=element;
		for(j=0;j<M;j+=4){
			if(i!=j){
				element=i+1.0+double(j*(M-j+1)+i*(M-i+1))-i*(M-j+1)-j*(M-i+1);
				A.insert(i,j)=element;
			}
		}
	}	
}
*/
