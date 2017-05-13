
//////HW10.cpp		 
//////Author:	YUNG-SHAN SU 
//////Date: 5/13 




#include<stdio.h>
#include<stdlib.h>
#include<string.h> 
#include"MAT.h"
#include<math.h>
#include<sys/time.h>
void NewtonMethod(VEC& x, double errorMax);
VEC function1(VEC x,double vdd);
MAT Jacobi1(VEC x,double vdd);
int main(int argc,char **argv){
	double errorMax=0.000000001;
	int dim=2;
	double vdd=1;
	VEC x(2);
	NewtonMethod( x, errorMax,vdd);
	x.print();
	return 0;
}	
VEC function1(VEC x,double v){
	double r2=0.01;
	double Is=1;
	double fi=0.026;
	VEC Fx(x.len());
	if (v>0){
		Fx[0]=(x[0]-x[1])/r2-Is*(exp((v-x[0])/fi)-1);
		Fx[1]=(x[0]-x[1])/r2-Is*(exp(x[1]/fi)-1);
	}
	else if (v<0){
		Fx[0]=(x[0]-x[1])/r2-Is*(exp((-x[0])/fi)-1);
		Fx[1]=(x[0]-x[1])/r2-Is*(exp((x[1]-v)/fi)-1);
	}
	else{
		Fx[0]=0;
		Fx[1]=0;
	}
	return Fx;
}
MAT Jacobi1(VEC x,double vdd){
	double r2=0.01;
	double Is=1;
	double fi=0.026;
	MAT Jf(x.len());
	if (v>0){
		Jf[0][0]=1.0/r2 + 1.0/fi*Is*(exp((vdd-x[0])/fi));
		Jf[0][1]=-1.0/r2;
		Jf[1][0]=1.0/r2 ;
		Jf[1][1]=-1.2/r2-1.0/fi*Is*exp(x[1]/fi);
	}
	else if (v<0){
		Jf[0][0]=1.0/r2 + 1.0/fi*Is*(exp((-x[0])/fi));
		Jf[0][1]=-1.0/r2;
		Jf[1][0]=1.0/r2 ;
		Jf[1][1]=-1.0/r2-1.0/fi*Is*exp((x[1]-vdd)/fi);
	}
	else{
		Jf=0;
	}
}

void NewtonMethod(VEC& x, double errorMax,double vdd){
	int k=0;
	error=errorMax+1;
	VEC Fx(x.len());
	MAX Jf(x.len());
	VEC deltaX(x.dim());
	while(error>=errorMax){
		Fx=function1(x,vdd);
		Jf=Jacobi1(x,vdd);
		getLinearSolution(Jf,deltaX,(Fx*-1));
		x=x+deltaX;
		error=max_norm(Fx(x));
	}
			
}
