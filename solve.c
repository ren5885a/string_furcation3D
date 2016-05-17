#include"string.h"
double distance(CONF *p1,CONF *p2,int NxNyNz){
int i,j,k;
double dis,distan;
int N1,N2;
distan=0;
dis=0;


	for(i=0;i<NxNyNz;i++){
		
			distan+=(p1->densityA[i]-p2->densityA[i])*(p1->densityA[i]-p2->densityA[i]);	
	
		}
		dis=sqrt(distan);	
	
return dis;
}
double splineValue (int n,double x0,
            double x[],double y[],
            double b[], double c[], double d[],double *y0
            )
{
double re,w;
int i,j,k;
int num;
for(i=1;i<n;i++){
	if((x0>=x[i-1])&&(x0<x[i])){
		num=i-1;
		break; 
	}
}
w=x0-x[num];
re=y[num] + b[num] * w + c[num] * w*w + d[num] * w*w*w;


*y0=re;
return re;

}

int chemicalPTrail(double *density,int N,double *chem){
int i,j,k;
double x,y;
x=2*density[0]-1;
y=2*density[1]-1;

	chem[0]=+4*(1-x*x-y*y)*x+(y*y*x*2)/((x*x+y*y)*(x*x+y*y));
	chem[1]=+4*(1-x*x-y*y)*y-(y*x*x*2)/((x*x+y*y)*(x*x+y*y));
//printf("x=%g y=%g chem %g %g\n",x,y,chem[0],chem[1]);

}


















