#include <stdio.h>
#include <stdlib.h>
#include"math.h"
typedef struct CONF
	{
		int    	Nx,Ny,Nz;   	/* Number of lattice */
		double 	*densityA;  	/* density A*/
		double	*densityB;
		double dx,dy,dz;	//length of lattice
		double freeEnergy[0];
	}CONF; 
double distance(CONF *p1,CONF *p2,int NxNyNz);


