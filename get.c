#include"stdio.h"
#include"stdlib.h"
int main(int argc,char **argv){
	FILE *dp;
	char filename[20];
	int i;
	double Energy[29],temp,index[29];
	sprintf(filename,"finial.dat");
	dp=fopen(filename,"r");
	
	for(i=0;i<29;i++){
		fscanf(dp,"%lg %lg %lg %lg %lg %lg\n",&index[i],&Energy[i],&temp,&temp,&temp,&temp);

	}
	fclose(dp);
	sprintf(filename,"final.dat");
	dp=fopen(filename,"w");
	for(i=0;i<29;i++){
		fprintf(dp,"%g %g\n",index[i]*29,Energy[i]+3.58657+0.009453503);

	}
	fclose(dp);
}
