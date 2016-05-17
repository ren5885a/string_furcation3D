#include"string.h"
#include"spline.h"
#include<mpi.h>
#include<time.h>
//#include "mpi.h"
#include</usr/include/fftw3.h>
#define MaxIT 100           //Maximum iteration steps
#define Nx 6 
#define Ny 96
#define Nz 96			//grid size

#define MaxLoop 500

long NxNyNz, NxNyNz1;
int Nxh1;

/*#define NxNyNz (Nx*Ny*Nz)

#define Nxh1 (Nx/2+1)
#define NxNyNz1 (Nz*Ny*Nxh1)*/

#define Pi 3.141592653589


double freeE(double *wA,double *wB,double *phA,double *phB);
double getConc(double *phlA,double *phlB,double *wA,
		double *wB);
void sovDifFft(double *g,double *w,double *qInt,
	        double z,int ns,int sign);
void write_ph(double *phA,double *phB,double *wA,double *wB);
double freeE_constrained(double *wA,double *wB,double *phA,double *phB,double *PhiA,double lambda);
void translation_good(double *density,double *target);
int translation(double *density,double *density_trans,int direction,int dis);
double distance_field(double *phi,double *phi_target);

int NsA, dNsB, Narm;
double kx[Nx],kz[Nz],ky[Ny],*kxyzdz,dx,dy,dz,*wdz;
double lx, ly, lz;
double hAB, fA, fB, dfB, ds0, ds2;
double *in;
fftw_complex *out;
fftw_plan p_forward, p_backward;
double temp;
char FEname[50], phname[50];
double splineValue (int n,double x0,
            double x[],double y[],
            double b[], double c[], double d[],double *y0
            );

int main(int argc,char **argv){
	MPI_Status status;
	MPI_Init(&argc, &argv);
	int my_node,totalnodes;
	double *wA_1,*wB_1,*phA_1,*phB_1;
	double e1,e2,e3,e4,ksq;
	double rjk,yj,zk,phat,phbt;
	int i,j,k,intag,iseed=-3,ntyp; //local_x_starti;
	long ijk;
	char comment[201];
	char density_name[20];
	CONF *p,*pj;
	double *PhiA;
	int loop;
	FILE *fp,*dp_out;
	time_t ts;
	double dis;
	double *X,*X_Inv;
	double *Y;
	double x0;
	double y0;
	double x[100],y[100],d[100],b[100],c[100];
	int iflag;
	double lambda;
	double *chem_1,*chem_2;
	double delta;
	double FreeEnergy_derive;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_node);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);	//MPI_Status status;
	
	Nxh1=Nx/2+1;
	NxNyNz=Nx*Ny*Nz;
	NxNyNz1=Nx*Ny*Nxh1;
	iseed=time(&ts);
	srand48(iseed);
	
	fp=fopen("ab.txt","r");
	fscanf(fp,"%d",&intag);		//in=1: inputing configuration is given;
	fscanf(fp,"%lf",&hAB);
	fscanf(fp, "%lf", &fA);
	fscanf(fp,"%lf, %lf, %lf",&lx, &ly, &lz);
	fscanf(fp,"%s",FEname);				//output file name for parameters;
	fscanf(fp,"%s",phname);				//output file name for configuration;
	fscanf(fp, "%lf", &ds0);
	fscanf(fp, "%d", &Narm);
	fclose(fp);

	ds2=ds0/2;

	fB=1.0-fA;

	dx=lx/Nx;
	dy=ly/Ny;
	dz=lz/Nz;
	dfB=fB/Narm;
	NsA = ((int)(fA/ds0+1.0e-8));
	dNsB = ((int)(dfB/ds0+1.0e-8));
	if(my_node!=0){
		

		wA_1=(double *)malloc(sizeof(double)*NxNyNz);
		wB_1=(double *)malloc(sizeof(double)*NxNyNz);
		phA_1=(double *)malloc(sizeof(double)*NxNyNz);
		phB_1=(double *)malloc(sizeof(double)*NxNyNz);
		
		kxyzdz=(double *)malloc(sizeof(double)*NxNyNz);
		wdz=(double *)malloc(sizeof(double)*NxNyNz);
		
		in=(double *)malloc(sizeof(double)*NxNyNz); /* for fftw3 */
		out=(fftw_complex *)malloc(sizeof(fftw_complex)*NxNyNz);

        	p_forward = fftw_plan_dft_r2c_3d(Nz, Ny, Nx, in, out, FFTW_ESTIMATE);
        	p_backward = fftw_plan_dft_c2r_3d(Nz, Ny, Nx, out, in, FFTW_ESTIMATE);


		
		for(i=0;i<=Nx/2-1;i++)kx[i]=2*Pi*i*1.0/Nx/dx;
		for(i=Nx/2;i<Nx;i++)kx[i]=2*Pi*(i-Nx)*1.0/dx/Nx;
		for(i=0;i<Nx;i++)kx[i]*=kx[i];

		for(i=0;i<=Ny/2-1;i++)ky[i]=2*Pi*i*1.0/Ny/dy;
		for(i=Ny/2;i<Ny;i++)ky[i]=2*Pi*(i-Ny)*1.0/dy/Ny;
		for(i=0;i<Ny;i++)ky[i]*=ky[i];

    		for(i=0;i<=Nz/2-1;i++)kz[i]=2*Pi*i*1.0/Nz/dz;
    		for(i=Nz/2;i<Nz;i++)kz[i]=2*Pi*(i-Nz)*1.0/dz/Nz;
    		for(i=0;i<Nz;i++)kz[i]*=kz[i];	

		for(k=0;k<Nz;k++)for(j=0;j<Ny;j++)for(i=0;i<Nx;i++)
		{
			ijk=(long)((k*Ny+j)*Nx+i);
			ksq=kx[i]+ky[j]+kz[k];
			kxyzdz[ijk]=exp(-ds0*ksq);
		}
		
		if(intag==1024){
			sprintf(density_name,"density_%d.dat",my_node);
			fp=fopen(density_name,"r");
			fgets(comment,200,fp);       
			fgets(comment,200,fp);
			PhiA=(double *)malloc(sizeof(double)*NxNyNz);
			
               		 for(ijk=0;ijk<NxNyNz;ijk++){
				fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&temp,&temp,&temp,&PhiA[ijk],&temp,&temp,&temp,&temp,&temp,&temp,&temp);
				
			}
                	fclose(fp);
		}
		else if(intag==1026){
			sprintf(density_name,"phi_%d.dat",my_node);
			fp=fopen(density_name,"r");
			fgets(comment,200,fp);       
			fgets(comment,200,fp);
			PhiA=(double *)malloc(sizeof(double)*NxNyNz);
			
               		 for(ijk=0;ijk<NxNyNz;ijk++){
				fscanf(fp,"%lg %lg %lg %lg\n",&PhiA[ijk],&temp,&wA_1[ijk],&wB_1[ijk]);
				
			}
                	fclose(fp);
		}
		else if(intag==1028){
			sprintf(density_name,"pha_%d.dat",my_node);
			fp=fopen(density_name,"r");
			fgets(comment,200,fp);       
			fgets(comment,200,fp);
			PhiA=(double *)malloc(sizeof(double)*NxNyNz);
			
               		 for(ijk=0;ijk<NxNyNz/6;ijk++){
				fscanf(fp,"%lf %lf %lf %lf\n",&PhiA[ijk*6],&temp,&wA_1[ijk*6],&wB_1[ijk*6]);
				for(i=1;i<6;i++){
					wA_1[ijk*6+i]=wA_1[ijk*6];
					wB_1[ijk*6+i]=wB_1[ijk*6];
					PhiA[ijk*6+i]=PhiA[ijk*6];
				}
			}
                	fclose(fp);
		}
		MPI_Send(PhiA,NxNyNz,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
		
		//MPI_Send(PhiB,NxNyNz,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
	}//end of node !=0
	else if(my_node==0){

		p=(CONF *)malloc((totalnodes-1)*sizeof(CONF));
		
		
		for(j=1;j<totalnodes;j=j+1){
				pj=&p[j-1];
				pj->densityA=(double *)malloc(NxNyNz*sizeof(double));
				//pj->densityB=(double *)malloc(NxNyNz*sizeof(double));
		}
		X=(double *)malloc((totalnodes-1)*sizeof(double));
		X_Inv=(double *)malloc((totalnodes-1)*sizeof(double));
		Y=(double *)malloc((totalnodes-1)*sizeof(double));
	}


	for(loop=0;loop<=MaxLoop;loop++){
		if(my_node==0){
			//printf("loop=%d\n",loop);
			
			for(j=1;j<totalnodes;j=j+1){
				pj=&p[j-1];
				
				MPI_Recv(pj->densityA,NxNyNz,MPI_DOUBLE,j,1,MPI_COMM_WORLD, &status);
				
		
			}
			//translation_good(pj->densityA,p[0].densityA);	
			dis=0;
			for(j=0;j<=totalnodes-3;j++){
				X_Inv[j+1]=distance(&p[j],&p[j+1],NxNyNz);
				dis+=X_Inv[j+1];
				X[j+1]=dis;
		
			}
			X[0]=0;
			printf("%d ",loop);
			for(j=0;j<totalnodes-1;j++){
				pj=&p[j-1];
				printf("%g ",X[j]);
			}
			printf("\n");
			for(j=1;j<=totalnodes-2;j++){
				X[j]=X[j]/dis;
			}
			
			for(i=0;i<NxNyNz;i++){
				for(j=0;j<=totalnodes-2;j++){
					pj=&p[j];
					Y[j]=pj->densityA[i];
					
				}
			
			
				
				spline (totalnodes-1,0, 0,0.0, 0.0,X,Y,b,c,d,&iflag);
				
				
				for(j=1;j<=totalnodes-3;j++){
					pj=&p[j];
					x0=(double)j/(double)(totalnodes-2);
					splineValue(totalnodes-1,x0,X,Y,b,c,d,&y0);
					pj->densityA[i]=y0;
					
				}
				
				

			}
			
			for(j=1;j<totalnodes;j=j+1){
				pj=&p[j-1];
				
				MPI_Send(pj->densityA,NxNyNz,MPI_DOUBLE,j,1,MPI_COMM_WORLD);
				
		
			}
			
			
		}// end if my_node==0
		if(my_node!=0){
			MPI_Recv(PhiA,NxNyNz,MPI_DOUBLE,0,1,MPI_COMM_WORLD, &status);
			//printf("my_node=%d %g\n",my_node,PhiA[1]);
			
			lambda=40;
			e1=freeE_constrained(wA_1,wB_1,phA_1,phB_1,PhiA,lambda);
			
			
			if(loop%5==0){
			
				
				sprintf(phname,"freeEnergy_%d.dat",loop);
				

				
				dp_out=fopen(phname,"a");
				fprintf(dp_out,"%d %g\n",my_node,e1);
				fclose(dp_out);
				sprintf(phname,"pha_%d.dat",my_node);
				write_ph(phA_1,phB_1,wA_1,wB_1);
			
			}
			
			for(ijk=0;ijk<NxNyNz;ijk++){
				PhiA[ijk]=phA_1[ijk];
				
			}
				
			
			
			
			if(loop<MaxLoop){
				MPI_Send(PhiA,NxNyNz,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
	
			}
		}// end if my_node!=0
		

		if(loop==MaxLoop){
			if(my_node==0){
				
			}
			else if(my_node!=0){
				sprintf(phname,"pha_%d.dat",my_node);
				write_ph(phA_1,phB_1,wA_1,wB_1);
				
			}
			break;
		}

	}//end loop


 	fftw_destroy_plan ( p_forward );
        fftw_destroy_plan ( p_backward );
	if(my_node!=0){
		free(wA_1);
		free(wB_1);
	
		free(phA_1);
		free(phB_1);
		
		free(kxyzdz);
		free(wdz);

		free(in);
		free(out);
	}
	MPI_Finalize(); 
	return 0;
}






double freeE_constrained(double *wA,double *wB,double *phA,double *phB,double *PhiA,double lambda)
{
	int i,j,k,iter,maxIter;
	long ijk;
	double freeEnergy,freeOld,qCab,eta;
	double freeW,freeAB,freeS,freeDiff,freeWsurf;
	double Sm1,Sm2,wopt,wcmp,beta,psum,fpsum;
	double waDiff,wbDiff,inCompMax,wa0, wb0, wc0;
	FILE *fp;
	
	
	Sm1=0.2e-7;
	Sm2=0.1e-10;
	maxIter=MaxIT;
	wopt=0.010;
	wcmp=0.10;
	beta=1.0;
	iter=0;	

	freeEnergy=0.0;
	
	do
	{
		iter=iter+1;

		wa0 = 0.0;
		wb0 = 0.0;

		for(ijk=0; ijk<NxNyNz;ijk++)
		{
			wa0 += wA[ijk];
			wb0 += wB[ijk];
		}

		wa0/=NxNyNz;
		wb0/=NxNyNz;
		
		for(ijk=0; ijk<NxNyNz;ijk++)
		{
			wA[ijk]-=wa0;
			wB[ijk]-=wb0;
		}
		
		qCab=getConc(phA,phB,wA,wB);
		
		freeW=0.0;
		freeAB=0.0;
		freeS=0.0;
		freeWsurf=0.0;
		inCompMax=0.0;
				
		for(ijk=0; ijk<NxNyNz; ijk++)
		{
			eta=(wA[ijk]+wB[ijk]-hAB)/2;

			psum=1.0-phA[ijk]-phB[ijk];
			fpsum=fabs(psum);
			if(fpsum>inCompMax)inCompMax=fpsum;
			waDiff=hAB*phB[ijk]+eta-lambda*(PhiA[ijk]-phA[ijk])-wA[ijk];
			wbDiff=hAB*phA[ijk]+eta-wB[ijk]-lambda*(1-PhiA[ijk]-phB[ijk]);
			waDiff-=wcmp*psum;
			wbDiff-=wcmp*psum;
			
			freeAB=freeAB+hAB*phA[ijk]*phB[ijk];
			freeW=freeW-(wA[ijk]*phA[ijk]+wB[ijk]*phB[ijk]);
			
			wA[ijk]+=wopt*waDiff;
			wB[ijk]+=wopt*wbDiff;
		}											
			
		freeAB/=NxNyNz;
		freeW/=NxNyNz;
		freeWsurf/=NxNyNz;
		freeS=-log(qCab);
			
		freeOld=freeEnergy;
		freeEnergy=freeAB+freeW+freeS;

		//**** print out the free energy and error results ****
			
		if(iter%5==0||iter>=maxIter)
        	{
			//fp=fopen("printout.txt","a");
			//fprintf(fp,"%d\n",iter);
			//fprintf(fp,"%10.8e, %10.8e, %10.8e, %10.8e, %10.8e, %e\n",
				//freeEnergy,freeAB,freeW,freeS,freeWsurf,inCompMax);
			//fclose(fp);
		}
		if(iter%500==0)//printf(" %5d : %.8e, %.8e\n", iter, freeEnergy, inCompMax);
		freeDiff=fabs(freeEnergy-freeOld);
        
		if(iter%50==0);//write_ph(phA,phB,wA,wB);
	}while(iter<maxIter&&(inCompMax>Sm1||freeDiff>Sm2));

	//fp=fopen("fe_end.dat","w");
       // fprintf(fp,"%d\n",iter);
	//fprintf(fp,"%10.8e, %10.8e, %10.8e, %10.8e, %10.8e, %e\n",freeEnergy,freeAB,freeW,freeS,freeWsurf,inCompMax);
	//fclose(fp);

	//write_ph(phA,phB,wA,wB);
	
	return freeEnergy;
}
//++++++++++++++++++++++++++++++++
double getConc(double *phlA,double *phlB,double *wA,double *wB)
{
	int i,j,k,iz,m;
	long ijk,ijkiz;
	double *qA,*qcA,*qB,*qcB;
	double ql,ffl,*qInt,qtmp;
	//MPI_Status status;
	
	qA=(double *)malloc(sizeof(double)*NxNyNz*(NsA+1));
	qcA=(double *)malloc(sizeof(double)*NxNyNz*(NsA+1));
        qB=(double *)malloc(sizeof(double)*NxNyNz*(dNsB+1));
        qcB=(double *)malloc(sizeof(double)*NxNyNz*(dNsB+1));
	qInt=(double *)malloc(sizeof(double)*NxNyNz);

	for(ijk=0;ijk<NxNyNz;ijk++)
	{
		qInt[ijk]=1.0;
	}
	
	sovDifFft(qA,wA,qInt,fA,NsA,1);  /* A(n-1)+A_star_B */
	sovDifFft(qcB,wB,qInt,dfB,dNsB,-1);

        for(ijk=0;ijk<NxNyNz;ijk++)
        {
                qInt[ijk]=qA[ijk*(NsA+1)+NsA];
		qtmp=qcB[ijk*(dNsB+1)];
		for(m=1;m<Narm;m++)qInt[ijk]*=qtmp;
        }

	sovDifFft(qB,wB,qInt,dfB,dNsB,1);      //fa to 0 for qcA

	for(ijk=0;ijk<NxNyNz;ijk++)
	{
		qInt[ijk]=qcB[ijk*(dNsB+1)];
		qtmp=qcB[ijk*(dNsB+1)];
		for(m=1;m<Narm;m++)qInt[ijk]*=qtmp;
	}

	sovDifFft(qcA,wA,qInt,fA,NsA,-1);      //fb to 1 for qB

	ql=0.0;
	for(ijk=0; ijk<NxNyNz; ijk++)
	{
		ql+=qB[ijk*(dNsB+1)+dNsB];
	}

	ql/=NxNyNz;
	ffl=ds0/ql;

	for(ijk=0; ijk<NxNyNz; ijk++)
	{
		phlA[ijk]=0.0;
		phlB[ijk]=0.0;

		for(iz=0;iz<=NsA;iz++)
		{
			ijkiz=ijk*(NsA+1)+iz;
			if(iz==0||iz==NsA)phlA[ijk]+=(0.50*qA[ijkiz]*qcA[ijkiz]);
			else phlA[ijk]+=(qA[ijkiz]*qcA[ijkiz]);
		}


		for(iz=0;iz<=dNsB;iz++)
		{
			ijkiz=ijk*(dNsB+1)+iz;
			if(iz==0||iz==dNsB)phlB[ijk]+=(0.50*qB[ijkiz]*qcB[ijkiz]);
			else phlB[ijk]+=(qB[ijkiz]*qcB[ijkiz]);
		}

		phlA[ijk]*=ffl;
		phlB[ijk]*=(Narm*ffl);
	}
	free(qA);
	free(qcA);
	free(qB);
	free(qcB);

	free(qInt);

	return ql;
}

void sovDifFft(double *g,double *w,double *qInt,
	double z,int ns,int sign)
{
    
	int i,j,k,iz,ns1;
	long ijk,ijkr;

	ns1=ns+1;	

	for(ijk=0;ijk<NxNyNz;ijk++)
    	{
        	wdz[ijk]=exp(-w[ijk]*ds2);
    	}
	
	if(sign==1)
	{
		for(ijk=0;ijk<NxNyNz;ijk++)
		{
			g[ijk*ns1]=qInt[ijk];
		}

		for(iz=1;iz<=ns;iz++)
		{
			for(ijk=0;ijk<NxNyNz;ijk++)
			{
				in[ijk]=g[ijk*ns1+iz-1]*wdz[ijk];
			}

			fftw_execute (p_forward);

        
			for(k=0;k<Nz;k++)for(j=0;j<Ny;j++)for(i=0;i<Nxh1;i++)
			{
				ijk=(long)((k*Ny+j)*Nxh1+i);
				ijkr=(long)((k*Ny+j)*Nx+i);
                
				out[ijk][0]*=kxyzdz[ijkr];	//out[].re or .im for fftw2
				out[ijk][1]*=kxyzdz[ijkr];	//out[][0] or [1] for fftw3
			}

			fftw_execute(p_backward);
			
			for(ijk=0;ijk<NxNyNz;ijk++)
				g[ijk*ns1+iz]=in[ijk]*wdz[ijk]/NxNyNz;
		}
	}
	else 
	{

		for(ijk=0;ijk<NxNyNz;ijk++)
		{
			g[ijk*ns1+ns] = qInt[ijk];
		}
		
		for(iz=ns-1;iz>=0;iz--)
		{
			for(ijk=0;ijk<NxNyNz;ijk++)
			{
				in[ijk]=g[ijk*ns1+iz+1]*wdz[ijk];
			}

			fftw_execute(p_forward);
			
			for(k=0;k<Nz;k++)for(j=0;j<Ny;j++)for(i=0;i<Nxh1;i++)
			{
				ijk=(long)((k*Ny+j)*Nxh1+i);
				ijkr=(long)((k*Ny+j)*Nx+i);
                
				out[ijk][0]*=kxyzdz[ijkr];
				out[ijk][1]*=kxyzdz[ijkr];
			}

			fftw_execute(p_backward);

			for(ijk=0;ijk<NxNyNz;ijk++)
				g[ijk*ns1+iz]=in[ijk]*wdz[ijk]/NxNyNz;
		}
	}
}

//********************Output configuration******************************

void write_ph(double *phA,double *phB,double *wA,double *wB)
{
	int i,j,k;
	long ijk;
	FILE *fp=fopen(phname,"w");

        fprintf(fp, "Nx=%d, Ny=%d, Nz=%d\n", Nx, Ny, Nz);
        fprintf(fp, "dx=%lf, dy=%lf, dz=%lf\n", dx, dy, dz);

	for(ijk=0;ijk<NxNyNz;ijk++)
	{
		fprintf(fp,"%lf %lf %lf %lf\n",phA[ijk],phB[ijk],wA[ijk],wB[ijk]);
	}

	fclose(fp);
}

int translation(double *density,double *density_trans,int direction,int dis){
int i,j,k;
	//printf("direction %d\n",direction);
	if(direction==1){
		for(k=0;k<Nz;k++)
		for(j=0;j<Ny;j++)
		for(i=0;i<Nx;i++){
			
			density_trans[i+Nx*j+k*Nx*Ny]=density[i+Nx*((j+dis+Ny)%Ny)+k*Nx*Ny];
			
		}

	}
	else if(direction==2){
		for(k=0;k<Nz;k++)
		for(j=0;j<Ny;j++)
		for(i=0;i<Nx;i++){
			density_trans[i+Nx*j+k*Nx*Ny]=density[i+Nx*((j-dis+Ny)%Ny)+k*Nx*Ny];
		}

	}
	else if(direction==3){
		for(k=0;k<Nz;k++)
		for(j=0;j<Ny;j++)
		for(i=0;i<Nx;i++){
			density_trans[i+Nx*j+k*Nx*Ny]=density[i+Nx*j+((k+dis+Nz)%Nz)*Nx*Ny];
		}

	}
	else if(direction==4){
		for(k=0;k<Nz;k++)
		for(j=0;j<Ny;j++)
		for(i=0;i<Nx;i++){
			density_trans[i+Nx*j+k*Nx*Ny]=density[i+Nx*j+((k-dis+Nz)%Nz)*Nx*Ny];
		}

	}
}
void translation_good(double *density,double *target){
	int i,j,k;
	double density_trans[NxNyNz];
	int tag,direction;
	double para_temp,para;
	para=distance_field(density,target);
	tag=0;
		while(tag==0){
			direction=1;
			translation(density,density_trans,direction,1);
			para_temp=distance_field(density_trans,target);
			
			if(para_temp<para){
				for(i=0;i<Nx*Ny*Nz;i++){density[i]=density_trans[i]; para= para_temp;}
			}
			else{
				tag=1;
			}
		}
		tag=0;
		while(tag==0){
			direction=2;
			translation(density,density_trans,direction,1);
			para_temp=distance_field(density_trans,target);
			//printf("%g \n",para_temp);
			if(para_temp<para){
				for(i=0;i<Nx*Ny*Nz;i++){density[i]=density_trans[i]; para= para_temp;}
			}
			else{
				tag=1;
			}
		}
		tag=0;
		while(tag==0){
			direction=3;
			translation(density,density_trans,direction,1);
			para_temp=distance_field(density_trans,target);
			//printf("%g \n",para_temp);
			if(para_temp<para){
				for(i=0;i<Nx*Ny*Nz;i++){density[i]=density_trans[i];para= para_temp; }
			}
			else{
				tag=1;
			}
		}
		tag=0;
		while(tag==0){
			direction=4;
			translation(density,density_trans,direction,1);
			para_temp=distance_field(density_trans,target);
			//printf("%g \n",para_temp);
			if(para_temp<para){
				for(i=0;i<Nx*Ny*Nz;i++){density[i]=density_trans[i];para= para_temp;}
			}
			else{
				tag=1;
			}
		}
	

}
double distance_field(double *phi,double *phi_target){
	int ijk;
	double distance,V;
	V=1;
	distance=0;
	for(ijk=0;ijk<NxNyNz;ijk++){
		distance+=(phi[ijk]-phi_target[ijk])*(phi[ijk]-phi_target[ijk])/V;
	}
	distance=sqrt(distance);
	return distance;

}

















