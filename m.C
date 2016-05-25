#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "m.h"
#define n 3
#define spacing 8
#define d 4
#define epsilon 0.24
#define beta 5.5
#define Ncor 50
#define Ncf 100
//storage of functions related to SU(3) matrix generator
#include "genRR.h"
//utility matrix functions
#include "utilM.h"
//misc. utility functions
#include "miscutil.h"
//update lattice, calculate S
#include "lattice.h"

int main()
{
  srand(time(NULL));
  //dimension of matrices (nxn)
  //int n=3;
  //dimensions of lattice d^4: 4x4x4x4 lattice
  //int d=4;
  //spacing=L/a
  //int spacing=4;
  
  ////intermediate result storage
  double **res1=(double**)malloc(n*n*sizeof(double*));
  double **res2=(double**)malloc(n*n*sizeof(double*));
  double **res3=(double**)malloc(n*n*sizeof(double*));
  double **newM=(double**)malloc(n*n*sizeof(double*));
  double **newM2=(double**)malloc(n*n*sizeof(double*));
  for(int i=0; i<n*n; i++)
    {
      newM[i]=(double*)calloc(2,sizeof(double));
      newM2[i]=(double*)calloc(2,sizeof(double));
      res1[i]=(double*)calloc(2,sizeof(double));
      res2[i]=(double*)calloc(2,sizeof(double));
      res3[i]=(double*)calloc(2,sizeof(double));
    }

  ////generate 100 SU(3) matrices
  double ***container=gen_rand_matrix();
  //return 0;
  //place the daggers of the random SU(3) matrices in the container as well
  container=(double***)realloc(container,2*100*sizeof(double**));
  for(int i=100; i<200; i++)
    {
      container[i]=(double**)malloc(n*n*sizeof(double*));
      newM=dagger(container[i-100],newM);
      for(int k=0; k<n*n;k++)
	{
	  container[i][k]=(double*)malloc(2*sizeof(double));
	  container[i][k][0]=newM[k][0];
	  container[i][k][1]=newM[k][1];
	}
    }
  //generate random matrices, check
  //double treal=0;
  //double timg=0;
  /*for(int i=0; i<100; i++)
    {
    treal=1.0/n*trace_real(container[i],n);
    timg=1.0/n*trace_img(container[i],n);
    printf("%.6f,%.6f\n",treal,timg);
    }*/
  //return 0;
  //print_mat(Times(dagger(container[148],res1,n),container[148],res2,n),n);
  //print_mat(Times(container[148],container[48],res2,n),n);
  //print_mat(container[48],n);
  //double *x=(double*)calloc(2,sizeof(double));
  //x=getdet(container[148],n,x);
  //printf("x[0]=%.6f\n", x[0]);
  //printf("x[1]=%.6f\n", x[1]);
  //x[0]=0;
  //x[1]=0;
  //x=getdet(container[48],n,x);
  //printf("x[0]=%.6f\n", x[0]);
  //printf("x[1]=%.6f\n", x[1]);
  //free(x);
  //return 0;

  //initialize the lattice
  double**** lattice=initialize_lat();
  //calculate the total number of plaquettes
  double s2=calculate_S(lattice);
  double total_plqts=d*(d-1)*pow(spacing,d)*.5;
  //printf("total plqts=%.6f\n", total_plqts);
  //printf("diff=%.6f\n",-s2/beta-total_plqts);
  //printf("axa avg=%.6f\n",calculate_S(lattice,-1,d,spacing,n)/total_plqts);
  //return 0;
  //printf("Well, S_init=%.6f\n",s2);

  //thermalize the lattice 10*Ncor times
  //for(int i=0; i<10*Ncor; i++)
  FILE * pFile;
  pFile=fopen ("results.txt","w");
  for(int i=0; i<10*Ncor; i++)
    {
      fprintf(pFile,"here we are=%d\n",10*Ncor-i);
      printf("here we are=%d\n",10*Ncor-i);
      fprintf(pFile,"axa avg=%.6f\n",calculate_S(lattice)/total_plqts/(-1*beta));
      printf("axa avg=%.6f\n",calculate_S(lattice)/total_plqts/(-1*beta));
      lattice=update(lattice,container);
    }


  //calculate the average plaquette, store it
  //double avg_plqt=0;
  double avg_plqt1=0;
  //update the links Ncor times, save S, repeat Ncf-1 times
  for(int i=0; i<Ncf; i++)
    {
      avg_plqt1=calculate_S(lattice)/total_plqts/(-1*beta);
      //avg_plqt+=avg_plqt/Ncf/total_plqts;
      fprintf(pFile,"avg_plqt(axa)(path number)=(%.6f)(%d)\n",avg_plqt1,i);
      printf("avg_plqt(axa)(path number)=(%.6f)(%d)\n",avg_plqt1,i);
      for(int j=0; j<Ncor; j++)
	{
	  lattice=update(lattice,container);
	}
    }

  //printf("average axa=%.6f\n",avg_plqt);

  //freeing of storage

  fclose(pFile);

  for(int i=0; i<200; i++)
    {
      for(int k=0; k<n*n; k++)
	{
	  free(container[i][k]);
	}
      free(container[i]);
    }
  free(container);
  for(int k=0; k<n*n; k++)
    {
      free(res1[k]);
      free(res2[k]);
      free(res3[k]);
      free(newM[k]);
      free(newM2[k]);
    }
  free(res1);
  free(res2);
  free(res3);
  free(newM);
  free(newM2);
  int pow=1;
  for(int k=0; k<d; k++)
    {
      pow*=(spacing);
    }
  for(int i=0; i<pow; i++)
    {
      for(int s=0; s<d; s++)
	{
	  for(int k=0; k<n*n; k++)
	    {
	      free(lattice[i][s][k]);
	    }
	  free(lattice[i][s]);
	}
      free(lattice[i]);
    }
  free(lattice);
  return 0;
}
