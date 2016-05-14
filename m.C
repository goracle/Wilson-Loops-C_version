#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double** minors(double **M, int n, int row, int column);

double* getdet(double** M,int n);

double* times(double* z, double* y, double* result);

double getnorm(double* R,double* I,int n, int col);

double*** gen_rand_matrix(int n, double epsilon);

double** minors(double **M, int n, int row, int column)
{
  double **minorr=(double**)malloc((n-1)*(n-1)*sizeof(double*));
  for(int g=0; g<(n-1)*(n-1); g++)
    {
      *(minorr+g)=(double*)malloc(2*sizeof(double));
    }
  int newstart=0;
  for(int i=0; i<(n-1)*(n-1); i++)
    {
      for(int j=newstart; j<n*n; j++)
	{
	  if(j%n==column || j/n==row)
	    {
	      continue;
	    }
	  minorr[i][1]=M[j][1];
	  minorr[i][0]=M[j][0];
	  newstart=j+1;
	  break;
	}
    }
  return minorr;
}

double* times(double* z, double* y)
{
  double* result=(double *)malloc(2*sizeof(double));
  result[0]=0;
  result[1]=0;
  double y0=y[0];
  double z0=z[0];
  result[0]+=z[0]*y[0]-z[1]*y[1];
  result[1]+=z[0]*y[1]+z[1]*y[0];
  return result;
}

double* getdet(double** M,int n)
{
  double* temp_fin=(double*)malloc(2*sizeof(double));
  temp_fin[0]=0;
  temp_fin[1]=0;
  for(int column=0; column<n; column++)
    {
      if(n==2)
	{
	  double* t=times(M[0],M[3]);
	  temp_fin[0]=t[0];
	  temp_fin[1]=t[1];
	  free(t);
	}
      else
	{
	  double y=-1.0;
	  double **mins=minors(M,n,0,column);
	  double* subdet=getdet(mins,n-1);
	  double *x=times(M[column],subdet);
	  free(mins);
	  temp_fin[0]+=pow(y,column)*x[0];
	  temp_fin[1]+=pow(y,column)*x[1];
	  free(x);
	}
    }
  return temp_fin;
}

double getnorm(double* R,double* I,int n, int col)
{
  double rsum=0;
  double sum;
  for(int row=0; row<n; row++)
    {
      printf("(a,b)=%.6f+I %.6f\n", (*(R+n*row+col)),*(I+n*row+col));
      rsum+=pow((*(R+n*row+col)),2)+pow((*(I+n*row+col)),2);
      printf("quantity we are adding=%.6f\n",pow((*(R+n*row+col)),2)+pow((*(I+n*row+col)),2));
    }
  sum=pow(rsum,.5);
  return sum;
}

double*** gen_rand_matrix(int n, double epsilon)
{
  //generate 100 SU(3) matrices
  //make a container
  double ***container=(double***)(malloc(100*sizeof(double**)));
  for(int i=0; i<100; i++)
    {
      *(container+i)=(double**)malloc(2*n*n*sizeof(double*));
      for(int j=0; j<n*n; j++)
	{
	  *(*(container+i)+j)=(double*)malloc(2*sizeof(double));
	}
    }
  //seed the random number generator
  srand(time(NULL));
  for(int i=0; i<100; i++)
    {
      //make a hermitian matrix
      //allocate matrix memory
      double *R=(double *)malloc(n*n*sizeof(double));
      double *I=(double *)malloc(n*n*sizeof(double));
      double *H=(double *)malloc(n*n*sizeof(double));
      for(int k=0; k<n*n; k++)
	{
	  //allocate memory for indices, random entries
	  double *r=(double *)malloc(sizeof(double));
	  int *a=(int *)malloc(sizeof(int));
	  int *b=(int *)malloc(sizeof(int));
	  int *k2=(int *)malloc(sizeof(int));
	  //generate the random number
	  *r=(double)rand()/((double)RAND_MAX/2)-1.0;
	  //get matrix element we are on (a,b)
	  *a=k/n;
	  *b=k-n*(*a);
	  printf("k=%d",k);
	  printf("a=%d",*a);
	  printf("b=%d",*b);
	  //get transposed matrix element, (b,a) then assign to k2
	  *k2=n*(*b)+*a;
	  printf("k2=%d\n",*k2);
	  //if this new k>old k, assign matrix element
	  if(*k2>=k)
	    {
	      *(H+k)=*r;
	      *(H+*k2)=*r;
	    }
	  free(r);
	  free(a);
	  free(b);
	  free(k2);
	}
      
      //do 1+i*epsilon H
      //make a real matrix 
      for(int j=0; j<n; j++)
	{
	  *(R+(n+1)*(j))=1; 
	}
      //and an imaginary matrix
      for(int j=0; j<n*n; j++)
	{
	  double *temp=(double *)malloc(sizeof(double));
	  *temp=*(H+j)*epsilon;
	  *(I+j)=*temp; 
	  free(temp);
	}
      free(H);

      //make it unitary
      
      printf("make it unitary\n");
      for(int col=0; col<n; col++)
	{
	  //implement gram-schmidt
	  //make it orthogonal to previous columns
	  printf("beginning gs procedure,%d\n",col);
	  for(int prev=col-1; prev>=0; prev--)
	    {
	      //get dot product
	      printf("getting dot product,prev=%d\n",prev);
	      double *real=(double *)(malloc(sizeof(double)));
	      double *imag=(double *)(malloc(sizeof(double)));
	      for(int row=0; row<n; row++)
		{
		  *imag+= -*(R+n*row+col)*(*(I+n*row+prev))+*(I+n*row+col)*(*(R+n*row+prev));
		  *real+= *(R+n*row+col)*(*(R+n*row+prev))+*(I+n*row+col)*(*(I+n*row+prev));
		}
	      //subtract the projection of the prev. col onto the current
	      //from the current
	      printf("subtracting projection\n");
	      for(int row=0; row<n; row++)
		{
		  *(R+n*row+col)-=*real*(*(R+n*row+prev))-*imag*(*(I+n*row+prev));
		  *(I+n*row+col)-=*real*(*(I+n*row+prev))+*imag*(*(R+n*row+prev));
		}
	      free(imag);
	      free(real);
	    }
	  ////normalize the resulting column
	  printf("normalize the result\n");
	  double rsum=0;
	  //get norm
	  double *sum=(double *)malloc(sizeof(double));
	  *sum=getnorm(R,I,n,col);
	  //divide by sum, aka, the norm
	  printf("divide by sum=%f\n",*sum);
	  for(int row=0; row<n; row++)
	    {
	      double *tempreal=(double *)malloc(sizeof(double));
	      double *tempim=(double *)malloc(sizeof(double));
	      *tempim=*(I+n*row+col);
	      *tempreal=*(R+n*row+col);
	      *(I+n*row+col)=*tempim/(*sum);
	      *(R+n*row+col)=*tempreal/(*sum);
	      free(tempreal);
	      free(tempim);
	    }
	  printf("new norm=%.6f\n",getnorm(R,I,n,col));
	  free(sum);
	}
      //get determinant
      double** M=(double**)malloc(n*n*sizeof(double*));
      for(int g=0; g<n*n; g++)
	{
	  *(M+g)=(double*)malloc(2*sizeof(double));
	}
      for(int ic=0; ic<n*n; ic++)
	{
	  M[ic][0]=R[ic];
	  M[ic][1]=I[ic];
	}
      double* determinant=getdet(M,n);
      printf("determinant=%.6f+i %.6f\n",determinant[0],determinant[1]);
      double *determinant_star=(double*)malloc(2*sizeof(double));
      determinant_star[0]=determinant[0];
      determinant_star[1]=-1*determinant[1];
      double powT=-1.0/n;
      double normdet=pow(times(determinant,determinant_star)[0],powT);
      double phase;
      if(determinant[0]!=0)
	{
	  phase=-1.0*acos(determinant[0]/normdet)/n;
	}
      else if(determinant[1]!=0)
	{
	  phase=-1.0*asin(determinant[1]/normdet)/n;
	}
      else
	{
	  printf("determinant=0, skipping this matrix\n");
	  i--;
	  continue;
	}
      double *c=(double*)malloc(2*sizeof(double));
      c[0]=normdet*cos(phase);
      c[1]=normdet*sin(phase);
      printf("c=(%.6f,%.6f)",c[0],c[1]);
      for(int d=0; d<n*n; d++)
	{
	  M[d]=times(M[d],c);
	}

     double* determinant2=getdet(M,n);
     printf("new determinant=%.6f+i %.6f\n",determinant2[0],determinant2[1]);
     free(determinant2);
     free(determinant);

     free(c);
 

      printf("\nmarker\n");
      printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", R[0],I[0],R[1],I[1],R[2],I[2]);
      printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", R[3],I[3],R[4],I[4],R[5],I[5]);
      printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", R[6],I[6],R[7],I[7],R[8],I[8]);
      //store the result
      for(int j=0;j<2;j++)
	{
	  for(int k=0; k<n*n; k++)
	    {
	      if(j==0)
		{
		  *(*(*(container+i)+j)+k)=*(R+k);
		}
	      if(j==1)
		{
		  *(*(*(container+i)+j)+k)=*(I+k);
		}

	    }
	}
      //free the memory
      free(I);
      free(R);
      
    }
  return container;
}

int main()
{
  //global constants
  double epsilon=.24;
  double beta=5.5;
  int Ncor=50;
  int Ncf=100;
  //dimension of matrices (nxn)
  int n=2;
  
  double ***container=gen_rand_matrix(n,epsilon);

  printf("\n imaginary part of 3,3 element %.6f\n",container[99][1][8]);

  free(container);

  double** M=(double**)malloc(n*n*sizeof(double*));


  for(int g=0; g<n*n; g++)
    {
      M[g]=(double*)malloc(2*sizeof(double));
    }
  printf("seg\n");
  int j=0;
  for(double i=0; i<n*n; i+=1.0)
    {
      printf("seg(%d)\n",j);
      M[j][0]=i;
      j++;
    }
  double **m=minors(M,n,1,1);

  for(double i=0; i<(n-1)*(n-1); i+=1.0)
    {
      int k=(int)i;
      printf("m[%d]=%d\n",i,m[k][0]); 
    }
  free(m);
  free(M);

  
  //thermalize the lattice 10*Ncor times

  //calculate the action, store it

  //update the links Ncor times, save S, repeat Ncf-1 times

  //loop over links

  ///update the link

  //calculate the action
  
  
  return 0;
}
