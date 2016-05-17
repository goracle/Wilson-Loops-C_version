#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double** minors(double **M, int n, int row, int column);

int* baseconv(int start, int basefrom, int baseto);

double* getdet(double** M,int n);

double* times(double* z, double* y, double* result);

double getnorm(double* R,double* I,int n, int col);

double*** gen_rand_matrix(int n, double epsilon);

double**** initialize_lat(int d, int spacing, int n);

double calculate_S(double ****lattice, double beta, int d, int spacing, int n);

int* findcoord(int d, int spacing, int coor, int* x);

double** dagger(double** M, int n);

double trace_real(double **M, int n);

double **Times(double** M, double** N, double** res, int n);

double plaq(double ****lattice, int coor, int nu, int mu, int n, int spacing);

double**** update(double**** lattice, double*** container, int n, int d, int spacing, double beta);

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
	  double* z=times(M[1],M[2]);
	  double* t=times(M[0],M[3]);
	  temp_fin[0]=t[0]-z[0];
	  temp_fin[1]=t[1]-z[1];
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
	  M[g][0]=0;
	  M[g][1]=0;
	}
      for(int ic=0; ic<n*n; ic++)
	{
	  M[ic][0]=R[ic];
	  M[ic][1]=I[ic];
	}
      double* determinant=getdet(M,n);
      //now that we have the determinant(aka the phase), we must divide the matrix by the determinant^n
      double phase=0;
      if(determinant[0]!=0)
	{
	  phase=-1.0*atan(determinant[1]/determinant[0])/n;
	}
      else
	{
	  phase=-1.0*asin(determinant[1])/n;
	}
      
      double *c=(double*)malloc(2*sizeof(double));
      c[0]=cos(phase);
      c[1]=sin(phase);
      //make the unitary matrix SU(N) with n=N
      for(int d=0; d<n*n; d++)
	{
	  M[d]=times(M[d],c);
	}

      free(c);
      free(determinant);
      
 

      printf("\nmarker\n");
      printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", M[0][0],M[0][1],M[1][0],M[1][1],M[2][0],M[2][1]);
      printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", M[3][0],M[3][1],M[4][0],M[4][1],M[5][0],M[5][1]);
      printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", M[6][0],M[6][1],M[7][0],M[7][1],M[8][0],M[8][1]);
      //store the result
      for(int j=0;j<2;j++)
	{
	  for(int k=0; k<n*n; k++)
	    {
	      if(j==0)
		{
		  *(*(*(container+i)+j)+k)=M[k][0];
		}
	      if(j==1)
		{
		  *(*(*(container+i)+j)+k)=M[k][1];
		}

	    }
	}
      //free the memory
      free(M);
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
  int n=3;
  //dimensions of lattice d^4: 4x4x4x4 lattice
  int d=4;
  //spacing=L/a
  int spacing=8;
  
  
  
  double ***container=gen_rand_matrix(n,epsilon);

  printf("\n imaginary part of 3,3 element %.6f\n",container[99][1][8]);

  //initialize the lattice
  double**** lattice=initialize_lat(d,spacing,n);
  double S=calculate_S(lattice, beta, d, spacing,n);

  printf("Well, S_init=%.6f\n",S);
  //indices are x,y,z,t,direction,matrix index, real part, imaginary part
  //loop over x



    
  //thermalize the lattice 10*Ncor times
  for(int i=0; i<10*Ncor; i++)
    {
      update(lattice,container,n,d,spacing,beta);
      printf("here we are\n");
    }


  //calculate the action, store it
  double avgS=0;
  //update the links Ncor times, save S, repeat Ncf-1 times
  for(int i=0; i<Ncf; i++)
    {
      avgS+=calculate_S(lattice, beta, d, spacing,n)/Ncf;
      for(int j=0; j<Ncor; j++)
	update(lattice,container,n,d,spacing,beta);
    }

  printf("\n average S=%.6f\n",avgS);
  //loop over links

  ///update the link

  //calculate the action
  free(lattice);
  free(container);
  return 0;
}



double**** update(double**** lattice, double*** container, int n, int d, int spacing, double beta)
{
  int* x=(int*)malloc(d*sizeof(int));
  int pow=1;
  for(int k=0; k<d; k++)
    {
      pow*=(spacing+1);
      x[d]=0;
    }
  for(int coor=0; coor<pow; coor++)
    {
      x=findcoord(d, spacing, coor, x);
      for(int mu=0; mu<d; mu++)
	{
	  if(x[mu]==spacing)
	    continue;
	  //compute the staples
	  int nstap=0;
	  double ***staples=(double***)malloc(nstap*sizeof(double**));
	  int mua=1;
	  for(int i=0; i<mu; i++)
	    {
	      mua*=spacing+1;
	    }
	  for(int nu=0; nu<d; nu++)
	    {
	      int stap_start=nstap;
	      if(mu==nu)
		continue;
	      int nua=1;
	      for(int j=0; j<nu; j++)
		{
		  nua*=spacing+1;
		}
	      if(x[nu]!=spacing && x[nu]!=0)
		nstap++;
	      nstap++;
	      printf("coor=%d\n",coor);
	      //increase the size of staples
	      staples=(double***)realloc(staples,nstap*sizeof(double**));
	      for(int j=stap_start; j<nstap; j++)
		{
		  double **res1=(double**)malloc(n*n*sizeof(double*));
		  double **res2=(double**)malloc(n*n*sizeof(double*));
		  for(int i=0; i<n*n; i++)
		    {
		      res1[i]=(double*)calloc(2,sizeof(double));
		      res2[i]=(double*)calloc(2,sizeof(double));
		    }

		  staples[j]=(double**)malloc(n*n*sizeof(double*));
		  if(x[nu]==spacing)
		    staples[j]=Times(Times(dagger(lattice[coor+mua-nua][nu],n),dagger(lattice[coor-nua][mu],n),res1,n),lattice[coor-nua][nu],res2,n);
		  if(x[nu]==0)
		    staples[j]=Times(lattice[coor+mua][nu],Times(dagger(lattice[coor+nua][mu],n),dagger(lattice[coor][nu],n),res1,n),res2,n);
		  else
		    {
		      staples[j]=Times(Times(dagger(lattice[coor+mua-nua][nu],n),dagger(lattice[coor-nua][mu],n),res1,n),lattice[coor-nua][nu],res2,n);
		      j++;
		      staples[j]=Times(lattice[coor+mua][nu],Times(dagger(lattice[coor+nua][mu],n),dagger(lattice[coor][nu],n),res1,n),res2,n);;
		    }
		  free(res1);
		  free(res2);
		}
	    }
	
	  //update the link

	  double **copy=(double**)malloc(n*n*sizeof(double*));
	  for(int j=0; j<n*n; j++)
	    {
	      copy[j]=(double*)calloc(2,sizeof(double));
	      copy[j][0]=lattice[coor][mu][j][0];
	      copy[j][1]=lattice[coor][mu][j][1];
	    }
	  //select 10 random SU(n) matrices
	  srand(time(NULL));
	  for(int i=0; i<10; i++)
	    {
	      double **res2=(double**)malloc(n*n*sizeof(double*));
	      for(int j=0; j<n*n; j++)
		{
		  res2[j]=(double*)calloc(2,sizeof(double));
		}
	      int ran=(int)((double)rand()/((double)RAND_MAX/100));
	      printf("ran=%d\n",ran);
	      copy=Times(container[ran],copy,res2,n);
	    }
	  //calculate dS (METROPOLIS PART)
	  double dS=0;
	  for(int j=0; j<nstap; j++)
	    {
	      double **res2=(double**)malloc(n*n*sizeof(double*));
	      for(int i=0; i<n*n; i++)
		{
		  res2[i]=(double*)calloc(2,sizeof(double));
		}
	      dS+=-beta*1.0/3.0*trace_real(Times(staples[j],copy,res2,n),n);
	      dS-=-beta*1.0/3.0*trace_real(Times(staples[j],lattice[coor][mu],res2,n),n);
	      free(res2);
	    }
	  free(staples);
	  if(dS<0)
	    {
	      for(int k=0; k<n*n; k++)
		{
		  lattice[coor][mu][k][0]=copy[k][0];
		  lattice[coor][mu][k][1]=copy[k][1];
		}
	    }
	  else
	    {
	      srand(time(NULL)); 
	      double ran2=(double)rand()/((double)RAND_MAX);
	      if(exp(-dS)>ran2)
		{
		for(int k=0; k<n*n; k++)
		  {
		    lattice[coor][mu][k][0]=copy[k][0];
		    lattice[coor][mu][k][1]=copy[k][1];
		  }
		}
	    }
	  free(copy);
	}
    }
  free(x);
  return lattice;
}


double calculate_S(double ****lattice, double beta, int d, int spacing, int n)
{
  double S=0;
  int* x=(int*)malloc(d*sizeof(int));
  int pow=1;
  for(int k=0; k<d; k++)
    {
      pow*=(spacing+1);
      x[d]=0;
    }
  for(int coor=0; coor<pow; coor++)
    {
      x=findcoord(d, spacing, coor, x);
      for(int nu=0; nu<d; nu++)
	{
	  if(x[nu]==spacing)
	    continue;
	  for(int mu=nu+1; mu<d; mu++)
	    {
	      if(x[mu]==spacing)
		continue;
	      S+=-beta*plaq(lattice, coor, nu, mu, n, spacing);
	    }
	}
    }
  free(x);
  return S;
}

//coordinate,direction of link, index of nxn matrix, real/imag part
double plaq(double ****lattice, int coor, int nu, int mu, int n, int spacing)
{
  int nua=1;
  int mua=1;
  for(int i=0; i<nu; i++)
    {
      nua*=spacing+1;
    }
  for(int i=0; i<mu; i++)
    {
      mua*=spacing+1;
    }
  double **res1=(double**)malloc(n*n*sizeof(double*));
  double **res2=(double**)malloc(n*n*sizeof(double*));
  double **res3=(double**)malloc(n*n*sizeof(double*));
  for(int i=0; i<n*n; i++)
    {
      res1[i]=(double*)calloc(2,sizeof(double));
      res2[i]=(double*)calloc(2,sizeof(double));
      res3[i]=(double*)calloc(2,sizeof(double));
    }
  res1=Times(lattice[coor][mu],lattice[coor+mua][nu],res1,n);
  res2=Times(dagger(lattice[coor+nua][mu], n),dagger(lattice[coor][nu], n),res2,n);
  double result=1.0/3.0*trace_real(Times(res1,res2,res3,n),n);
  free(res1);
  free(res2);
  free(res3);
  return result;
}

double **Times(double** M, double** N, double** res, int n)
{
  double ***A=(double***)malloc(n*sizeof(double**));
  double ***B=(double***)malloc(n*sizeof(double**));
  for(int i=0; i<n; i++)
    {
      A[i]=(double**)malloc(n*sizeof(double*));
      B[i]=(double**)malloc(n*sizeof(double*));
      for(int j=0; j<n; j++)
	{
	  B[i][j]=(double*)calloc(2,sizeof(double));
	  A[i][j]=(double*)calloc(2,sizeof(double));
	}
    }
  for(int k=0; k<n*n; k++)
    {
      int i=k/n; 
      int j=k-i*n;
      A[i][j][0]=M[k][0];
      A[i][j][1]=M[k][1];
      B[i][j][0]=N[k][0];
      B[i][j][1]=N[k][1];
    }
  for(int i=0; i<n; i++)
    {
      for(int j=0; j<n; j++)
	{
	  for(int k=0; k<n; k++)
	    {
	      res[i*n+j][0]+=times(A[i][k],B[k][j])[0]; 
	      res[i*n+j][1]+=times(A[i][k],B[k][j])[1]; 
	    } 
	}
    }
  free(A);
  free(B);
  return res;
}

double trace_real(double **M, int n)
{
  double sum=0;
  for(int k=0; k<n*n; k+=n+1)
    {
      sum+=M[k][0];
    }
  return sum;
}

double** dagger(double** M, int n)
{
  for(int k=0; k<n*n; k++)
    {
      int i=k/n;
      int j=k-i*n;
      int kprime=j*n+i;
      if(kprime==k)
	M[k][1]*=-1;
      else if(kprime>k)
	{
	  double temp0=M[k][0];
	  double temp1=M[k][1];
	  M[k][0]=M[kprime][0];
	  M[k][1]=-1*M[kprime][1];
	  M[kprime][0]=temp0;
	  M[kprime][1]=-1*temp1;
	}
      else
	break;
    }
  return M;
}

double**** initialize_lat(int d, int spacing, int n)
{
  int pow=1;
  for(int k=0; k<d; k++)
    {
      pow*=spacing+1;
    }
  double**** lattice=(double****)malloc(pow*sizeof(double***));
  //coordinate,direction of link, index of nxn matrix, real/imag part
  int* x=(int*)malloc(d*sizeof(int));
  for(int coor=0; coor<pow; coor++)
    {
      lattice[coor]=(double***)malloc(d*sizeof(double**));
      for(int dir=0; dir<d; dir++)
	{
	  lattice[coor][dir]=(double**)malloc(n*n*sizeof(double*));
	  for(int index=0; index<n*n; index++)
	    {
	      x=findcoord(d, spacing, coor, x);
	      lattice[coor][dir][index]=(double*)malloc(2*sizeof(double));
	      //below: first condition sees if we are on diagonal
	      //second condition tests if we are on the border of the lattice
	      //in which case we create no dangling links
	      if(index%(n+1)==0 && x[dir]!=(spacing))
		{
		  lattice[coor][dir][index][0]=1;
		  lattice[coor][dir][index][1]=0;
		}
	      else
		{
		  lattice[coor][dir][index][0]=0;
		  lattice[coor][dir][index][1]=0;
		}
	    }
	}
    }
  free(x);
  return lattice;
}

int* findcoord(int d, int spacing, int coor, int* x)
{
  int temp=0;
  int coors=coor;
  while(coors>=0 && temp<d)
    {
      x[temp]=coors-(coors/(1+spacing))*(spacing+1);
      coors/=(spacing+1);
      temp++;
    }
  return x;
}
