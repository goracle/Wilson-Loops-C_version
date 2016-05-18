#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double** minors(double **M, int n, int row, int column);

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

double plaq(double ****lattice, int coor, int nu, int mu, int n, int spacing, double result);

double**** update(double**** lattice, double*** container, int n, int d, int spacing, double beta);

int print_mat(double** M,int n);

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

double* times(double* z, double* y,double* result)
{
  result[0]=z[0]*y[0]-z[1]*y[1];
  result[1]=z[0]*y[1]+z[1]*y[0];
  return result;
}

double* getdet(double** M,int n, double* temp_fin)
{
  for(int column=0; column<n; column++)
    {
      if(n==2)
	{
	  double* z=(double*)calloc(2,sizeof(double));
	  double* t=(double*)calloc(2,sizeof(double));
	  z=times(M[1],M[2],z);
	  t=times(M[0],M[3],t);
	  temp_fin[0]=t[0]-z[0];
	  temp_fin[1]=t[1]-z[1];
	  free(z);
	  free(t);
	}
      else
	{
	  double y=-1.0;
	  double **mins=minors(M,n,0,column);
	  double* subdet=(double*)calloc(2,sizeof(double));
	  double* x=(double*)calloc(2,sizeof(double));
	  subdet=getdet(mins,n-1,subdet);
	  x=times(M[column],subdet,x);
	  free(subdet);
	  free(mins);
	  temp_fin[0]+=pow(y,column)*x[0];
	  temp_fin[1]+=pow(y,column)*x[1];
	  free(x);
	}
    }
  return temp_fin;
}

double getnorm(double* R,double* I,int n, int col=0)
{
  double rsum=0;
  double sum=0;
  for(int row=0; row<n; row++)
    {
      //printf("(a,b)=%.6f+I %.6f\n", R[n*row+col],I[n*row+col]);
      rsum+=pow(R[n*row+col],2)+pow(I[n*row+col],2);
      //printf("quantity we are adding=%.6f\n",(pow(R[n*row+col],2)+pow(I[n*row+col],2)));
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
      container[i]=(double**)malloc(n*n*sizeof(double*));
      for(int j=0; j<n*n; j++)
	{
	  container[i][j]=(double*)calloc(2,sizeof(double));
	}
    }
  //seed the random number generator
  srand(clock());
  for(int i=0; i<100; i++)
    {
      //make a hermitian matrix
      //allocate matrix memory
      double *R=(double *)calloc(n*n,sizeof(double));
      double *I=(double *)calloc(n*n,sizeof(double));
      double *H=(double *)calloc(n*n,sizeof(double));
      for(int k=0; k<n*n; k++)
	{
	  //allocate memory for indices, random entries
	  double r=0;
	  int a=0;
	  int b=0;
	  int k2=0;
	  //generate the random number
	  r=(double)rand()/((double)RAND_MAX/2)-1.0;
	  //get matrix element we are on (a,b)
	  a=k/n;
	  b=k-n*a;
	  //printf("k=%d",k);
	  //printf("a=%d",a);
	  //printf("b=%d",b);
	  //get transposed matrix element, (b,a) then assign to k2
	  k2=n*b+a;
	  //printf("k2=%d\n",k2);
	  //if this new k>old k, assign matrix element
	  if(k2>=k)
	    {
	      H[k]=r;
	      H[k2]=r;
	    }
	}
      
      //do 1+i*epsilon H
      //make a real matrix 
      for(int j=0; j<n; j++)
	{
	  R[(n+1)*(j)]=1; 
	}
      //and an imaginary matrix
      for(int j=0; j<n*n; j++)
	{
	  double temp=0;
	  temp=H[j]*epsilon;
	  I[j]=temp; 
	}
      free(H);

      //make it unitary
      
      //printf("make it unitary\n");
      for(int col=0; col<n; col++)
	{
	  //implement gram-schmidt
	  //make it orthogonal to previous columns
	  //printf("beginning gs procedure,%d\n",col);
	  for(int prev=col-1; prev>=0; prev--)
	    {
	      //get dot product
	      //printf("getting dot product,prev=%d\n",prev);
	      double real=0;
	      double imag=0;
	      for(int row=0; row<n; row++)
		{
		  imag+= -R[n*row+col]*(I[n*row+prev])+I[n*row+col]*(R[n*row+prev]);
		  real+= R[n*row+col]*(R[n*row+prev])+I[n*row+col]*(I[n*row+prev]);
		}
	      //subtract the projection of the prev. col onto the current
	      //from the current
	      //printf("subtracting projection\n");
	      for(int row=0; row<n; row++)
		{
		  R[n*row+col]-=real*(R[n*row+prev])-imag*(I[n*row+prev]);
		  I[n*row+col]-=real*(I[n*row+prev])+imag*(R[n*row+prev]);
		}
	    }
	  ////normalize the resulting column
	  //printf("normalize the result\n");
	  double rsum=0;
	  //get norm
	  double sum=0;
	  sum=getnorm(R,I,n,col);
	  //divide by sum, aka, the norm
	  //printf("divide by sum=%f\n",sum);
	  for(int row=0; row<n; row++)
	    {
	      double *tempreal=(double *)malloc(sizeof(double));
	      double *tempim=(double *)malloc(sizeof(double));
	      *tempim=*(I+n*row+col);
	      *tempreal=*(R+n*row+col);
	      *(I+n*row+col)=*tempim/(sum);
	      *(R+n*row+col)=*tempreal/(sum);
	      free(tempreal);
	      free(tempim);
	    }
	  //printf("new norm=%.6f\n",getnorm(R,I,n,col));
	}
      //get determinant
      double** M=(double**)malloc(n*n*sizeof(double*));
      for(int g=0; g<n*n; g++)
	{
	  M[g]=(double*)malloc(2*sizeof(double));
	  M[g][0]=0;
	  M[g][1]=0;
	}
      for(int ic=0; ic<n*n; ic++)
	{
	  M[ic][0]=R[ic];
	  M[ic][1]=I[ic];
	}
      double* determinant=(double*)calloc(2,sizeof(double));
      determinant=getdet(M,n,determinant);
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
	  double* res=(double*)malloc(2*sizeof(double));
	  res=times(M[d],c,res);
	  M[d][0]=res[0];
	  M[d][1]=res[1];
	  free(res);
	}

      free(c);
      free(determinant);
      //prints matrix we just generated 
      //print_mat(M,n); 

      //store the result
      for(int k=0; k<n*n; k++)
	{
	  for(int j=0;j<2;j++)
	    {
	      if(j==0)
		{
		  container[i][k][j]=M[k][0];
		}
	      if(j==1)
		{
		  container[i][k][j]=M[k][1];
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

int print_mat(double** M,int n)
{
  if(n!=3)
    return -1;
  printf("\nmarker\n");
  printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", M[0][0],M[0][1],M[1][0],M[1][1],M[2][0],M[2][1]);
  printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", M[3][0],M[3][1],M[4][0],M[4][1],M[5][0],M[5][1]);
  printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", M[6][0],M[6][1],M[7][0],M[7][1],M[8][0],M[8][1]);
  return 0;
}



int factorial(int d)
{
  if(d==1 || d<0)
    return d;
  return d*factorial(d-1);
}


double**** update(double**** lattice, double*** container, int n, int d, int spacing, double beta)
{
  //storage for coordinates
  int* x=(int*)calloc(d,sizeof(int));
  //storage for a copy of the link we are currently updating
  double **copy=(double**)malloc(n*n*sizeof(double*));
  //storage for intermediate results
  double **res1=(double**)malloc(n*n*sizeof(double*));
  double **res2=(double**)malloc(n*n*sizeof(double*));
  for(int i=0; i<n*n; i++)
    {
      copy[i]=(double*)calloc(2,sizeof(double));
      res1[i]=(double*)calloc(2,sizeof(double));
      res2[i]=(double*)calloc(2,sizeof(double));
    }
  //storage for staples
  double ***staples=(double***)malloc((d-1)*2*sizeof(double**));
  for(int j=0; j<(d-1)*2; j++)
    {
      staples[j]=(double**)malloc(n*n*sizeof(double*));
      for(int i=0; i<n*n; i++)
	staples[j][i]=(double*)calloc(2,sizeof(double));
    }
  //pow=number of coordinates in the lattice
  int pow=1;
  for(int k=0; k<d; k++)
    {
      pow*=(spacing+1);
    }
  //loop over coordinates
  for(int coor=0; coor<pow; coor++)
    {
      x=findcoord(d, spacing, coor, x);
      for(int mu=0; mu<d; mu++)
	{
	  //check if link even exists (none extend beyond edge of lattice)
	  if(x[mu]==spacing)
	    continue;
	  //number of staples
	  int nstap=0;
	  //increment of the coordinate in the mu direction by one unit
	  int mua=1;
	  for(int i=0; i<mu; i++)
	    {
	      mua*=spacing+1;
	    }
	  //loop over planes (mu,nu plane)
	  for(int nu=0; nu<d; nu++)
	    {
	      //index of staples to store the new staples of the mu,nu plane
	      int stap_start=nstap;
	      //can't have xx or yy plane, clearly nonsense
	      if(mu==nu)
		continue;
	      //increment of the coordinate in the nu direction by one unit
	      int nua=1;
	      for(int j=0; j<nu; j++)
		{
		  nua*=spacing+1;
		}
	      //check if we are on the boundary of the lattice,
	      //if not, we have an extra staple to compute
	      if(x[nu]!=spacing && x[nu]!=0)
		nstap++;
	      nstap++;
	      //j is our starting index
	      int j=stap_start;
	      //boundary, one staple
	      if(x[nu]==spacing)
		{
		  for(int v=0; v<n*n; v++)
		    {
		      res1=Times(dagger(lattice[coor+mua-nua][nu],n),dagger(lattice[coor-nua][mu],n),res1,n);
		      res2=Times(res1,lattice[coor-nua][nu],res2,n);
		      staples[j][v][0]=res2[v][0];
		      staples[j][v][1]=res2[v][1];
		    }
		}
	      //boundary, one staples
	      else if(x[nu]==0)
		{
		  res1=Times(dagger(lattice[coor+nua][mu],n),dagger(lattice[coor][nu],n),res1,n);
		  res2=Times(lattice[coor+mua][nu],res1,res2,n);
		  for(int v=0; v<n*n; v++)
		    {
		      staples[j][v][0]=res2[v][0];
		      staples[j][v][1]=res2[v][1];
		    }
		}
	      //else, midplane, 2 staples
	      else
		{
		  res1=Times(dagger(lattice[coor+mua-nua][nu],n),dagger(lattice[coor-nua][mu],n),res1,n);
		  res2=Times(res1,lattice[coor-nua][nu],res2,n);
		  for(int v=0; v<n*n; v++)
		    {
		      staples[j][v][0]=res2[v][0];
		      staples[j][v][1]=res2[v][1];
		    }
		  j++;
		  res1=Times(dagger(lattice[coor+nua][mu],n),dagger(lattice[coor][nu],n),res1,n);
		  res2=Times(lattice[coor+mua][nu],res1,res2,n);
		  for(int v=0; v<n*n; v++)
		    {
		      staples[j][v][0]=res2[v][0];
		      staples[j][v][1]=res2[v][1];
		    }
		}
	    }
	
	  //update the link
	  for(int tr=0; tr<10; tr++)
	    {

	      //select a random SU(n) matrix
	      srand(clock());
	      int ran=(int)((double)rand()/((double)RAND_MAX/100));
	      //form M-1
	      for(int v=0; v<n*n; v++)
		{
		  copy[v][0]=container[ran][v][0]; 
		  copy[v][1]=container[ran][v][1]; 
		  if(v%(n+1)==0)
		    copy[v][0]-=1;
		}
	      res2=Times(copy,lattice[coor][mu],res2,n);
	      	      //calculate dS (METROPOLIS PART)
	      double dS=0;
	      for(int j=0; j<nstap; j++)
		dS+=-beta*1.0/n*trace_real(Times(res2,staples[j],res1,n),n);
	      if(dS<0)
		{
		  res1=Times(container[ran],lattice[coor][mu],res1,n);
		  for(int k=0; k<n*n; k++)
		    {
		      lattice[coor][mu][k][0]=res1[k][0];
		      lattice[coor][mu][k][1]=res1[k][1];
		    }
		}
	      else
		{
		  srand(clock());
		  double ran2=(double)rand()/((double)RAND_MAX);
		  if(exp(-dS)>ran2)
		    {
		      res1=Times(container[ran],lattice[coor][mu],res1,n);
		      for(int k=0; k<n*n; k++)
			{
			  lattice[coor][mu][k][0]=res1[k][0];
			  lattice[coor][mu][k][1]=res1[k][1];
			}
		    }
		}
	    }
	}
    }
  free(copy);
  free(res1);
  free(res2);
  free(staples);
  free(x);
  return lattice;
}

double calculate_S(double ****lattice, double beta, int d, int spacing, int n)
{
  double S=0;
  int* x=(int*)calloc(d,sizeof(int));
  int pow=1;
  for(int k=0; k<d; k++)
    {
      pow*=(spacing+1);
    }
  for(int coor=0; coor<pow; coor++)
    {
      x=findcoord(d, spacing, coor, x);
      for(int mu=0; mu<d; mu++)
	{
	  if(x[mu]==spacing)
	    continue;
	  for(int nu=mu+1; nu<d; nu++)
	    {
	      if(x[nu]==spacing)
		continue;
	      double res=0;
	      S+=-beta*plaq(lattice, coor, mu, nu, n, spacing, res);
	    }
	}
    }
  free(x);
  return S;
}
//coordinate,direction of link, index of nxn matrix, real/imag part
double plaq(double ****lattice, int coor, int mu, int nu, int n, int spacing, double result=0)
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
  result=1.0/n*trace_real(Times(res1,res2,res3,n),n);
  free(res1);
  free(res2);
  free(res3);
  return result;
}

double **Times(double** M, double** N, double** res, int n)
{
  double* res4=(double*)malloc(2*sizeof(double));
  for(int i=0; i<n; i++)
    {
      for(int j=0; j<n; j++)
	{
	  res[i*n+j][0]=0;
	  res[i*n+j][1]=0;
	  for(int k=0; k<n; k++)
	    {
	      res4=times(M[i*n+k],N[k*n+j],res4);
	      res[i*n+j][0]+=res4[0]; 
	      res[i*n+j][1]+=res4[1]; 
	    } 
	}
    }
  free(res4);
  return res;
}
//confident
double trace_real(double **M, int n)
{
  double sum=0;
  for(int k=0; k<n; k++)
    {
      sum+=M[k*(n+1)][0];
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
	continue;
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

int main()
{
  //global constants
  double epsilon=.24;
  double beta=5.5;
  //double beta=1;
  int Ncor=50;
  int Ncf=100;
  //dimension of matrices (nxn)
  int n=3;
  //dimensions of lattice d^4: 4x4x4x4 lattice
  int d=4;
  //spacing=L/a
  int spacing=8;
  
  
  
  double ***container=gen_rand_matrix(n,epsilon);

  printf("\n imaginary part of 3,3 element %.6f\n",container[99][8][1]);

  //initialize the lattice
  double**** lattice=initialize_lat(d,spacing,n);
  double s2=calculate_S(lattice, beta, d, spacing,n);
  double total_plqts=s2;

  printf("Well, S_init=%.6f\n",s2);
  //indices are x,y,z,t,direction,matrix index, real part, imaginary part
  //loop over x



    
  //thermalize the lattice 10*Ncor times
  for(int i=0; i<10*Ncor; i++)
    {
      update(lattice,container,n,d,spacing,beta);
      printf("here we are=%d\n",10*Ncor-i);
      printf("axa avg=%.6f\n",calculate_S(lattice,beta,d,spacing,n)/s2);
    }


  //calculate the average plaquette, store it
  double avg_plqt=0;
  //update the links Ncor times, save S, repeat Ncf-1 times
  for(int i=0; i<Ncf; i++)
    {
      avg_plqt+=calculate_S(lattice, 1/total_plqts, d, spacing,n)/Ncf;
      printf("avg_plqt(axa)=%.6f\n",avg_plqt);
      for(int j=0; j<Ncor; j++)
	update(lattice,container,n,d,spacing,beta);
    }

  printf("\n average axa=%.6f\n",avg_plqt);
  //loop over links

  ///update the link

  //calculate the action
  free(lattice);
  free(container);
  return 0;
}
