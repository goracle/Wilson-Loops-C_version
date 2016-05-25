#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double** minors(double **M, int n, int row, int column, double** minorr);

double* getdet(double** M,int n, double* temp_fin);

double trace_img(double **M, int n);
  
double* times(double* z, double* y, double* result);

int testarg(int arg, int d, int spacing);

double getnorm(double** M,int n, int col);

double getnorm2(double* a, double *result);

double*** gen_rand_matrix(int n, double epsilon);

double**** initialize_lat(int d, int spacing, int n);

double calculate_S(double ****lattice, double beta, int d, int spacing, int n);

double calculate_S21(double ****lattice, double beta, int d, int spacing, int n);

int* findcoord(int d, int spacing, int coor, int* x);

double** dagger(double** M, double** newM, int n);

double trace_real(double **M, int n);

double **Times(double** M, double** N, double** res, int n);

double plaq(double ****lattice, int coor, int nu, int mu, int n, int spacing,int d, double result );

double plaq21(double ****lattice, int coor, int nu, int mu, int n, int spacing,int d, double result );

double**** update(double**** lattice, double*** container, int n, int d, int spacing, double beta);

int print_mat(double** M,int n);

double* star(double* x, double* result);

double** minors(double **M, int n, int row, int column, double** minorr)
{

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

double* times(double* y, double* z,double* result)
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
	  double **minorr=(double**)malloc((n-1)*(n-1)*sizeof(double*));
	  for(int g=0; g<(n-1)*(n-1); g++)
	    {
	      minorr[g]=(double*)malloc(2*sizeof(double));
	    }
	  minorr=minors(M,n,0,column,minorr);
	  double* subdet=(double*)calloc(2,sizeof(double));
	  double* x=(double*)calloc(2,sizeof(double));
	  subdet=getdet(minorr,n-1,subdet);
	  x=times(M[column],subdet,x);
	  free(subdet);
	  for(int g=0; g<(n-1)*(n-1); g++)
	    {
	      free(minorr[g]);
	    }
	  free(minorr);
	  temp_fin[0]+=pow(y,column)*x[0];
	  temp_fin[1]+=pow(y,column)*x[1];
	  free(x);
	}
    }
  return temp_fin;
}

double getnorm(double** M,int n, int col=0)
{
  double rsum=0;
  double sum=0;
  for(int row=0; row<n; row++)
    {
      //printf("(a,b)=%.6f+I %.6f\n", R[n*row+col],I[n*row+col]);
      rsum+=pow(M[n*row+col][0],2)+pow(M[n*row+col][1],2);
      //printf("quantity we are adding=%.6f\n",(pow(R[n*row+col],2)+pow(I[n*row+col],2)));
    }
  sum=pow(rsum,.5);
  return sum;
}

double getnorm2(double* a, double *result)
{
  return times(star(a,result),a,result)[0];
}

double* star(double* x, double* result)
{
  result[0]=x[0];
  result[1]=-x[1];
  return result;
}

double*** gen_rand_matrix(int n, double epsilon)
{
  //generate 100 SU(3) matrices
  //make a container
  double ***container=(double***)(malloc(100*sizeof(double**)));
  //seed the random number generator
  //storage of intermediate results
  double* product=(double*)malloc(2*sizeof(double));
  double* produ=(double*)malloc(2*sizeof(double));
  double* produ1=(double*)malloc(2*sizeof(double));
  double *c=(double*)malloc(2*sizeof(double));
  double *x=(double*)malloc(2*sizeof(double));
  double *y=(double*)malloc(2*sizeof(double));
  double *z=(double*)malloc(2*sizeof(double));
  double* determinant=(double*)calloc(2,sizeof(double));
  double* res=(double*)malloc(2*sizeof(double));
  for(int i=0; i<100; i++)
    {
      //make a hermitian matrix with entries between -1 and 1
      double** M=(double **)malloc(n*n*sizeof(double*));
      for(int k=0; k<n*n; k++)
	{
	  //get matrix element we are on (a,b)
	  int a=k/n;
	  int b=k-a*n;
	  //transposed k
	  int k2=b*n+a;
	  if(k2>=k)
	    {
	      M[k]=(double*)calloc(2,sizeof(double));
	      //generate the random number
	      //srand(clock());
	      double r=(double)rand()/((double)RAND_MAX/2)-1.0;
	      if(k!=k2)
		{
		  M[k2]=(double*)calloc(2,sizeof(double));
		  M[k2][0]=r; 
		  M[k2][1]=0;
		}
	      M[k][0]=r; 
	      M[k][1]=0;
	    }
	}
      
      //do 1+i*epsilon H
      //make a real matrix 
      for(int j=0; j<n*n; j++)
	{
	  M[j][1]=M[j][0]*epsilon;
	  M[j][0]=0;
	  if(j%(n+1)==0)
	    M[j][0]+=1;
	}
      //make it unitary
      //printf("make it unitary\n");
      for(int col=0; col<n; col++)
	{
	  //implement gram-schmidt
	  //make it orthogonal to previous columns
	  //printf("beginning gs procedure,%d\n",col);
	  if(col>0 && (n!=3 || col!=2))//then we aren't on the first
	    //column, so we must subtract the projection
	    //we also aren't on the third column (which is col1 cross col2)
	    {
	      for(int prev=col-1; prev>=0; prev--)
		{
		  //get dot product
		  //printf("getting dot product,prev=%d\n",prev);
		  product[0]=0;
		  product[1]=0;
		  for(int row=0; row<n; row++)
		    {
		      produ=times(star(M[n*row+prev],produ1),
				  M[n*row+col],produ);
		      product[1]+=produ[1];
		      product[0]+=produ[0];
		    }
		  //subtract the projection of the prev.
		  //col onto the current
		  //from the current
		  //printf("subtracting projection\n");
		  //printf("prod[1]=%.6f\n",product[1]);
		  //printf("prod[0]=%.6f\n",product[0]);
		  for(int row=0; row<n; row++)
		    {
		      produ=times(product,M[n*row+prev],produ);
		      M[n*row+col][0]-=produ[0];
		      M[n*row+col][1]-=produ[1];
		    }
		}
	    }
	  if(n==3 && col==2) //then we are on the third column 
	    //implement cross product
	    {
	      y=times(star(M[3],c),star(M[7],x),y);
	      z=times(star(M[4],c),star(M[6],x),z);
	      M[2][0]=y[0]-z[0];
	      M[2][1]=y[1]-z[1];
	      y=times(star(M[6],c),star(M[1],x),y);
	      z=times(star(M[0],c),star(M[7],x),z);
	      M[5][0]=y[0]-z[0];
	      M[5][1]=y[1]-z[1];
	      y=times(star(M[0],c),star(M[4],x),y);
	      z=times(star(M[1],c),star(M[3],x),z);
	      M[8][0]=y[0]-z[0];
	      M[8][1]=y[1]-z[1];
	    }
	  ////normalize the resulting column
	  //printf("normalize the result\n");
	  if((n==3 && col!=2) || n!=3)
	    {
	      double sum=getnorm(M,n,col);
	      //divide by sum, aka, the norm
	      //printf("divide by sum=%f\n",sum);
	      for(int row=0; row<n; row++)
		{
		  M[n*row+col][0]/=sum;
		  M[n*row+col][1]/=sum;
		}
	    }
	}
      //gram schmidt done, move on to sending determinant to one
      if(n!=3)
	{
	  //get determinant
	  determinant[0]=0;
	  determinant[1]=0;
	  determinant=getdet(M,n,determinant);
	  //printf("norm of det=%.6f\n", getnorm2(determinant,c));
	  //now that we have the determinant(aka the phase),
	  //we must divide the matrix by the determinant^n
	  double phase=-1.0*atan2(determinant[1],determinant[0])/n;
	  c[0]=cos(phase);
	  c[1]=sin(phase);
	  //make the unitary matrix SU(N) with n=N
	  for(int d=0; d<n*n; d++)
	    {
	      res=times(M[d],c,res);
	      M[d][0]=res[0];
	      M[d][1]=res[1];
	    }
	}
      container[i]=M;
      if(trace_real(M,n)/n>-.8)
	{
	  //printf("trace real=%.6f\n",trace_real(M,n));
	  container[i]=M;
	}
      else
	i--;
    }
  free(produ);
  free(produ1);
  free(res);
  free(c);
  free(x);
  free(y);
  free(z);
  free(product);
  free(determinant);
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

int testarg(int arg, int d, int spacing)
{
  int pow=1;
  for(int k=0; k<d; k++)
    {
      pow*=(spacing);
    } 
  if(arg>=pow || arg<0)
    {
      printf("arg that failed=%d\n", arg);
      exit(EXIT_FAILURE);
    }
  else
    return 0;
}


double**** update(double**** lattice, double*** container, int n, int d, int spacing, double beta)
{
  //int accept=0;
  //int reject=0;
  //storage for coordinates
  int* x=(int*)calloc(d,sizeof(int));
  //storage for a copy of the link we are currently updating
  double **copy=(double**)malloc(n*n*sizeof(double*));
  //storage for intermediate results
  double **res1=(double**)malloc(n*n*sizeof(double*));
  double **res2=(double**)malloc(n*n*sizeof(double*));
  double **newM=(double**)malloc(n*n*sizeof(double*));
  double **newM2=(double**)malloc(n*n*sizeof(double*));
  for(int i=0; i<n*n; i++)
    {
      copy[i]=(double*)calloc(2,sizeof(double));
      newM[i]=(double*)calloc(2,sizeof(double));
      newM2[i]=(double*)calloc(2,sizeof(double));
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
      pow*=(spacing);
    }
  //loop over coordinates
  for(int coor=0; coor<pow; coor++)
    {
      x=findcoord(d, spacing, coor, x);
      for(int mu=0; mu<d; mu++)
	{
	  //increment of the coordinate in the mu direction by one unit
	  int mua=1;
	  for(int i=0; i<mu; i++)
	    {
	      mua*=spacing;
	    }
	  //index of staples
	  int nstap=0;
	  //loop over planes (mu,nu plane)
	  for(int nu=0; nu<d; nu++)
	    {
	      //can't have xx or yy plane, clearly nonsense
	      if(mu==nu)
		continue;
	      //increment of the coordinate in the nu direction by one unit
	      int nua=1;
	      for(int j=0; j<nu; j++)
		{
		  nua*=spacing;
		}
	      //index of staples to store the new staples of the mu,nu plane
	      int stap_start=nstap;
	      //due to periodic boundary conditions we always have 2 staples to compute
	      nstap+=2;
	      //j is our starting index
	      int j=stap_start;
	      //loopy staples part
	      int modup=0;
	      int moddown=0;
	      int sidemod=0;
	      if(x[nu]==spacing-1)
		{
		  moddown=0;
		  modup=-nua*(spacing);
		}
	      if(x[nu]==0)
		{
		  modup=0;
		  moddown=nua*(spacing);
		}
	      if(x[mu]==(spacing-1))
		{
		  sidemod=-(spacing)*mua;
		}
	      //down staple
	      testarg(coor+mua-nua+moddown+sidemod,d,spacing);
	      testarg(coor-nua+moddown,d,spacing);
	      res1=Times(dagger(lattice[coor+mua-nua+moddown+sidemod][nu],newM,n),dagger(lattice[coor-nua+moddown][mu],newM2,n),res1,n);
	      res2=Times(res1,lattice[coor-nua+moddown][nu],res2,n);
	      for(int v=0; v<n*n; v++)
		{
		  staples[j][v][0]=res2[v][0];
		  staples[j][v][1]=res2[v][1];
		}
	      j++;
	      //up staple
	      testarg(coor+nua+modup,d,spacing);
	      testarg(coor+mua+sidemod,d,spacing);
	      res1=Times(dagger(lattice[coor+nua+modup][mu],newM,n),dagger(lattice[coor][nu],newM2,n),res1,n);
	      res2=Times(lattice[coor+mua+sidemod][nu],res1,res2,n);
	      for(int v=0; v<n*n; v++)
		{
		  staples[j][v][0]=res2[v][0];
		  staples[j][v][1]=res2[v][1];
		}
	    }
	
	  //update the link
	  //srand(clock());
	  for(int tr=0; tr<10; tr++)
	    {
	      //select a random SU(n) matrix
	      int ran=(int)((double)rand()/((double)RAND_MAX/200));
	      double ran2=(double)rand()/((double)RAND_MAX);
	      //printf("ran=%d\n",ran);
	      //printf("ran2=%.6f\n",ran2);
	      //form M-1
	      for(int v=0; v<n*n; v++)
		{
		  copy[v][0]=container[ran][v][0]; 
		  copy[v][1]=container[ran][v][1]; 
		  //the -1 part
		  if(v%(n+1)==0)
		    copy[v][0]-=1;
		}
	      //printf("container=\n");
	      //print_mat(container[ran],n);
	      //printf("copy=\n");
	      //print_mat(copy,n);
	      //M-1 * link
	      res2=Times(copy,lattice[coor][mu],res2,n);
	      //calculate dS (METROPOLIS PART)
	      double dS=0;
	      for(int y=0; y<nstap; y++)
		{
		  dS+=-beta*1.0/((double)n*1.0)*trace_real(Times(res2,staples[y],res1,n),n);
		}
	      //metropolis condition
	      //printf("dS=%.6f\n",dS);
	      if(dS<0)
		{
		  //accept++;
		  //printf("accept=%d\n",i);
		  //printf("ds=%.6f\n",dS);
		  //double init=calculate_S(lattice,beta,d,spacing,n);
		  res1=Times(container[ran],lattice[coor][mu],res1,n);
		  for(int k=0; k<n*n; k++)
		    {
		      lattice[coor][mu][k][0]=res1[k][0];
		      lattice[coor][mu][k][1]=res1[k][1];
		    }
		  //printf("ds22222=%.6f\n",calculate_S(lattice,beta,d,spacing,n)-init);
		}
	      else if(dS>0 && exp(-dS)>ran2)
		{
		  //accept++;
		  //printf("exp(-dS)=%.6f\n",exp(-dS));
		  //printf("dS=%.6f\n",dS);
		  //printf("accept=%d\n",i);
		  //double init=calculate_S(lattice,beta,d,spacing,n);
		  res1=Times(container[ran],lattice[coor][mu],res1,n);
		  for(int k=0; k<n*n; k++)
		    {
		      lattice[coor][mu][k][0]=res1[k][0];
		      lattice[coor][mu][k][1]=res1[k][1];
		    }
		  //printf("ds22222=%.6f\n",calculate_S(lattice,beta,d,spacing,n)-init);
		}
	      else
		{
		  //reject++;
		}
	      //end 10 hits of metropolis
	    }
	  //end link loop
	}
      //end coordinate loop
    }
  for(int j=0;j<(d-1)*2;j++)
    {
      for(int i=0; i<n*n; i++)
	{
	  free(staples[j][i]);  
	}
      free(staples[j]);
    }
  for(int i=0; i<n*n; i++)
    {
      free(newM[i]);
      free(newM2[i]);
      free(res1[i]);
      free(res2[i]);
      free(copy[i]);
    }
  free(copy);
  free(res1);
  free(res2);
  free(newM2);
  free(newM);
  free(staples);
  free(x);
  //printf("ratio=(accept,reject)=%d,%d\n",accept,reject);
  return lattice;
}

double calculate_S21(double ****lattice, double beta, int d, int spacing, int n)
{
  double S=0;
  int pow=1;
  for(int k=0; k<d; k++)
    {
      pow*=(spacing);
    }
  for(int coor=0; coor<pow; coor++)
    {
      for(int mu=0; mu<d; mu++)
	{
	  for(int nu=mu+1; nu<d; nu++)
	    {
	      double res=0;
	      res=-beta*plaq21(lattice, coor, mu, nu, n, spacing, d, res);
	      S+=res;
	      if(res/beta>1.01 || -res/beta>1.01)
		{
		  printf("res=%.6f\n",res);
		  printf("fuck\n");
		  exit(EXIT_FAILURE);
		}
	    }
	}
    }
  return S;
}

double calculate_S(double ****lattice, double beta, int d, int spacing, int n)
{
  double S=0;
  int pow=1;
  for(int k=0; k<d; k++)
    {
      pow*=(spacing);
    }
  //printf("pow=%d\n", pow);
  for(int coor=0; coor<pow; coor++)
    {
      for(int mu=0; mu<d; mu++)
	{
	  for(int nu=mu+1; nu<d; nu++)
	    {
	      double res=0;
	      //check1,2,1,2
	      res=-beta*plaq(lattice, coor, nu, mu, n, spacing, d, res);
	      S+=res;
	      if(res/beta>1.01 || -res/beta>1.01)
		{
		  printf("res=%.6f\n",res);
		  printf("fuck\n");
		  exit(EXIT_FAILURE);
		}
	    }
	}
    }
  return S;
}

//coordinate,direction of link, index of nxn matrix, real/imag part
double plaq21(double ****lattice, int coor, int mu, int nu, int n, int spacing,int d, double result=0)
{
  //storage for intermediate results
  int *x=(int*)malloc(d*sizeof(int));
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
  //coordinate
  x=findcoord(d, spacing, coor, x);
  //coordinate increments
  int nua=1;
  int mua=1;
  for(int i=0; i<nu; i++)
    {
      nua*=spacing;
    }
  for(int i=0; i<mu; i++)
    {
      mua*=spacing;
    }
  //periodic boundary conditions means the link that sticks out loops around
  if(x[mu]==spacing-1 && x[nu]!=spacing-1)
    {
      testarg(coor+mua-mua*spacing,d,spacing);
      testarg(coor+nua,d,spacing);
      res1=Times(lattice[coor][mu],lattice[coor+mua-mua*spacing][nu],res1,n);
      res2=Times(dagger(lattice[coor+nua][mu],newM,n),dagger(lattice[coor][nu],newM2, n),res2,n);
    }
  else if(x[mu]!=spacing-1 && x[nu]==spacing-1)
    {
      testarg(coor+mua,d,spacing);
      testarg(coor+nua-nua*spacing,d,spacing);
      res1=Times(lattice[coor][mu],lattice[coor+mua][nu],res1,n);
      res2=Times(dagger(lattice[coor+nua-nua*spacing][mu],newM,n),dagger(lattice[coor][nu],newM2, n),res2,n);
    }
  else if(x[mu]==spacing-1 && x[nu]==spacing-1)
    {
      testarg(coor+nua-nua*spacing,d,spacing);
      testarg(coor+mua-mua*spacing,d,spacing);
      res1=Times(lattice[coor][mu],lattice[coor+mua-mua*spacing][nu],res1,n);
      res2=Times(dagger(lattice[coor+nua-nua*spacing][mu],newM,n),dagger(lattice[coor][nu],newM2, n),res2,n);
    }
  else
    {
      testarg(coor+nua,d,spacing);
      testarg(coor+mua,d,spacing);
      res1=Times(lattice[coor][mu],lattice[coor+mua][nu],res1,n);
      res2=Times(dagger(lattice[coor+nua][mu],newM,n),dagger(lattice[coor][nu],newM2, n),res2,n);
    }
  //result
  result=1.0/((double)n*1.0)*trace_real(Times(res1,res2,res3,n),n);
  for(int i=0; i<n*n; i++)
    {
      free(newM[i]);
      free(newM2[i]);
      free(res1[i]);
      free(res2[i]);
      free(res3[i]);
    }
  free(newM);
  free(newM2);
  free(res1);
  free(res2);
  free(res3);
  free(x);
  return result;
}

//coordinate,direction of link, index of nxn matrix, real/imag part
double plaq(double ****lattice, int coor, int mu, int nu, int n, int spacing,int d, double result=0)
{
  //storage for intermediate results
  int *x=(int*)malloc(d*sizeof(int));
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
  //coordinate
  x=findcoord(d, spacing, coor, x);
  //coordinate increments
  int nua=1;
  int mua=1;
  for(int i=0; i<nu; i++)
    {
      nua*=spacing;
    }
  for(int i=0; i<mu; i++)
    {
      mua*=spacing;
    }
  //periodic boundary conditions means the link that sticks out loops around
  if(x[mu]==spacing-1 && x[nu]!=spacing-1)
    {
      testarg(coor+mua-mua*spacing,d,spacing);
      testarg(coor+nua,d,spacing);
      res1=Times(lattice[coor][mu],lattice[coor+mua-mua*spacing][nu],res1,n);
      res2=Times(dagger(lattice[coor+nua][mu],newM,n),dagger(lattice[coor][nu],newM2, n),res2,n);
    }
  else if(x[mu]!=spacing-1 && x[nu]==spacing-1)
    {
      testarg(coor+mua,d,spacing);
      testarg(coor+nua-nua*spacing,d,spacing);
      res1=Times(lattice[coor][mu],lattice[coor+mua][nu],res1,n);
      res2=Times(dagger(lattice[coor+nua-nua*spacing][mu],newM,n),dagger(lattice[coor][nu],newM2, n),res2,n);
    }
  else if(x[mu]==spacing-1 && x[nu]==spacing-1)
    {
      testarg(coor+nua-nua*spacing,d,spacing);
      testarg(coor+mua-mua*spacing,d,spacing);
      res1=Times(lattice[coor][mu],lattice[coor+mua-mua*spacing][nu],res1,n);
      res2=Times(dagger(lattice[coor+nua-nua*spacing][mu],newM,n),dagger(lattice[coor][nu],newM2, n),res2,n);
    }
  else
    {
      testarg(coor+nua,d,spacing);
      testarg(coor+mua,d,spacing);
      res1=Times(lattice[coor][mu],lattice[coor+mua][nu],res1,n);
      res2=Times(dagger(lattice[coor+nua][mu],newM,n),dagger(lattice[coor][nu],newM2, n),res2,n);
    }
  //result
  result=1.0/((double)n*1.0)*trace_real(Times(res1,res2,res3,n),n);
  for(int i=0; i<n*n; i++)
    {
      free(newM[i]);
      free(newM2[i]);
      free(res1[i]);
      free(res2[i]);
      free(res3[i]);
    }
  free(newM);
  free(newM2);
  free(res1);
  free(res2);
  free(res3);
  free(x);
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

double trace_img(double **M, int n)
{
  double sum=0;
  for(int k=0; k<n; k++)
    {
      sum+=M[k*(n+1)][1];
    }
  return sum;
}

double trace_real(double **M, int n)
{
  double sum=0;
  for(int k=0; k<n; k++)
    {
      sum+=M[k*(n+1)][0];
    }
  return sum;
}

double** dagger(double** M, double** newM, int n)
{
  for(int k=0; k<n*n; k++)
    {
      int i=k/n;
      int j=k-i*n;
      int kprime=j*n+i;
      if(kprime==k)
	{
	  newM[k][0]=M[k][0];
	  newM[k][1]=-1*M[k][1];
	}
      else if(kprime>k)
	{
	  newM[k][0]=M[kprime][0];
	  newM[k][1]=-1*M[kprime][1];
	  newM[kprime][0]=M[k][0];
	  newM[kprime][1]=-1*M[k][1];
	}
    }
  return newM;
}

double**** initialize_lat(int d, int spacing, int n)
{
  int pow=1;
  for(int k=0; k<d; k++)
    {
      pow*=spacing;
    }
  double**** lattice=(double****)malloc(pow*sizeof(double***));
  //coordinate,direction of link, index of nxn matrix, real/imag part
  int* x=(int*)malloc(d*sizeof(int));
  for(int coor=0; coor<pow; coor++)
    {
      x=findcoord(d, spacing, coor, x);
      lattice[coor]=(double***)malloc(d*sizeof(double**));
      for(int dir=0; dir<d; dir++)
	{
	  lattice[coor][dir]=(double**)malloc(n*n*sizeof(double*));
	  for(int index=0; index<n*n; index++)
	    {
	      lattice[coor][dir][index]=(double*)malloc(2*sizeof(double));
	      //below: first condition sees if we are on diagonal
	      //second condition tests if we are on the border of the lattice
	      //in which case we create no dangling links
	      if(index%(n+1)==0 && x[dir]!=(spacing))
		{
		  lattice[coor][dir][index][0]=1;
		  lattice[coor][dir][index][1]=0;
		  if((index==0 || index==n+1)&& n%2!=0 && dir==0)
		    {
		      //slightly disordered starting point
		      //lattice[coor][dir][index][0]=-1;
		      //ordered
		      lattice[coor][dir][index][0]=1;
		    }
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
      x[temp]=coors-(coors/(spacing))*(spacing);
      coors/=(spacing);
      temp++;
    }
  return x;
}

int main()
{
  srand(time(NULL));
  ////global constants
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
  double ***container=gen_rand_matrix(n,epsilon);
  //return 0;
  //place the daggers of the random SU(3) matrices in the container as well
  container=(double***)realloc(container,2*100*sizeof(double**));
  for(int i=100; i<200; i++)
    {
      container[i]=(double**)malloc(n*n*sizeof(double*));
      newM=dagger(container[i-100],newM,n);
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
  double**** lattice=initialize_lat(d,spacing,n);
  //calculate the total number of plaquettes
  double s2=calculate_S(lattice, beta, d, spacing,n);
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
      fprintf(pFile,"axa avg=%.6f\n",calculate_S(lattice,-1,d,spacing,n)/total_plqts);
      printf("axa avg=%.6f\n",calculate_S(lattice,-1,d,spacing,n)/total_plqts);
      lattice=update(lattice,container,n,d,spacing,beta);
    }


  //calculate the average plaquette, store it
  //double avg_plqt=0;
  double avg_plqt1=0;
  //update the links Ncor times, save S, repeat Ncf-1 times
  for(int i=0; i<Ncf; i++)
    {
      avg_plqt1=calculate_S(lattice, -1, d, spacing,n)/total_plqts;
      //avg_plqt+=avg_plqt/Ncf/total_plqts;
      fprintf(pFile,"avg_plqt(axa)(path number)=(%.6f)(%d)\n",avg_plqt1,i);
      printf("avg_plqt(axa)(path number)=(%.6f)(%d)\n",avg_plqt1,i);
      for(int j=0; j<Ncor; j++)
	{
	  lattice=update(lattice,container,n,d,spacing,beta);
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
