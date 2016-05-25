double** minors(double **M, int nn, int row, int column, double** minorr)
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

double* getdet(double** M, int nn, double* temp_fin)
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

double getnorm(double** M, int col=0)
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
double*** gen_rand_matrix()
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
	      double sum=getnorm(M,col);
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
	  for(int r=0; r<n*n; r++)
	    {
	      res=times(M[r],c,res);
	      M[r][0]=res[0];
	      M[r][1]=res[1];
	    }
	}
      container[i]=M;
      if(trace_real(M)/n>-.8)
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
