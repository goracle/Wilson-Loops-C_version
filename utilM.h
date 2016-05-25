double* times(double* y, double* z,double* result)
{
  result[0]=z[0]*y[0]-z[1]*y[1];
  result[1]=z[0]*y[1]+z[1]*y[0];
  return result;
}

double* star(double* x, double* result)
{
  result[0]=x[0];
  result[1]=-x[1];
  return result;
}


int print_mat(double** M)
{
  if(n!=3)
    return -1;
  printf("\nmarker\n");
  printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", M[0][0],M[0][1],M[1][0],M[1][1],M[2][0],M[2][1]);
  printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", M[3][0],M[3][1],M[4][0],M[4][1],M[5][0],M[5][1]);
  printf("%.6f + i %.6f, %.6f + i %.6f, %.6f + i %.6f \n", M[6][0],M[6][1],M[7][0],M[7][1],M[8][0],M[8][1]);
  return 0;
}
double trace_img(double **M)
{
  double sum=0;
  for(int k=0; k<n; k++)
    {
      sum+=M[k*(n+1)][1];
    }
  return sum;
}

double trace_real(double **M)
{
  double sum=0;
  for(int k=0; k<n; k++)
    {
      sum+=M[k*(n+1)][0];
    }
  return sum;
}

double** dagger(double** M, double** newM)
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
double **Times(double** M, double** N, double** res)
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

