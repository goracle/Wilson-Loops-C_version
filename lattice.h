

double**** update(double**** lattice, double*** container)
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
      x=findcoord(coor, x);
      for(int mu=0; mu<d; mu++)
	{
	  //increment of the coordinate in the mu direction by one unit
	  int mua=increment(mu);
	  //index of staples
	  int nstap=0;
	  //loop over planes (mu,nu plane)
	  for(int nu=0; nu<d; nu++)
	    {
	      //can't have xx or yy plane, clearly nonsense
	      if(mu==nu)
		continue;
	      int stap_start=nstap;
	      //due to periodic boundary conditions
	      //we always have 2 staples to compute
	      nstap+=2;
	      //index of staples to store the new staples of the mu,nu plane
	      //j is our starting index
	      int j=stap_start;
	      staples=findstap(mu,mua,nu,x,coor,j,res1,res2,newM,newM2,lattice,staples,0);
	      staples=findstap(mu,mua,nu,x,coor,j+1,res1,res2,newM,newM2,lattice,staples,1);
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
	      res2=Times(copy,lattice[coor][mu],res2);
	      //calculate dS (METROPOLIS PART)
	      double dS=0;
	      for(int y=0; y<nstap; y++)
		{
		  dS+=-beta*1.0/((double)n*1.0)*trace_real(Times(res2,staples[y],res1));
		}
	      //metropolis condition
	      //printf("dS=%.6f\n",dS);
	      if(dS<0)
		{
		  //accept++;
		  //printf("accept=%d\n",i);
		  //printf("ds=%.6f\n",dS);
		  //double init=calculate_S(lattice,beta,d,spacing,n);
		  res1=Times(container[ran],lattice[coor][mu],res1);
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
		  res1=Times(container[ran],lattice[coor][mu],res1);
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

double*** findstap(int mu, int mua, int nu, int* x, int coor, int j, double** res1, double** res2, double** newM, double** newM2, double**** lattice, double ***staples, int updown)
{
  if(mu==nu)
    {
      printf("mu==nu;  fuck\n");
      exit(EXIT_FAILURE);
    }
  //increment of the coordinate in the nu direction by one unit
  int nua=increment(nu);
  //loopy (periodic bound. cond.) staples part
  int* modsv=(int*)malloc(3*sizeof(int));
  modsv=mods(coor,mu,nu,modsv);
  //down staple
  if(updown==0)
    {
      testarg(coor+mua-nua+modsv[1]+modsv[2]);
      testarg(coor-nua+modsv[1]);
      res1=Times(dagger(lattice[coor+mua-nua+modsv[1]+modsv[2]][nu],newM),dagger(lattice[coor-nua+modsv[1]][mu],newM2),res1);
      res2=Times(res1,lattice[coor-nua+modsv[1]][nu],res2);
    }
  //up staple
  else if(updown==1)
    {
      testarg(coor+nua+modsv[0]);
      testarg(coor+mua+modsv[2]);
      res1=Times(dagger(lattice[coor+nua+modsv[0]][mu],newM),dagger(lattice[coor][nu],newM2),res1);
      res2=Times(lattice[coor+mua+modsv[2]][nu],res1,res2);
    }
  else
    {
      printf("nonsense argument to findstap(...,updown=%d)\n",updown);
      exit(EXIT_FAILURE);
    }
  //store the result
  for(int v=0; v<n*n; v++)
    {
      staples[j][v][0]=res2[v][0];
      staples[j][v][1]=res2[v][1];
    }
  free(modsv);
  return staples;
}



double calculate_S(double ****lattice)
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
	      res=-beta*plaq(lattice, coor, nu, mu, res);
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
double plaq(double ****lattice, int coor, int mu, int nu, double result=0)
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
  x=findcoord(coor, x);
  //coordinate increments
  int nua=increment(nu);
  int mua=increment(mu);
  //periodic boundary conditions means the link that sticks out loops around
  if(x[mu]==spacing-1 && x[nu]!=spacing-1)
    {
      testarg(coor+mua-mua*spacing);
      testarg(coor+nua);
      res1=Times(lattice[coor][mu],lattice[coor+mua-mua*spacing][nu],res1);
      res2=Times(dagger(lattice[coor+nua][mu],newM),dagger(lattice[coor][nu],newM2),res2);
    }
  else if(x[mu]!=spacing-1 && x[nu]==spacing-1)
    {
      testarg(coor+mua);
      testarg(coor+nua-nua*spacing);
      res1=Times(lattice[coor][mu],lattice[coor+mua][nu],res1);
      res2=Times(dagger(lattice[coor+nua-nua*spacing][mu],newM),dagger(lattice[coor][nu],newM2),res2);
    }
  else if(x[mu]==spacing-1 && x[nu]==spacing-1)
    {
      testarg(coor+nua-nua*spacing);
      testarg(coor+mua-mua*spacing);
      res1=Times(lattice[coor][mu],lattice[coor+mua-mua*spacing][nu],res1);
      res2=Times(dagger(lattice[coor+nua-nua*spacing][mu],newM),dagger(lattice[coor][nu],newM2),res2);
    }
  else
    {
      testarg(coor+nua);
      testarg(coor+mua);
      res1=Times(lattice[coor][mu],lattice[coor+mua][nu],res1);
      res2=Times(dagger(lattice[coor+nua][mu],newM),dagger(lattice[coor][nu],newM2),res2);
    }
  //result
  result=1.0/((double)n*1.0)*trace_real(Times(res1,res2,res3));
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


double**** initialize_lat()
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
      x=findcoord(coor, x);
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
