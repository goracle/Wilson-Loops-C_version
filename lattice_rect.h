//update with the improved action
double**** update_rect(double**** lattice, double*** container)
{
  //int accept=0;
  //int reject=0;
  //storage for coordinates
  int* x=(int*)calloc(d,sizeof(int));
  //storage for a copy of the link we are currently updating
  double **copy=(double**)malloc(n*n*sizeof(double*));
  //storage for intermediate results
  double **res1=(double**)malloc(n*n*sizeof(double*));
  double **res11=(double**)malloc(n*n*sizeof(double*));
  double **res2=(double**)malloc(n*n*sizeof(double*));
  double **newM=(double**)malloc(n*n*sizeof(double*));
  double **newM2=(double**)malloc(n*n*sizeof(double*));
  //storage for findstap_rect test
  /*  //mods for periodic boundary conditions
   int* modsv=(int*)malloc(3*sizeof(int));
  //misc stack vars
  double resr=0;
  double diff=0;*/
  for(int i=0; i<n*n; i++)
    {
      copy[i]=(double*)calloc(2,sizeof(double));
      newM[i]=(double*)calloc(2,sizeof(double));
      newM2[i]=(double*)calloc(2,sizeof(double));
      res1[i]=(double*)calloc(2,sizeof(double));
      res11[i]=(double*)calloc(2,sizeof(double));
      res2[i]=(double*)calloc(2,sizeof(double));
    }
  //storage for staples
  double ***staples=(double***)malloc((d-1)*2*sizeof(double**));
  double ***staples1=(double***)malloc((d-1)*6*sizeof(double**));
  for(int j=0; j<(d-1)*6; j++)
    {
      staples1[j]=(double**)malloc(n*n*sizeof(double*));
      for(int i=0; i<n*n; i++)
	{
	  staples1[j][i]=(double*)calloc(2,sizeof(double));
	}
    }
  for(int j=0; j<(d-1)*2; j++)
    {
      staples[j]=(double**)malloc(n*n*sizeof(double*));
      for(int i=0; i<n*n; i++)
	{
	staples[j][i]=(double*)calloc(2,sizeof(double));
	}
    }
  //plow=number of coordinates in the lattice
  int plow=1;
  for(int k=0; k<d; k++)
    {
      plow*=(spacing);
    }
  //loop over coordinates
  for(int coor=0; coor<plow; coor++)
    {
      x=findcoord(coor, x);
      for(int mu=0; mu<d; mu++)
	{
	  //increment of the coordinate in the mu direction by one unit
	  int mua=increment(mu);
	  //index of staples
	  int nstap=0;
	  int nstap1=0;
	  //loop over planes (mu,nu plane)
	  for(int nu=0; nu<d; nu++)
	    {
	      //can't have xx or yy plane, clearly nonsense
	      if(mu==nu)
		continue;
	      //modsv=mods(coor,mu,nu,modsv);
	      int stap_start=nstap;
	      int stap_start1=nstap1;
	      //due to periodic boundary conditions
	      //we always have 6 staples to compute
	      nstap1+=6;
	      nstap+=2;
	      //index of staples to store the new staples of the mu,nu plane
	      //j is our starting index
	      int nua=increment(nu);
	      for(int j=stap_start1; j<nstap1; j++)
		{
		  //printf("j=%d, stap_start1=%d, nstap1=%d \n",j,stap_start1,nstap1);
		  staples1=findstap_rect(mu,mua,nu,x,coor,j,res1,res11,res2,newM,newM2,lattice,staples1,j-stap_start1);
		  //test for rect staples
		  /*		  if(j-stap_start1==0)
		    {
		      diff=1.0/((double)n)*trace_real(Times(staples1[j],lattice[coor][mu],res1))-plaq_rect(lattice,coor-nua+modsv[1],nu,mu,resr);
		    }
		  if(j-stap_start1==1)
		    {
		      diff=1.0/((double)n)*trace_real(Times(staples1[j],lattice[coor][mu],res1))-plaq_rect(lattice,coor+mua+modsv[2]-nua+modsv[1],nu,mu,resr);
		    }
		  if(j-stap_start1==2)
		    {
		      diff=1.0/((double)n)*trace_real(Times(staples1[j],lattice[coor][mu],res1))-plaq_rect(lattice,coor-nua+modsv[1],mu,nu,resr);
		    }
		  if(j-stap_start1==3)
		    {
		      diff=1.0/((double)n)*trace_real(Times(staples1[j],lattice[coor][mu],res1))-plaq_rect(lattice,coor,nu,mu,resr);
		    }
		  if(j-stap_start1==4)
		    {
		      diff=1.0/((double)n)*trace_real(Times(staples1[j],lattice[coor][mu],res1))-plaq_rect(lattice,coor+mua+modsv[2],nu,mu,resr);
		    }
		  if(j-stap_start1==5)
		    {
		      diff=1.0/((double)n)*trace_real(Times(staples1[j],lattice[coor][mu],res1))-plaq_rect(lattice,coor+nua+modsv[0],mu,nu,resr);
		    }
		  if((diff>0.0000000001 || diff<-0.0000000001) && (j-stap_start1!=4))
			{
			  printf("diff=%.6f, staple=%d\n", diff,j-stap_start1);
			  exit(EXIT_FAILURE);
			}
		  */
		}
	      for(int j=stap_start; j<nstap; j++)
		{
		  //printf("j=%d, stap_start=%d, nstap=%d \n",j,stap_start,nstap);
		  staples=findstap(mu,mua,nu,x,coor,j,res1,res2,newM,newM2,lattice,staples,j-stap_start);
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
	      res2=Times(copy,lattice[coor][mu],res2);
	      //calculate dS (METROPOLIS PART)
	      double dS=0;
	      for(int y=0; y<nstap1; y++)
		{
		 dS+=beta*1.0/((double)n*1.0)/((double)pow((double)u0,6.0))/12.0*trace_real(Times(res2,staples1[y],res1));
		}
	      for(int y=0; y<nstap; y++)
		{
		  dS+=-beta*5.0/((double)n*1.0)/(pow(u0,4.0))/3.0*trace_real(Times(res2,staples[y],res1));
		}
	      //metropolis condition
	      //printf("dS=%.6f\n",dS);
	      if(dS<0)
		{
		  //printf("dS,from staples=%.6f\n",dS);
		  //double init=calculate_S_rect(lattice);
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
		  //printf("dS,from S_calc=%.6f\n",calculate_S_rect(lattice)-init);
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
  for(int j=0; j<(d-1)*6; j++)
    {
      for(int i=0; i<n*n; i++)
	{
	  free(staples1[j][i]);
	}
      free(staples1[j]);
    }
  for(int i=0; i<n*n; i++)
    {
      free(newM[i]);
      free(newM2[i]);
      free(res1[i]);
      free(res11[i]);
      free(res2[i]);
      free(copy[i]);
    }
  // free(modsv);
  free(copy);
  free(res1);
  free(res11);
  free(res2);
  free(newM2);
  free(newM);
  free(staples);
  free(staples1);
  free(x);
  //printf("ratio=(accept,reject)=%d,%d\n",accept,reject);
  //exit(EXIT_SUCCESS);
  return lattice;
}


double*** findstap_rect(int mu, int mua, int nu, int* x, int coor, int j, double** res1, double** res11, double** res2, double** newM, double** newM2, double**** lattice, double ***staples, int updown)
{
  int* y=(int*)malloc(d*sizeof(int));
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
  //down staples
  if(updown==0)
    {
      testarg(coor+mua-nua+modsv[2]+modsv[1]);
      testarg(coor-nua+modsv[1]);
      res1=Times(dagger(lattice[coor+mua-nua+modsv[1]+modsv[2]][nu],newM),dagger(lattice[coor-nua+modsv[1]][mu],newM2),res1);
      y=findcoord(coor-nua+modsv[1],y);
      staples=findstap(nu,nua,mu,y,coor-nua+modsv[1],j,res11,res2,newM,newM2,lattice,staples,0);
      res2=Times(res1,dagger(staples[j],newM2),res2);
    }
  else if(updown==1)
    {
      testarg(coor+mua-nua+modsv[1]+modsv[2]);
      testarg(coor-nua+modsv[1]);
      //end testargs
      res1=Times(dagger(lattice[coor-nua+modsv[1]][mu],newM),lattice[coor-nua+modsv[1]][nu],res1);
      y=findcoord(coor-nua+modsv[1]+mua+modsv[2],y);
      staples=findstap(nu,nua,mu,y,coor-nua+mua+modsv[1]+modsv[2],j,res11,res2,newM,newM2,lattice,staples,1);
      res2=Times(staples[j],res1,res2);
    }
  else if(updown==2)
    {
      testarg(coor+mua-nua+modsv[1]+modsv[2]);
      testarg(coor-nua+modsv[1]);
      y=findcoord(coor-nua+modsv[1],y);
      staples=findstap(mu,mua,nu,y,coor-nua+modsv[1],j,res11,res2,newM,newM2,lattice,staples,0);
      res1=Times(dagger(lattice[coor+mua-nua+modsv[2]+modsv[1]][nu],newM2),staples[j],res1);
      res2=Times(res1,lattice[coor-nua+modsv[1]][nu],res2);
    }
  //up staples
  else if(updown==3)
    {
      testarg(coor+nua+modsv[0]);
      testarg(coor+mua+modsv[2]);
      y=findcoord(coor,y);
      staples=findstap(nu,nua,mu,y,coor,j,res1,res2,newM,newM2,lattice,staples,0);
      res1=Times(dagger(lattice[coor+nua+modsv[0]][mu],newM),staples[j],res1);
      res2=Times(lattice[coor+mua+modsv[2]][nu],res1,res2);
    }
 else if(updown==4)
   {
      testarg(coor+nua+modsv[0]);
      testarg(coor+mua+modsv[2]);
      res1=Times(dagger(lattice[coor+nua+modsv[0]][mu],newM),dagger(lattice[coor][nu],newM2),res1);
      y=findcoord(coor+mua+modsv[2],y);
      staples=findstap(nu,nua,mu,y,coor+mua+modsv[2],j,res11,res2,newM,newM2,lattice,staples,1);
      res2=Times(dagger(staples[j],newM),res1,res2);
   }
 else if(updown==5)
   {
      testarg(coor+nua+modsv[0]);
      testarg(coor+mua+modsv[2]);
      y=findcoord(coor+nua+modsv[0],y);
      staples=findstap(mu,mua,nu,y,coor+nua+modsv[0],j,res1,res2,newM,newM2,lattice,staples,1);
      res1=Times(lattice[coor+mua+modsv[2]][nu],staples[j],res1);
      res2=Times(res1,dagger(lattice[coor][nu],newM),res2);
   }
  else
    {
      printf("nonsense argument to findstap_rect(...,updown=%d)\n",updown);
      exit(EXIT_FAILURE);
    }
  //store the result
  for(int v=0; v<n*n; v++)
    {
      staples[j][v][0]=res2[v][0];
      staples[j][v][1]=res2[v][1];
    }
  free(y);
  free(modsv);
  return staples;
}

double calculate_S_rect(double ****lattice)
{
  double S=0;
  int plow=1;
  double res=0;
  for(int k=0; k<d; k++)
    {
      plow*=(spacing);
    }
  //printf("plow=%d\n", plow);
  for(int coor=0; coor<plow; coor++)
    {
      for(int mu=0; mu<d; mu++)
	{
	  for(int nu=mu+1; nu<d; nu++)
	    {
	      //check1,2,1,2
	      S+=beta/pow(u0,6.0)/12.0*plaq_rect(lattice, coor, nu, mu, res);
	      S+=beta/pow(u0,6.0)/12.0*plaq_rect(lattice, coor, mu, nu, res);
	      S+=-5.0*beta/pow(u0,4.0)/3.0*plaq(lattice, coor, nu, mu, res);
	      /*if(res/beta>1.01 || -res/beta>1.01)
		{
		  printf("res=%.6f\n",res);
		  printf("fuck\n");
		  exit(EXIT_FAILURE);
		}
	      */
	    }
	}
    }
  return S;
}

double plaq_rect(double ****lattice, int coor, int mu, int nu, double result)
{
  //storage for intermediate results
  int mua=increment(mu);
  int *x=(int*)malloc(d*sizeof(int));
  x=findcoord(coor, x);
  double ***staples=(double***)malloc(2*sizeof(double**));
  staples[0]=(double**)malloc(n*n*sizeof(double*));
  staples[1]=(double**)malloc(n*n*sizeof(double*));
  double **res1=(double**)malloc(n*n*sizeof(double*));
  double **res2=(double**)malloc(n*n*sizeof(double*));
  double **res3=(double**)malloc(n*n*sizeof(double*));
  double **newM=(double**)malloc(n*n*sizeof(double*));
  double **newM2=(double**)malloc(n*n*sizeof(double*));
  for(int i=0; i<n*n; i++)
    {
      staples[0][i]=(double*)calloc(2,sizeof(double));
      staples[1][i]=(double*)calloc(2,sizeof(double));
      newM[i]=(double*)calloc(2,sizeof(double));
      newM2[i]=(double*)calloc(2,sizeof(double));
      res1[i]=(double*)calloc(2,sizeof(double));
      res2[i]=(double*)calloc(2,sizeof(double));
      res3[i]=(double*)calloc(2,sizeof(double));
    }
  findstap(mu,mua,nu,x,coor,0,res1,res2,newM,newM2,lattice,staples,0);
  findstap(mu,mua,nu,x,coor,1,res1,res2,newM,newM2,lattice,staples,1);
  newM=dagger(staples[0],newM);
  result=1.0/((double)n*1.0)*trace_real(Times(staples[1],newM,res2));
  for(int i=0; i<n*n; i++)
    {
      free(staples[0][i]);
      free(staples[1][i]);
      free(newM[i]);
      free(newM2[i]);
      free(res1[i]);
      free(res2[i]);
      free(res3[i]);
    }
  free(newM);
  free(staples[0]);
  free(staples[1]);
  free(staples);
  free(newM2);
  free(res1);
  free(res2);
  free(res3);
  free(x);
  return result;
}

//unused code below
//coordinate,direction of link, index of nxn matrix, real/imag part
/*
double plaq_rect1(double ****lattice, int coor, int mu, int nu, double result)
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
  int sidemod2=0;
  int sidemod1=0;
  int sidemod=0;
  int modup1=0;
  int modup=0;
  //sidemod
  if(x[mu]==spacing-1)
    {
      sidemod=-spacing*mua;
    }
  if(x[mu]==spacing-2)
    {
      sidemod1=-spacing*mua;
    }
  if(x[nu]==spacing-1)
    {
      modup=-spacing*nua;
    }
  if(x[nu]==spacing-1)
    {
      modup1=-spacing*nua;
    }
  if(x[mu]==0 && x[mu]==spacing-1)
    {
      sidemod2=-spacing*mua;
      printf("spacing=%d\n",spacing);
      exit(EXIT_FAILURE);
    }
  //periodic boundary conditions means the link that sticks out loops around
  res1=Times(lattice[coor][mu],lattice[coor+mua+sidemod][mu],res1);
  res2=Times(lattice[coor+2*mua+sidemod1+sidemod2+sidemod][nu],dagger(lattice[coor+mua+nua+modup+sidemod][nu],newM),res2);
  res3=Times(res1,res2,res3);
  res1=Times(dagger(lattice[coor+nua+modup][mu],newM),dagger(lattice[coor][nu],newM2),res1);
  //result
  result=1.0/((double)n*1.0)*trace_real(Times(res3,res1,res2));
  //print_mat(res2);
  //printf("result=%.6f\n",result);
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
*/
