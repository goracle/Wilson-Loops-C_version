int testarg(int arg)
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

double getnorm2(double* a, double *result)
{
  return times(star(a,result),a,result)[0];
}

int increment(int mu)
{
  int res=1;
  for(int r=0; r<mu; r++)
    {
      res*=spacing;
    }
  return res;
}
int* findcoord(int coor, int* x)
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


int* mods(int coor, int mu, int nu, int* results)
{
  int* x=(int*)malloc(d*sizeof(int)); 
  x=findcoord(coor,x);
  int mua=increment(mu);
  int nua=increment(nu);
  //mnemonic=Konami code
  //upmod
  if(x[nu]==spacing-1)
    {
      results[0]=-nua*spacing;
    }
  else
    results[0]=0;
  //downmod
  if(x[nu]==0)
    {
      results[1]=nua*spacing; 
    }
  else
    results[1]=0;
  //sidemod
  if(x[mu]==spacing-1)
    {
      results[2]=-mua*spacing; 
    }
  else
    results[2]=0;
  //b a...
  free(x);
  return results;
}
