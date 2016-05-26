double avg_plqt_rect(double**** lattice)
{
  double result=0; 
  double res=0;
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
	  for(int nu=0; nu<d; nu++)
	    {
	      if(mu==nu)
		continue;
	      res=plaq_rect(lattice, coor, mu, nu, res);
	      result+=res;
	    }
	}
    }
  return 1.0*result/((double)pow)/((double)d)/((double)(d-1));
}

double avg_plqt(double**** lattice)
{
  double result=0; 
  int pow=1;
  double res=0;
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
	      res=plaq(lattice, coor, nu, mu, res);
	      result+=res;
	    }
	}
    }
  return 2*result/((double)pow)/((double)d)/((double)d-1.0);
}
