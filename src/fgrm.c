#include "common.h"

static flouble freq_evolve(spec_type,nu_0,beta,temp,nu)
{
  flouble x_to,x_from,ex;
  switch(spec_type)  {
  case 0 :
    x_to=0.017611907*nu;
    ex=exp(x);
    x_to=x_to/(ex-1);
    return ex*x_to*x_to;
    break;
  case 1 :
    return pow(nu/nu_0,beta-2.);
    break;
  case 2 :
    x_to=0.0479924466*nu/temp; //DAM: possible optimization, use 1/T instead of T
    x_from=0.0479924466*nu_0/temp;
    return pow(nu/nu_0,beta+1.)*(exp(x_from)-1)/(exp(x_to)-1);
  }
}

