#include "common.h"

typedef struct {
  flouble *f_matrix;
  flouble *cov_inv;
  flouble *vec_mean;
} PixelState;

static flouble freq_evolve(int spec_type,double nu_0,double beta,double temp,double nu)
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

static void compute_f_matrix(ParamFGRM *par,flouble *x_spec,flouble *f_matrix)
{
  int inu;
  for(inu=0;inu<par->n_nu;inu++) {
    int ipol;
    flouble nu=par->freqs[inu];
    if(par->flag_include_cmb) {
      for(ipol=0;ipol<par->n_pol;ipol++) {
	f_matrix[par->index_cmb+par->n_comp*(inu+ipol*par->n_nu)]=
	  freq_evolve(0,-1,-1,-1,nu);
      }
    }
    if(par->flag_include_synchrotron) {
      f_matrix[par->index_synchrotron+par->n_comp*(inu+0*par->n_nu)]=
	freq_evolve(1,par->nu0_s,x_spec[par->index_beta_s_t],-1,nu);
      for(ipol=1;ipol<par->n_pol;ipol++) {
	f_matrix[par->index_synchrotron+par->n_comp*(inu+ipol*par->n_nu)]=
	  freq_evolve(1,par->nu0_s,x_spec[par->index_beta_s_p],-1,nu);
      }
    }
    if(par->flag_include_dust) {
      f_matrix[par->index_dust+par->n_comp*(inu+0*par->n_nu)]=
	freq_evolve(2,par->nu0_d,x_spec[par->index_beta_d_t],
		    x_spec[par->index_temp_d_t],nu);
      for(ipol=1;ipol<par->n_pol;ipol++) {
	f_matrix[par->index_dust+par->n_comp*(inu+ipol*par->n_nu)]=
	  freq_evolve(2,par->nu0_d,x_spec[par->index_beta_d_p],
		      x_spec[par->index_temp_d_p],nu);
      }
    }
  }
}

static flouble compute_chi2(ParamFGRM *par,flouble *data,flouble *noise_w,
			    flouble *amps,flouble *x_spec,PixelState *pst)
{
  int ipix;
  flouble chi2=0;

  compute_f_matrix(par,x_spec,pst->f_matrix);

  for(ipix=0;ipix<par->n_sub;ipix++) {
    int ipol;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int inu;
      for(inu=0;inu<par->n_nu;inu++) {
	int icomp;
	flouble res=data[inu+par->n_nu*(ipol+par->n_pol*ipix)];
	for(icomp=0;icomp<par->n_comp;icomp++) {
	  res-=pst->f_matrix[icomp+par->n_comp*(inu+par->n_nu*ipol)]*
	    amps[icomp+par->n_comp*(ipol+par->n_pol*ipix)];
	}
	chi2+=res*res*noise_w[inu+par->n_nu*(ipol+par->n_pol*ipix)];
      }
    }
  }

  return chi2;
}

static void analyze_linear_chi2(ParamFGRM *par,flouble *data,flouble *noise_w,
				flouble *x_spec,PixelState *pst)
{
  int ipix;

  compute_f_matrix(par,x_spec,pst->f_matrix);

  for(ipix=0;ipix<par->n_sub;ipix++) {
    int ipol;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int ic1;
      int index_pix=ipol+par->n_pol*ipix;
      for(ic1=0;ic1<par->n_comp;ic1++) {
	int ic2;
	for(ic2=0;ic2<par->n_comp;ic2++) {
	  int inu;
	  flouble icov=0;
	  for(inu=0;inu<par->n_nu;inu++) {
	    int i_lo=inu+par->n_nu*ipol;
	    icov+=pst->f_matrix[ic1+par->n_comp*i_lo]*pst->f_matrix[ic2+par->n_comp*i_lo]*
	      noise_w[inu+par->n_nu*index_pix];
	  }
	  pst->cov_inv[ic2+par->n_comp*(ic1+par->n_comp*index_pix)]=icov;
	}
      }
    }
  }

  for(ipix=0;ipix<par->n_sub;ipix++) {
    int ipol;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int ic;
      int index_pix=ipol+par->n_pol*ipix;
      for(ic=0;ic<par->n_comp;ic++) {
	int inu;
	flouble vec=0;
	for(inu=0;inu<par->n_nu;inu++) {
	  int i_hi=inu+par->n_nu*index_pix;
	  vec+=pst->f_matrix[ic+par->n_comp*(inu+par->n_nu*ipol)]*data[i_hi]*noise_w[i_hi];
	}
	pst->vec_mean[ic+par->n_comp*index_pix]=vec;
      }
      //We should solve the linear system here (FT*Nw*F)*x=FT*Nw*d
    }
  }
  //HERE
}
