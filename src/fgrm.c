#include "common.h"

PixelState *pixel_state_new(ParamFGRM *par,unsigned long seed)
{
  int ip;
  PixelState *pst=my_malloc(sizeof(PixelState));
  pst->f_matrix=my_malloc(par->n_pol*par->n_nu*par->n_comp*sizeof(flouble));
  pst->rng=init_rng(seed);
  pst->cov_inv=my_malloc(par->n_pix*par->n_pol*sizeof(gsl_matrix *));
  pst->vec_mean=my_malloc(par->n_pix*par->n_pol*sizeof(gsl_vector *));

  for(ip=0;ip<par->n_pix*par->n_pol;ip++) {
    pst->cov_inv[ip]=gsl_matrix_alloc(par->n_comp,par->n_comp);
    pst->vec_mean[ip]=gsl_vector_alloc(par->n_comp);
  }

  return pst;
}

void pixel_state_free(PixelState *pst,ParamFGRM *par)
{
  int ip;
  free(pst->f_matrix);
  end_rng(pst->rng);
  for(ip=0;ip<par->n_pix*par->n_pol;ip++) {
    gsl_matrix_free(pst->cov_inv[ip]);
    gsl_vector_free(pst->vec_mean[ip]);
  }
  free(pst->cov_inv);
  free(pst->vec_mean);
  free(pst);
}

static flouble freq_evolve(int spec_type,double nu_0,double beta,double temp,double nu)
{
  flouble x_to,x_from,ex;
  switch(spec_type)  {
  case 0 :
    x_to=0.017611907*nu;
    ex=exp(x_to);
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
  return -1;
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

static double compute_chi2(ParamFGRM *par,flouble *data,flouble *noise_w,
			    flouble *amps,flouble *x_spec,PixelState *pst)
{
  int ipix;
  double chi2=0;

  compute_f_matrix(par,x_spec,pst->f_matrix);

  for(ipix=0;ipix<par->n_sub;ipix++) {
    int ipol;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int inu;
      int index_pix=ipol+par->n_pol*ipix;
      for(inu=0;inu<par->n_nu;inu++) {
	int icomp;
	double res=(double)(data[inu+par->n_nu*index_pix]);
	for(icomp=0;icomp<par->n_comp;icomp++) {
	  res-=pst->f_matrix[icomp+par->n_comp*(inu+par->n_nu*ipol)]*
	    amps[icomp+par->n_comp*index_pix];
	}
	chi2+=res*res*noise_w[inu+par->n_nu*index_pix];
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
	int ic2,inu;
	flouble vec=0;
	for(inu=0;inu<par->n_nu;inu++) {
	  vec+=
	    pst->f_matrix[ic1+par->n_comp*(inu+par->n_nu*ipol)]*
	    data[inu+par->n_nu*index_pix]*noise_w[inu+par->n_nu*index_pix]; //v = F^T N^-1 d
	}
	gsl_vector_set(pst->vec_mean[index_pix],ic1,vec);
	for(ic2=0;ic2<par->n_comp;ic2++) {
	  flouble icov=0;
	  for(inu=0;inu<par->n_nu;inu++) {
	    icov+=
	      pst->f_matrix[ic1+par->n_comp*(inu+par->n_nu*ipol)]*
	      pst->f_matrix[ic2+par->n_comp*(inu+par->n_nu*ipol)]*
	      noise_w[inu+par->n_nu*index_pix]; //C^-1 = F^T N^-1 F
	  }
	  gsl_matrix_set(pst->cov_inv[index_pix],ic1,ic2,icov);
	}
      }
      gsl_linalg_cholesky_decomp(pst->cov_inv[index_pix]);
      gsl_linalg_cholesky_svx(pst->cov_inv[index_pix],pst->vec_mean[index_pix]); //v = (F^T N^-1 F)^-1 F^T N^-1 d
    }
  }
}

//static void solve_lower_triangular(int n,gsl_matrix *mat,flouble *v)
//{
//  int i;
//  for(i=0;i<n;i++) {
//    int j;
//    flouble res=v[i];
//    for(j=0;j<i;j++)
//      res-=v[j]*gsl_matrix_get(mat,i,j);
//    v[i]=res/gsl_matrix_get(mat,i,i);
//  }
//}

static void solve_upper_triangular(int n,gsl_matrix *mat,flouble *v)
{
  int i;
  for(i=n-1;i>=0;i--) {
    int j;
    flouble res=v[i];
    for(j=n-1;j>i;j--)
      res-=gsl_matrix_get(mat,i,j)*v[j];
    v[i]=res/gsl_matrix_get(mat,i,i);
  }
}

static void draw_amplitudes(ParamFGRM *par,flouble *data,flouble *noise_w,
			    flouble *x_spec,PixelState *pst,flouble *amps)
{
  int ipix;

  analyze_linear_chi2(par,data,noise_w,x_spec,pst);
  
  for(ipix=0;ipix<par->n_sub;ipix++) {
    int ipol;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int icomp;
      int index_pix=ipol+par->n_pol*ipix;
      flouble *amp_here=&(amps[par->n_comp*index_pix]);
      for(icomp=0;icomp<par->n_comp;icomp++)
	amp_here[icomp]=rand_gauss(pst->rng);
      solve_upper_triangular(par->n_comp,pst->cov_inv[index_pix],amp_here); //t_rand = L^T^-1 u
      for(icomp=0;icomp<par->n_comp;icomp++)
	amp_here[icomp]+=gsl_vector_get(pst->vec_mean[index_pix],icomp); //t = t_mean + t_rand
    }
  }
}

static int draw_spectral_indices(ParamFGRM *par,flouble *data,flouble *noise_w,
				 flouble *amps,flouble *x_spec_old,PixelState *pst,
				 flouble *sigma_step,flouble *x_spec_new)
{
  int ipar;
  double chi2_old,chi2_new,ratio;
  memcpy(x_spec_new,x_spec_old,par->n_param_max*sizeof(flouble));
  for(ipar=0;ipar<par->n_spec_vary;ipar++)
    x_spec_new[ipar]+=sigma_step[ipar]*rand_gauss(pst->rng);

  chi2_new=compute_chi2(par,data,noise_w,amps,x_spec_new,pst);
  chi2_old=compute_chi2(par,data,noise_w,amps,x_spec_old,pst);
      
  ratio=exp(-0.5*(chi2_new-chi2_old));
  
  if(ratio<1) {
    if(rand_real01(pst->rng)>ratio) {
      memcpy(x_spec_new,x_spec_old,par->n_spec_vary*sizeof(flouble));
      return 0;
    }
  }

  //  for(ipar=0;ipar<par->n_spec_vary;ipar++)
  //    printf(" (%lf %lf)",x_spec_old[ipar],x_spec_new[ipar]);
  //  printf(" | %lf %lf %lf\n",chi2_old,chi2_new,ratio);

  return 1;
}

#define N_CHECK 1000
#define RATE_TARGET 0.1

flouble *clean_pixel(ParamFGRM *par,PixelState *pst,int ipix_big,int n_samples)
{
  int i_sample,ic1,ic2,ipix;
  int ip_spc=par->n_spec_vary*ipix_big;
  flouble *data=&(par->maps_data[ipix_big*par->n_sub*par->n_pol*par->n_nu]);
  flouble *noise_w=&(par->maps_noise_weight[ipix_big*par->n_sub*par->n_pol*par->n_nu]);
  flouble *amps_mean=&(par->map_components_mean[ipix_big*par->n_sub*par->n_pol*par->n_comp]);
  flouble *amps_covar=&(par->map_components_covar[ipix_big*par->n_sub*par->n_pol*par->n_comp*par->n_comp]);
  flouble *amps_dum=my_malloc(par->n_sub*par->n_pol*par->n_comp*sizeof(flouble));
#ifdef _DEBUG
  flouble *x_spec_samples=my_malloc(par->n_spec_vary*n_samples*sizeof(flouble));
#endif //_DEBUG
  flouble *x_spec_old=my_malloc(par->n_param_max);
  flouble *x_spec_new=my_malloc(par->n_param_max);
  flouble *sigma_step=my_malloc(par->n_param_max);

  memset(amps_mean,0,par->n_sub*par->n_pol*par->n_comp*sizeof(flouble));
  memset(amps_covar,0,par->n_sub*par->n_pol*par->n_comp*par->n_comp*sizeof(flouble));

  x_spec_old[par->index_beta_s_t]=par->beta_s_0;
  x_spec_old[par->index_beta_d_t]=par->beta_d_0;
  x_spec_old[par->index_temp_d_t]=par->temp_d_0;
  sigma_step[par->index_beta_s_t]=par->sigma_beta_s;
  sigma_step[par->index_beta_d_t]=par->sigma_beta_d;
  sigma_step[par->index_temp_d_t]=par->sigma_temp_d;
  if(par->flag_include_polarization && par->flag_independent_polarization) {
    x_spec_old[par->index_beta_s_t]=par->beta_s_0;
    x_spec_old[par->index_beta_d_t]=par->beta_d_0;
    x_spec_old[par->index_temp_d_t]=par->temp_d_0;
    sigma_step[par->index_beta_s_t]=par->sigma_beta_s;
    sigma_step[par->index_beta_d_t]=par->sigma_beta_d;
    sigma_step[par->index_temp_d_t]=par->sigma_temp_d;
  }
    
  flouble ratio_accepted=0;
  for(i_sample=0;i_sample<n_samples;i_sample++) {
    int accepted;
    draw_amplitudes(par,data,noise_w,x_spec_old,pst,amps_dum); //A_{n+1}(b_n)
    accepted=draw_spectral_indices(par,data,noise_w,amps_dum,x_spec_old,pst,sigma_step,x_spec_new); //b_{n+1}(A_{n+1})

    ratio_accepted+=accepted;
    if(i_sample%N_CHECK==0) {
      flouble factor_rescale;
      ratio_accepted/=N_CHECK;
      printf("%d Acceptance ratio %lE\n",i_sample,ratio_accepted);
      if(ratio_accepted>RATE_TARGET*0.1)
	factor_rescale=ratio_accepted/RATE_TARGET;
      else
	factor_rescale=0.5;
      for(ic1=0;ic1<par->n_spec_vary;ic1++)
	sigma_step[ic1]*=factor_rescale;
      ratio_accepted=0;
    }

    //#ifdef _DEBUG
    //    printf("%d %d\n",i_sample,accepted);
    //#endif //_DEBUG
    memcpy(x_spec_old,x_spec_new,par->n_param_max*sizeof(flouble));

    //Add amplitudes
#ifdef _DEBUG
    memcpy(&(x_spec_samples[i_sample*par->n_spec_vary]),x_spec_old,par->n_spec_vary*sizeof(flouble));
#endif //_DEBUG
    for(ipix=0;ipix<par->n_sub*par->n_pol;ipix++) {
      int index_here=ipix*par->n_comp;
      for(ic1=0;ic1<par->n_comp;ic1++) {
	amps_mean[ic1+index_here]+=amps_dum[ic1+index_here];
	for(ic2=0;ic2<par->n_comp;ic2++)
	  amps_covar[ic2+par->n_comp*(ic1+index_here)]+=amps_dum[ic1+index_here]*amps_dum[ic2+index_here];
      }
    }

    //Add indices
    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
      par->map_indices_mean[ic1+ip_spc]+=x_spec_old[ic1];
      for(ic2=0;ic2<par->n_spec_vary;ic2++)
	par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]+=x_spec_old[ic1]*x_spec_old[ic2];
    }
  }

  //Compute averages
  for(ipix=0;ipix<par->n_sub*par->n_pol;ipix++) {
    int index_here=ipix*par->n_comp;
    for(ic1=0;ic1<par->n_comp;ic1++)
      amps_mean[ic1+index_here]/=par->n_samples;
  }
  for(ipix=0;ipix<par->n_sub*par->n_pol;ipix++) {
    int index_here=ipix*par->n_comp;
    for(ic1=0;ic1<par->n_comp;ic1++) {
      for(ic2=0;ic2<par->n_comp;ic2++) {
	amps_covar[ic2+par->n_comp*(ic1+index_here)]/=par->n_samples;
	amps_covar[ic2+par->n_comp*(ic1+index_here)]-=amps_mean[ic1+index_here]*amps_mean[ic2+index_here];
      }
    }
  }
  for(ic1=0;ic1<par->n_spec_vary;ic1++)
    par->map_indices_mean[ic1+ip_spc]/=par->n_samples;
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    for(ic2=0;ic2<par->n_spec_vary;ic2++) {
      par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]/=par->n_samples;
      par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]-=
	par->map_indices_mean[ic1+ip_spc]*par->map_indices_mean[ic1+ip_spc];
    }
  }

  free(amps_dum);
  free(x_spec_old);
  free(x_spec_new);
  free(sigma_step);

#ifdef _DEBUG
  return x_spec_samples;
#else //_DEBUG
  return NULL;
#endif //_DEBUG
}
