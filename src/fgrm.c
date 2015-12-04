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
  pst->rand_spec=my_malloc(par->n_spec_vary*sizeof(flouble));
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
  free(pst->rand_spec);
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
  /*
    #define BS_M -1.1
    #define BD_M 1.55
    #define BS_S 0.03
    #define BD_S 0.03
    #define XICR 0.9
    flouble bs=x_spec[par->index_beta_s_t];
    flouble bd=x_spec[par->index_beta_d_t];
    double den=1./(BS_S*BS_S*BD_S*BD_S*(1-XICR*XICR));
    double chi2=((bs-BS_M)*(bs-BS_M)*BD_S*BD_S+
    (bd-BD_M)*(bd-BD_M)*BS_S*BS_S-
    2*(bs-BS_M)*(bd-BD_M)*XICR*BS_S*BD_S)*den;
  */

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
				 gsl_matrix *mat_step,flouble *x_spec_new)
{ //DAM: possible optimization: demote mat_step to flouble *
  int ipar;
  double chi2_old,chi2_new,ratio;
  memcpy(x_spec_new,x_spec_old,par->n_param_max*sizeof(flouble));
  for(ipar=0;ipar<par->n_spec_vary;ipar++)
    pst->rand_spec[ipar]=rand_gauss(pst->rng);
  for(ipar=0;ipar<par->n_spec_vary;ipar++) {
    int ipar2;
    for(ipar2=0;ipar2<=ipar;ipar2++)
      x_spec_new[ipar]+=gsl_matrix_get(mat_step,ipar,ipar2)*pst->rand_spec[ipar2];
  }

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

void clean_pixel(ParamFGRM *par,PixelState *pst,int ipix_big)
{
  int i_sample,ic1,ic2,ipix,n_updated,err;
  flouble ratio_accepted;
  int ip_spc=par->n_spec_vary*ipix_big;
  int id_cell=ipix_big*par->n_sub*par->n_pol;
  flouble *data=&(par->maps_data[id_cell*par->n_nu]);
  flouble *noise_w=&(par->maps_noise_weight[id_cell*par->n_nu]);
  flouble *amps_mean=&(par->map_components_mean[id_cell*par->n_comp]);
  flouble *amps_covar=&(par->map_components_covar[id_cell*par->n_comp*par->n_comp]);
  flouble *amps_dum=my_malloc(par->n_sub*par->n_pol*par->n_comp*sizeof(flouble));
  flouble *x_spec_old=my_malloc(par->n_param_max*sizeof(flouble));
  flouble *x_spec_new=my_malloc(par->n_param_max*sizeof(flouble));
  gsl_matrix *mat_step=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_save=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_spec=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  flouble *mean_spec=my_calloc(par->n_spec_vary,sizeof(flouble));
  flouble factor_rescale=2.4/sqrt((double)(par->n_spec_vary));

  memset(amps_mean,0,par->n_sub*par->n_pol*par->n_comp*sizeof(flouble));
  memset(amps_covar,0,par->n_sub*par->n_pol*par->n_comp*par->n_comp*sizeof(flouble));

  gsl_matrix_set_zero(mat_step);
  x_spec_old[par->index_beta_s_t]=par->beta_s_0;
  x_spec_old[par->index_beta_d_t]=par->beta_d_0;
  x_spec_old[par->index_temp_d_t]=par->temp_d_0;
  if(par->flag_beta_s_free)
    gsl_matrix_set(mat_step,par->index_beta_s_t,par->index_beta_s_t,par->sigma_beta_s);
  if(par->flag_beta_d_free)
    gsl_matrix_set(mat_step,par->index_beta_d_t,par->index_beta_d_t,par->sigma_beta_d);
  if(par->flag_temp_d_free)
    gsl_matrix_set(mat_step,par->index_temp_d_t,par->index_temp_d_t,par->sigma_temp_d);
  if(par->flag_include_polarization && par->flag_independent_polarization) {
    x_spec_old[par->index_beta_s_p]=par->beta_s_0;
    x_spec_old[par->index_beta_d_p]=par->beta_d_0;
    x_spec_old[par->index_temp_d_p]=par->temp_d_0;
    if(par->flag_beta_s_free)
      gsl_matrix_set(mat_step,par->index_beta_s_p,par->index_beta_s_p,par->sigma_beta_s);
    if(par->flag_beta_d_free)
      gsl_matrix_set(mat_step,par->index_beta_d_p,par->index_beta_d_p,par->sigma_beta_d);
    if(par->flag_temp_d_free)
      gsl_matrix_set(mat_step,par->index_temp_d_p,par->index_temp_d_p,par->sigma_temp_d);
  }

  dbg_printf("Burning\n");
  ratio_accepted=0;
  gsl_matrix_set_zero(cov_spec);
  gsl_matrix_set_zero(cov_save);
  n_updated=0;
  for(i_sample=0;i_sample<par->n_samples_burn;i_sample++) {
    int accepted;

    if(i_sample%par->n_update_covar==0)
      draw_amplitudes(par,data,noise_w,x_spec_old,pst,amps_dum); //A_{n+1}(b_n)
    accepted=draw_spectral_indices(par,data,noise_w,amps_dum,x_spec_old,pst,mat_step,x_spec_new); //b_{n+1}(A_{n+1})
    memcpy(x_spec_old,x_spec_new,par->n_param_max*sizeof(flouble));
#ifdef _DEBUG
    memcpy(&(par->dbg_extra[i_sample*par->n_spec_vary]),x_spec_old,par->n_spec_vary*sizeof(flouble));
#endif //_DEBUG
    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
      mean_spec[ic1]+=x_spec_old[ic1];
      for(ic2=0;ic2<par->n_spec_vary;ic2++) {
	flouble cov_plus=x_spec_old[ic1]*x_spec_old[ic2];
	gsl_matrix_set(cov_spec,ic1,ic2,gsl_matrix_get(cov_spec,ic1,ic2)+cov_plus);
      }
    }

    ratio_accepted+=accepted;
    if(i_sample%par->n_update_covar==par->n_update_covar-1) { //Update covariance
      ratio_accepted/=par->n_update_covar;
      dbg_printf("%d Acceptance ratio %.2lf\n",i_sample,ratio_accepted);

      //Compute mean in this batch
      dbg_printf("Current mean :");
      for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	mean_spec[ic1]/=par->n_update_covar;
	dbg_printf(" %lf",mean_spec[ic1]);
      }
      dbg_printf("\n");
      //Compute covariance in this batch
      dbg_printf("Current covariance :\n");
      for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	dbg_printf("   |");
      	for(ic2=0;ic2<par->n_spec_vary;ic2++) {
      	  flouble cov=gsl_matrix_get(cov_spec,ic1,ic2)/par->n_update_covar-
	    mean_spec[ic1]*mean_spec[ic2];
      	  gsl_matrix_set(cov_spec,ic1,ic2,cov);
	  dbg_printf(" %lE",gsl_matrix_get(cov_spec,ic1,ic2));
      	}
	dbg_printf("|\n");
      }
      //Add to mean covariance
      if(i_sample>par->n_samples_burn/2) {
	n_updated++;
	for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	  for(ic2=0;ic2<par->n_spec_vary;ic2++)
	    gsl_matrix_set(cov_save,ic1,ic2,gsl_matrix_get(cov_save,ic1,ic2)+gsl_matrix_get(cov_spec,ic1,ic2));
	}
      }
      //Save diagonal
      for(ic1=0;ic1<par->n_spec_vary;ic1++)
      	mean_spec[ic1]=sqrt(gsl_matrix_get(cov_spec,ic1,ic1));
      //Cholesky decomposition of covariance
      err=gsl_linalg_cholesky_decomp(cov_spec);
      gsl_matrix_set_zero(mat_step);
      if(err==GSL_EDOM) { //If covariance is not positive definite just save standard deviations
	dbg_printf("Covariance is not positive definite\n");
      	for(ic1=0;ic1<par->n_spec_vary;ic1++)
      	  gsl_matrix_set(mat_step,ic1,ic1,factor_rescale*mean_spec[ic1]);
      }
      else { //Store cholesky decomposition in stepping function
      	for(ic1=0;ic1<par->n_spec_vary;ic1++) {
      	  for(ic2=0;ic2<=ic1;ic2++)
      	    gsl_matrix_set(mat_step,ic1,ic2,factor_rescale*gsl_matrix_get(cov_spec,ic1,ic2));
      	}
      }

      gsl_matrix_set_zero(cov_spec);
      for(ic1=0;ic1<par->n_spec_vary;ic1++)
	mean_spec[ic1]=0;
      ratio_accepted=0;
    }
  }
  dbg_printf("Final covariance :\n");
  for(ic1=0;ic1<par->n_spec_vary;ic1++) { //Compute covariance
    dbg_printf("   |");
    for(ic2=0;ic2<par->n_spec_vary;ic2++) {
      flouble cov=gsl_matrix_get(cov_save,ic1,ic2)/n_updated;
      gsl_matrix_set(cov_save,ic1,ic2,cov);
      dbg_printf(" %lE",cov);
    }
    dbg_printf("|\n");
    mean_spec[ic1]=sqrt(gsl_matrix_get(cov_save,ic1,ic1)); //Save diagonal
  }
  err=gsl_linalg_cholesky_decomp(cov_save); //Cholesky decomposition of covariance
  gsl_matrix_set_zero(mat_step);
  if(err==GSL_EDOM) { //If covariance is not positive definite just save standard deviations
    dbg_printf("Final covariance is not positive definite!!!!\n");
    for(ic1=0;ic1<par->n_spec_vary;ic1++)
      gsl_matrix_set(mat_step,ic1,ic1,factor_rescale*mean_spec[ic1]);
  }
  else { //Store cholesky decomposition in stepping function
    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
      for(ic2=0;ic2<=ic1;ic2++)
	gsl_matrix_set(mat_step,ic1,ic2,factor_rescale*gsl_matrix_get(cov_save,ic1,ic2));
    }
  }

  dbg_printf("Starting actual sampling\n");
#ifdef _DEBUG
  ratio_accepted=0;
#endif //_DEBUG
  for(i_sample=0;i_sample<par->n_samples;i_sample++) {
    int accepted;
    if(i_sample%par->n_spec_resample==0)
      draw_amplitudes(par,data,noise_w,x_spec_old,pst,amps_dum); //A_{n+1}(b_n)
    accepted=draw_spectral_indices(par,data,noise_w,amps_dum,x_spec_old,pst,mat_step,x_spec_new); //b_{n+1}(A_{n+1})
    memcpy(x_spec_old,x_spec_new,par->n_param_max*sizeof(flouble));

#ifdef _DEBUG
    memcpy(&(par->dbg_extra[(i_sample+par->n_samples_burn)*par->n_spec_vary]),x_spec_old,par->n_spec_vary*sizeof(flouble));
    ratio_accepted+=accepted;
    if(i_sample && (i_sample%par->n_update_covar==par->n_update_covar-1)) {
      ratio_accepted/=par->n_update_covar;
      dbg_printf("%d Acceptance ratio %.2lf\n",i_sample,ratio_accepted);
      ratio_accepted=0;
    }
#endif //_DEBUG

    //Add amplitudes
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
  free(mean_spec);
  gsl_matrix_free(mat_step);
  gsl_matrix_free(cov_save);
  gsl_matrix_free(cov_spec);
}
