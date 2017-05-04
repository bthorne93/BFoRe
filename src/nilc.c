#include "common.h"

static void output_ilc_band(ParamBFoRe *par,he_needlet_params *ntp)
{
  int ii,lmx;
  char fname[256];
  FILE *fo;

  lmx=3*ntp->nside0-1;
  sprintf(fname,"%s_ILC_bands.txt",par->output_prefix);
  fo=my_fopen(fname,"w");
  for(ii=0;ii<=lmx;ii++) {
    int jj;
    double bsum=0;
    fprintf(fo,"%d ",ii);
    for(jj=0;jj<ntp->nj;jj++) {
      double b=ntp->b_arr[jj][ii];
      bsum+=b*b;
      fprintf(fo,"%lE ",b);
    }
    fprintf(fo,"%lE\n",bsum);
  }
  fclose(fo);

  /*
#ifdef _DEBUG
  flouble *mp,*mpb,***nt;//,*mp_rec;
  flouble *cl=my_malloc((lmx+1)*sizeof(flouble));
  FILE *fi=my_fopen("../../cmb/cls/planck1_r0p00_lensedCls.dat","r");
  for(ii=0;ii<lmx-1;ii++) {
    double l,dum,c;
    int stat=fscanf(fi,"%lf %lf %lf %lf %lf",&l,&c,&dum,&dum,&dum);
    if(stat!=5)
      report_error(1,"WTH\n");
    cl[ii]=c*2*M_PI/(l*(l+1));
  }
  fclose(fi);
  mpb=my_calloc(12*ntp->nside0*ntp->nside0,sizeof(flouble));
  mp=he_synfast(cl,ntp->nside0,lmx,1234);
  nt=he_alloc_needlet(ntp,0);
  he_map2needlet(ntp,&mp,nt,0,0,0,0);
  he_needlet2map(ntp,&mpb,nt,0,0,0,0);

  sprintf(fname,"!%s_ILC_mapsim.fits",par->output_prefix);
  he_write_healpix_map(&mp,1,ntp->nside0,fname);
  sprintf(fname,"!%s_ILC_mapsim_back.fits",par->output_prefix);
  he_write_healpix_map(&mpb,1,ntp->nside0,fname);
  for(ii=0;ii<ntp->nj;ii++) {
    sprintf(fname,"!%s_ILC_nt%03d.fits",par->output_prefix,ii+1);
    he_write_healpix_map(&(nt[ii][0]),1,ntp->nside_arr[ii],fname);
  }
  he_free_needlet(ntp,0,nt);
  free(mp);
  free(mpb);
#endif //_DEBUG
*/
}

static void do_nilc_discs(ParamBFoRe *par,he_needlet_params *ntp,flouble ****nt_in,flouble ***nt_out)
{
  int ii;
#ifdef _DEBUG
  char fname[256];
#endif //_DEBUG
  gsl_vector *vec_a=gsl_vector_alloc(par->n_nu);

  //Compute CMB frequency dependence
  for(ii=0;ii<par->n_nu;ii++)
    gsl_vector_set(vec_a,ii,freq_evolve(0,-1.,-1.,-1.,par->freqs[ii]));

#ifdef _DEBUG
  for(ii=0;ii<par->n_nu;ii++) {
    int jj;
    for(jj=0;jj<ntp->nj;jj++) {
      long nside_j=ntp->nside_arr[jj];
      sprintf(fname,"!%s_nt_nu%03d_j%03d.fits",par->output_prefix,ii+1,jj+1);
      he_write_healpix_map(nt_in[ii][jj],par->n_pol,nside_j,fname);
    }
  }
#endif //_DEBUG

  //Do ILC for each needlet scale
  for(ii=0;ii<ntp->nj;ii++) {
    int ipol;
    int nside_j=ntp->nside_arr[ii];
    int npix_j=12*nside_j*nside_j;
    double dlmx=pow(ntp->b,ii+1+ntp->jmax_min);
    double theta_cutoff=M_PI/dlmx;
    double theta_disc_1=MIN((0.5*M_PI),(par->nilc_ndim*theta_cutoff));
    double theta_disc_2=acos(1-(double)(par->nilc_ndim*par->nilc_ndim)/(6.*nside_j*nside_j));
    double theta_disc=MAX(theta_disc_1,theta_disc_2);
    int npix_disc=(int)(1.2*12*nside_j*nside_j*0.5*(1-cos(theta_disc)));
    for(ipol=0;ipol<par->n_pol;ipol++) {
      printf("%d %d %lf %lf %d\n",ii,ipol,theta_cutoff*180/M_PI,theta_disc*180/M_PI,npix_disc);
#pragma omp parallel default(none)					\
  shared(npix_disc,ntp,par,nt_in,nt_out,nside_j,npix_j,theta_disc,ii,ipol,vec_a)
      {
	int ip;
	double inv_a_Ni_a;
	int *listpix=my_malloc(2*npix_disc*sizeof(int));
	gsl_vector *vec_d=gsl_vector_alloc(par->n_nu);
	gsl_vector *vec_Ni_a=gsl_vector_alloc(par->n_nu);
	gsl_matrix *ncov=gsl_matrix_alloc(par->n_nu,par->n_nu);

#pragma omp for
	for(ip=0;ip<npix_j;ip++) {
	  int ipd,in1,in2,nplist_here;
	  double cth,phi,map_value;
	  nplist_here=npix_disc;
	  
	  pix2ang_ring(nside_j,ip,&cth,&phi);
	  cth=cos(cth);
	  he_query_disc(nside_j,cth,phi,theta_disc,listpix,&nplist_here,0);
	  
	  //Compute covariance matrix
	  for(in1=0;in1<par->n_nu;in1++) {
	    for(in2=0;in2<=in1;in2++) {
	      double cov_element=0;
	      for(ipd=0;ipd<nplist_here;ipd++) {
		int ipix=listpix[ipd];
		cov_element+=nt_in[in1][ii][ipol][ipix]*nt_in[in2][ii][ipol][ipix];
	      }
	      cov_element/=nplist_here;
	      gsl_matrix_set(ncov,in1,in2,cov_element);
	      if(in1!=in2)
		gsl_matrix_set(ncov,in2,in1,cov_element);
	    }
	  }
	  
	  //Invert covariance
	  gsl_linalg_cholesky_decomp(ncov);
	  gsl_linalg_cholesky_solve(ncov,vec_a,vec_Ni_a); //a^T N^-1
	  gsl_blas_ddot(vec_Ni_a,vec_a,&inv_a_Ni_a); //a^T N^-1 a
	  inv_a_Ni_a=1./inv_a_Ni_a;
	  
	  //Compute ILC
	  for(in1=0;in1<par->n_nu;in1++)
	    gsl_vector_set(vec_d,in1,nt_in[in1][ii][ipol][ip]);
	  
	  //(a^T N^-1 d) / (a^T N^-1 a)
	  gsl_blas_ddot(vec_Ni_a,vec_d,&map_value);
	  nt_out[ii][ipol][ip]=map_value*inv_a_Ni_a;
	}
	free(listpix);
	gsl_matrix_free(ncov);
	gsl_vector_free(vec_d);
	gsl_vector_free(vec_Ni_a);
      }
    }
  }
  gsl_vector_free(vec_a);
  
#ifdef _DEBUG
  for(ii=0;ii<ntp->nj;ii++) {
    long nside_j=ntp->nside_arr[ii];
    sprintf(fname,"!%s_nt_out_j%03d.fits",par->output_prefix,ii+1);
    he_write_healpix_map(nt_out[ii],par->n_pol,nside_j,fname);
  }
#endif //_DEBUG
}

static void do_nilc_nested(ParamBFoRe *par,he_needlet_params *ntp,flouble ****nt_in,flouble ***nt_out)
{
  int ii;
#ifdef _DEBUG
  char fname[256];
#endif //_DEBUG
  double inv_a_Ni_a;
  gsl_vector *vec_d=gsl_vector_alloc(par->n_nu);
  gsl_vector *vec_Ni_a=gsl_vector_alloc(par->n_nu);
  gsl_matrix *ncov=gsl_matrix_alloc(par->n_nu,par->n_nu);
  gsl_vector *vec_a=gsl_vector_alloc(par->n_nu);

  //Compute CMB frequency dependence
  for(ii=0;ii<par->n_nu;ii++)
    gsl_vector_set(vec_a,ii,freq_evolve(0,-1.,-1.,-1.,par->freqs[ii]));

  //Ring2nest
  for(ii=0;ii<par->n_nu;ii++) {
    int jj;
    for(jj=0;jj<ntp->nj;jj++) {
      int ipol;
      long nside_j=ntp->nside_arr[jj];
#ifdef _DEBUG
      sprintf(fname,"!%s_nt_nu%03d_j%03d.fits",par->output_prefix,ii+1,jj+1);
      he_write_healpix_map(nt_in[ii][jj],par->n_pol,nside_j,fname);
#endif //_DEBUG
      for(ipol=0;ipol<par->n_pol;ipol++)
	he_ring2nest_inplace(nt_in[ii][jj][ipol],nside_j);
    }
  }
    
  //Do ILC for each needlet scale
  for(ii=0;ii<ntp->nj;ii++) {
    int ipol;
    int nside_j=ntp->nside_arr[ii];
    int nside_sub=nside_j/par->nilc_ndim;
    int npix_sub=12*nside_sub*nside_sub;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int ip_sub;
      for(ip_sub=0;ip_sub<npix_sub;ip_sub++) {
	int ip,in1,in2;
	
	//Compute covariance matrix
	for(in1=0;in1<par->n_nu;in1++) {
	  for(in2=0;in2<=in1;in2++) {
	    double cov_element=0;
	    for(ip=0;ip<par->nilc_ndim*par->nilc_ndim;ip++) {
	      int ipix=par->nilc_ndim*par->nilc_ndim*ip_sub+ip;
	      cov_element+=nt_in[in1][ii][ipol][ipix]*nt_in[in2][ii][ipol][ipix];
	    }
	    cov_element/=(par->nilc_ndim*par->nilc_ndim);
	    gsl_matrix_set(ncov,in1,in2,cov_element);
	    if(in1!=in2)
	      gsl_matrix_set(ncov,in2,in1,cov_element);
	  }
	}
	  
	//Invert covariance
	gsl_linalg_cholesky_decomp(ncov);
	gsl_linalg_cholesky_solve(ncov,vec_a,vec_Ni_a); //a^T N^-1
	gsl_blas_ddot(vec_Ni_a,vec_a,&inv_a_Ni_a); //a^T N^-1 a
	inv_a_Ni_a=1./inv_a_Ni_a;
	  
	//Compute ILC
	for(ip=0;ip<par->nilc_ndim*par->nilc_ndim;ip++) {
	  double map_value;
	  int ipix=par->nilc_ndim*par->nilc_ndim*ip_sub+ip;
	  for(in1=0;in1<par->n_nu;in1++)
	    gsl_vector_set(vec_d,in1,nt_in[in1][ii][ipol][ipix]);
	    
	  //(a^T N^-1 d) / (a^T N^-1 a)
	  gsl_blas_ddot(vec_Ni_a,vec_d,&map_value);
	  nt_out[ii][ipol][ipix]=map_value*inv_a_Ni_a;
	}
      }
    }
  }
  gsl_vector_free(vec_d);
  gsl_vector_free(vec_Ni_a);
  gsl_matrix_free(ncov);
  gsl_vector_free(vec_a);
  
  //Nest2ring
  for(ii=0;ii<ntp->nj;ii++) {
    int ipol;
    long nside_j=ntp->nside_arr[ii];
    for(ipol=0;ipol<par->n_pol;ipol++)
      he_nest2ring_inplace(nt_out[ii][ipol],nside_j);
#ifdef _DEBUG
    sprintf(fname,"!%s_nt_out_j%03d.fits",par->output_prefix,ii+1);
    he_write_healpix_map(nt_out[ii],par->n_pol,nside_j,fname);
#endif //_DEBUG
  }
}

void do_nilc(ParamBFoRe *par)
{
  int ii;
  flouble ****nt_in,***nt_out;
  he_needlet_params *ntp;

  printf("NILC\n");
  ntp=he_nt_init(par->b_nt,par->nside,par->nilc_niter);
  output_ilc_band(par,ntp);

  //Compute needlet coefficients
  printf("  Computing needlet coefficients\n");
  nt_out=he_alloc_needlet(ntp,par->flag_include_polarization);
  nt_in=my_malloc(par->n_nu*sizeof(flouble ***));
  for(ii=0;ii<par->n_nu;ii++) {
    nt_in[ii]=he_alloc_needlet(ntp,par->flag_include_polarization);
    he_map2needlet(ntp,par->map_nilc_in[ii],nt_in[ii],0,par->flag_include_polarization,0,1);
  }

  //NILC
  printf("Computing ILC\n");
  int do_nested=0;
  if(do_nested)
    do_nilc_nested(par,ntp,nt_in,nt_out);
  else 
    do_nilc_discs(par ,ntp,nt_in,nt_out);

  //Synthesize needlets
  printf("  Synthesizing needlets\n");
  he_needlet2map(ntp,par->map_nilc_out,nt_out,0,par->flag_include_polarization,1,1);

  for(ii=0;ii<par->n_nu;ii++)
    he_free_needlet(ntp,par->flag_include_polarization,nt_in[ii]);
  he_free_needlet(ntp,par->flag_include_polarization,nt_out);
  free(nt_in);
  he_nt_end(ntp);
}
