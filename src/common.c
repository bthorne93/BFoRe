#include "common.h"

int NNodes=1;
int NodeThis=0;
int IThread0=0;

size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  if(fwrite(ptr,size,nmemb,stream)!=nmemb)
    report_error(1,"Error fwriting\n");

  return nmemb;
}

int my_linecount(FILE *f)
{
  int i0=0;
  char ch[1000];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

void dbg_printf(char *fmt,...)
{
#ifdef _DEBUG
  va_list args;
  char msg[256];

  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);
  printf("%s",msg);
#endif //_DEBUG
}

void report_error(int level,char *fmt,...)
{
  va_list args;
  char msg[256];

  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);
  
  if(level) {
    fprintf(stderr," Fatal error: %s",msg);
    exit(level);
  }
  else
    fprintf(stderr," Warning: %s",msg);
}

void *my_malloc(size_t size)
{
  void *outptr=malloc(size);
  if(outptr==NULL) report_error(1,"Out of memory\n");

  return outptr;
}

void *my_calloc(size_t nmemb,size_t size)
{
  void *outptr=calloc(nmemb,size);
  if(outptr==NULL)
    report_error(1,"Out of memory\n");

  return outptr;
}

FILE *my_fopen(const char *path,const char *mode)
{
  FILE *fout=fopen(path,mode);
  if(fout==NULL)
    report_error(1,"Couldn't open file %s\n",path);

  return fout;
}

static ParamBFoRe *param_bfore_new(void)
{
  ParamBFoRe *par=(ParamBFoRe *)my_malloc(sizeof(ParamBFoRe));

  par->nside=256;
  par->nside_spec=32;
  par->n_side_sub=8;
  par->n_sub=64;
  par->n_pix=786432;
  
  sprintf(par->input_data_prefix,"default");
  par->maps_data=NULL;

  sprintf(par->input_noise_prefix,"default");
  par->maps_noise_weight=NULL;

  sprintf(par->output_prefix,"default");
  par->map_components_mean=NULL;
  par->map_components_covar=NULL;
  par->map_indices_mean=NULL;
  par->map_indices_covar=NULL;
  par->map_chi2=NULL;

  sprintf(par->fname_nulist,"default");
  par->n_nu=-1;
  par->freqs=NULL;
  
  par->flag_include_polarization=0;
  par->n_pol=1;

  par->flag_include_cmb=0;
  par->flag_include_synchrotron=0;
  par->flag_include_dust=0;
  par->n_comp=0;
  par->index_cmb=-1;
  par->index_synchrotron=-1;
  par->index_dust=-1;

  par->flag_independent_polarization=0;
  par->flag_beta_s_free=0;
  par->flag_beta_d_free=0;
  par->flag_temp_d_free=0;
  par->n_param_max=0;
  par->n_spec_vary=0;
  par->n_dof_pix=0;
  par->index_beta_s_t=-1;
  par->index_beta_s_p=-1;
  par->index_beta_d_t=-1;
  par->index_beta_d_p=-1;
  par->index_temp_d_t=-1;
  par->index_temp_d_p=-1;
  
  par->beta_s_0=-1.;
  par->beta_d_0=1.54;
  par->temp_d_0=20.;
  par->beta_s_step=0.1;
  par->beta_d_step=0.1;
  par->temp_d_step=0.1;
  par->beta_s_prior=-1;
  par->beta_d_prior=-1;
  par->temp_d_prior=-1;
  par->nu0_s=23.;
  par->nu0_d=353.;

  par->seed=1234;
  par->n_samples=100000;
  par->frac_samples_burn=0.2;
  par->n_update_covar=1000;
  par->n_samples_burn=20000;
  par->n_spec_resample=1;

  par->dbg_ipix=0;
  par->dbg_extra=NULL;
  
  return par;
}

static void param_bfore_print(ParamBFoRe *par)
{
  printf("Read parameters:\n");
  printf(" - Nside = %d\n",par->nside);
  printf(" - Nside_spec = %d\n",par->nside_spec);
  printf(" - Input maps: %s\n",par->input_data_prefix);
  printf(" - Input noise: %s\n",par->input_noise_prefix);
  printf(" - Output prefix: %s\n",par->output_prefix);
  printf(" - Frequency list: %s\n",par->fname_nulist);
  if(par->flag_include_polarization)
    printf(" - Will include T, Q and U\n");
  else 
    printf(" - Will only include T\n");
  if(par->flag_independent_polarization)
    printf(" - Independent spectral parameters for polarization\n");
  else
    printf(" - Same spectral parameters for T, Q and U\n");
  printf(" - %d components include:",par->n_comp);
  if(par->flag_include_cmb) printf(" CMB(%d)",par->index_cmb);
  if(par->flag_include_synchrotron) printf(" Synchrotron(%d)",par->index_synchrotron);
  if(par->flag_include_dust) printf(" Dust(%d)",par->index_dust);
  printf("\n");
  printf(" - Free spectral parameters:");
  if(par->flag_include_synchrotron && par->flag_beta_s_free) {
    if(par->flag_independent_polarization)
      printf(" beta_s(%d,%d)",par->index_beta_s_t,par->index_beta_s_p);
    else
      printf(" beta_s(%d)",par->index_beta_s_t);
  }      
  if(par->flag_include_dust && par->flag_beta_d_free) {
    if(par->flag_independent_polarization)
      printf(" beta_d(%d,%d)",par->index_beta_d_t,par->index_beta_d_p);
    else
      printf(" beta_d(%d)",par->index_beta_d_t);
  }      
  if(par->flag_include_dust && par->flag_temp_d_free) {
    if(par->flag_independent_polarization)
      printf(" temp_d(%d,%d)",par->index_temp_d_t,par->index_temp_d_p);
    else
      printf(" temp_d(%d)",par->index_temp_d_t);
  }
  printf(". %d out of %d\n",par->n_spec_vary,par->n_param_max);
  printf(" - %d DOF per pixel\n",par->n_dof_pix);
  printf(" - beta_s_0 = %.3lf, ",par->beta_s_0);
  if(par->beta_s_prior>0)
    printf(" S(beta_s) = %.3lf",par->beta_s_prior);
  printf(" D(beta_s) = %.3lf\n",par->beta_s_step);
  printf(" - beta_d_0 = %.3lf, ",par->beta_d_0);
  if(par->beta_d_prior>0)
    printf(" S(beta_d) = %.3lf",par->beta_d_prior);
  printf(" D(beta_d) = %.3lf\n",par->beta_d_step);
  printf(" - temp_d_0 = %.3lf, ",par->temp_d_0);
  if(par->temp_d_prior>0)
    printf(" S(temp_d) = %.3lf",par->temp_d_prior);
  printf(" D(temp_d) = %.3lf\n",par->temp_d_step);
  printf(" - Seed = %lu\n",par->seed);
  printf(" - Will take %d samples\n",par->n_samples);
  printf("   after %d burning steps\n",par->n_samples_burn);
  printf(" - Proposal covariance will be updated every %d burning steps\n",
	 par->n_update_covar);
  printf(" - Amplitudes will be resampled every %d steps\n",par->n_spec_resample);
#ifdef _DEBUG
  printf(" - Will print debug information for pixel %d\n",par->dbg_ipix);
#endif //_DEBUG
}

void param_bfore_free(ParamBFoRe *par)
{
  free(par->maps_data);
  free(par->maps_noise_weight);
  free(par->map_components_mean);
  free(par->map_components_covar);
  free(par->map_indices_mean);
  free(par->map_indices_covar);
  free(par->map_chi2);
  free(par->freqs);
  free(par->prior_mean);
  free(par->prior_isigma);
  free(par);
}

ParamBFoRe *read_params(char *fname)
{
  FILE *fi;
  char fname_in[256];
  int n_lin,ii;
  long nside_dum;
  ParamBFoRe *par=param_bfore_new();
  flouble *map_dum;

  //Read parameters from file
  if(NodeThis==0)
    printf("*** Reading run parameters \n");
  fi=my_fopen(fname,"r");
  n_lin=my_linecount(fi); rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      report_error(1,"Error reading file %s, line %d\n",fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      report_error(1,"Error reading file %s, line %d\n",fname,ii+1);

    if(!strcmp(s1,"input_data_prefix="))
      sprintf(par->input_data_prefix,"%s",s2);
    else if(!strcmp(s1,"input_noise_prefix="))
      sprintf(par->input_noise_prefix,"%s",s2);
    else if(!strcmp(s1,"output_prefix="))
      sprintf(par->output_prefix,"%s",s2);
    else if(!strcmp(s1,"fname_nulist="))
      sprintf(par->fname_nulist,"%s",s2);

    else if(!strcmp(s1,"include_polarization="))
      par->flag_include_polarization=atoi(s2);
    else if(!strcmp(s1,"include_cmb="))
      par->flag_include_cmb=atoi(s2);
    else if(!strcmp(s1,"include_synchrotron="))
      par->flag_include_synchrotron=atoi(s2);
    else if(!strcmp(s1,"include_dust="))
      par->flag_include_dust=atoi(s2);
    else if(!strcmp(s1,"independent_polarization="))
      par->flag_independent_polarization=atoi(s2);
    else if(!strcmp(s1,"beta_s_free="))
      par->flag_beta_s_free=atoi(s2);
    else if(!strcmp(s1,"beta_d_free="))
      par->flag_beta_d_free=atoi(s2);
    else if(!strcmp(s1,"temp_d_free="))
      par->flag_temp_d_free=atoi(s2);
    else if(!strcmp(s1,"nu0_synchrotron="))
      par->nu0_s=atof(s2);
    else if(!strcmp(s1,"nu0_dust="))
      par->nu0_d=atof(s2);
    else if(!strcmp(s1,"beta_s_0="))
      par->beta_s_0=atof(s2);
    else if(!strcmp(s1,"beta_d_0="))
      par->beta_d_0=atof(s2);
    else if(!strcmp(s1,"temp_d_0="))
      par->temp_d_0=atof(s2);
    else if(!strcmp(s1,"beta_s_step="))
      par->beta_s_step=atof(s2);
    else if(!strcmp(s1,"beta_d_step="))
      par->beta_d_step=atof(s2);
    else if(!strcmp(s1,"temp_d_step="))
      par->temp_d_step=atof(s2);
    else if(!strcmp(s1,"beta_s_prior="))
      par->beta_s_prior=atof(s2);
    else if(!strcmp(s1,"beta_d_prior="))
      par->beta_d_prior=atof(s2);
    else if(!strcmp(s1,"temp_d_prior="))
      par->temp_d_prior=atof(s2);
    else if(!strcmp(s1,"nside="))
      par->nside=atoi(s2); 
    else if(!strcmp(s1,"nside_spec="))
      par->nside_spec=atoi(s2);
    else if(!strcmp(s1,"seed="))
      par->seed=atoi(s2);
    else if(!strcmp(s1,"n_samples="))
      par->n_samples=atoi(s2);
    else if(!strcmp(s1,"burning_fraction="))
      par->frac_samples_burn=atof(s2);
    else if(!strcmp(s1,"n_update_covar="))
      par->n_update_covar=atoi(s2);
    else if(!strcmp(s1,"n_spec_resample="))
      par->n_spec_resample=atoi(s2);
#ifdef _DEBUG
    else if(!strcmp(s1,"debug_pixel="))
      par->dbg_ipix=atoi(s2);
#endif //_DEBUG
   else
      fprintf(stderr,"BFoRe: Unknown parameter %s\n",s1);
  }
  fclose(fi);

  par->n_side_sub=par->nside/par->nside_spec;
  par->n_sub=par->n_side_sub*par->n_side_sub;
  par->n_pix=12*par->nside*par->nside;
  par->n_pix_spec=12*par->nside_spec*par->nside_spec;

  par->n_samples_burn=(int)(par->frac_samples_burn*par->n_samples);
  if(par->n_samples_burn>0.5*par->n_samples)
    report_error(1,"Wrong burning fraction %lE\n",par->frac_samples_burn);
  if(par->n_samples_burn<2*par->n_update_covar)
    report_error(1,"#burning samples should be larger than the covariance update period\n");

  if(par->flag_include_polarization==0)
    par->flag_independent_polarization=0;

  fi=my_fopen(par->fname_nulist,"r");
  par->n_nu=my_linecount(fi); rewind(fi);
  par->freqs=my_malloc(par->n_nu*sizeof(flouble));
  for(ii=0;ii<par->n_nu;ii++) {
    flouble dum;
#ifdef _SPREC
    int stat=fscanf(fi,"%f %f %f %f %f\n",&dum,&(par->freqs[ii]),&dum,&dum,&dum);
#else //_SPREC
    int stat=fscanf(fi,"%lf %lf %lf %lf %lf\n",&dum,&(par->freqs[ii]),&dum,&dum,&dum);
#endif //_SPREC
    if(stat!=5)
      report_error(1,"Error reading file %s, line %d\n",par->fname_nulist,ii+1);
  }
  fclose(fi);

  if(par->flag_include_polarization)
    par->n_pol=3;
  else
    par->n_pol=1;

  par->n_comp=0;
  if(par->flag_include_cmb) {
    par->index_cmb=par->n_comp;
    par->n_comp++;
  }
  if(par->flag_include_synchrotron) {
    par->index_synchrotron=par->n_comp;
    par->n_comp++;
  }
  if(par->flag_include_dust) {
    par->index_dust=par->n_comp;
    par->n_comp++;
  }

  par->n_param_max=3;
  if(par->flag_independent_polarization)
    par->n_param_max*=2;
    
  int index_novary=par->n_param_max;
  par->n_spec_vary=0;
  if(par->flag_include_synchrotron) {
    if(par->flag_beta_s_free)
      par->index_beta_s_t=par->n_spec_vary++;
    else
      par->index_beta_s_t=--index_novary;
  }
  if(par->flag_include_dust) {
    if(par->flag_beta_d_free)
      par->index_beta_d_t=par->n_spec_vary++;
    else
      par->index_beta_d_t=--index_novary;
    if(par->flag_temp_d_free)
      par->index_temp_d_t=par->n_spec_vary++;
    else
      par->index_temp_d_t=--index_novary;
  }
  if(par->flag_include_polarization && par->flag_independent_polarization) {
    if(par->flag_include_synchrotron) {
      if(par->flag_beta_s_free)
	par->index_beta_s_p=par->n_spec_vary++;
      else
	par->index_beta_s_p=--index_novary;
    }
    if(par->flag_include_dust) {
      if(par->flag_beta_d_free)
	par->index_beta_d_p=par->n_spec_vary++;
      else
	par->index_beta_d_p=--index_novary;
      if(par->flag_temp_d_free)
	par->index_temp_d_p=par->n_spec_vary++;
      else
	par->index_temp_d_p=--index_novary;
    }
  }
  else {
    par->index_beta_s_p=par->index_beta_s_t;
    par->index_beta_d_p=par->index_beta_d_t;
    par->index_temp_d_p=par->index_temp_d_t;
  }
  par->n_dof_pix=par->n_sub*par->n_pol*(par->n_nu-par->n_comp)-par->n_spec_vary;

  //Set prior
  par->prior_mean=my_malloc(par->n_param_max*sizeof(flouble));
  par->prior_isigma=my_malloc(par->n_param_max*sizeof(flouble));
  par->prior_mean[par->index_beta_s_t]=par->beta_s_0;
  par->prior_mean[par->index_beta_s_p]=par->beta_s_0;
  if(par->beta_s_prior<=0) {
    par->prior_isigma[par->index_beta_s_t]=-1;
    par->prior_isigma[par->index_beta_s_p]=-1;
  }
  else {
    par->prior_isigma[par->index_beta_s_t]=1./par->beta_s_prior;
    par->prior_isigma[par->index_beta_s_p]=1./par->beta_s_prior;
  }
  par->prior_mean[par->index_beta_d_t]=par->beta_d_0;
  par->prior_mean[par->index_beta_d_p]=par->beta_d_0;
  if(par->beta_d_prior<=0) {
    par->prior_isigma[par->index_beta_d_t]=-1;
    par->prior_isigma[par->index_beta_d_p]=-1;
  }
  else {
    par->prior_isigma[par->index_beta_d_t]=1./par->beta_d_prior;
    par->prior_isigma[par->index_beta_d_p]=1./par->beta_d_prior;
  }
  par->prior_mean[par->index_temp_d_t]=par->temp_d_0;
  par->prior_mean[par->index_temp_d_p]=par->temp_d_0;
  if(par->temp_d_prior<=0) {
    par->prior_isigma[par->index_temp_d_t]=-1;
    par->prior_isigma[par->index_temp_d_p]=-1;
  }
  else {
    par->prior_isigma[par->index_temp_d_t]=1./par->temp_d_prior;
    par->prior_isigma[par->index_temp_d_p]=1./par->temp_d_prior;
  }

  //Allocate maps
  if(NodeThis==0)
    printf("Allocating memory for maps\n");
  par->maps_data=my_malloc(par->n_pix*par->n_pol*par->n_nu*sizeof(flouble));
  par->maps_noise_weight=my_malloc(par->n_pix*par->n_pol*par->n_nu*sizeof(flouble));
  par->map_components_mean=my_calloc(par->n_pix*par->n_pol*par->n_comp,sizeof(flouble));
  par->map_components_covar=my_calloc(par->n_pix*par->n_pol*par->n_comp*par->n_comp,
				      sizeof(flouble));
  par->map_indices_mean=my_calloc(par->n_pix_spec*par->n_spec_vary,sizeof(flouble));
  par->map_indices_covar=my_calloc(par->n_pix_spec*par->n_spec_vary*par->n_spec_vary,
				   sizeof(flouble));
  par->map_chi2=my_calloc(par->n_pix,sizeof(flouble));
  par->dbg_extra=my_malloc(par->n_spec_vary*(par->n_samples+par->n_samples_burn)*
			   sizeof(flouble));

  //Read input maps
  if(NodeThis==0)
    printf("Reading data from %sXXX.fits\n",par->input_data_prefix);
  for(ii=0;ii<par->n_nu;ii++) {
    int jj;
    sprintf(fname_in,"%s%03d.fits",par->input_data_prefix,ii+1);
    for(jj=0;jj<par->n_pol;jj++) {
      int ip;
      map_dum=he_read_healpix_map(fname_in,&nside_dum,jj);
      if(nside_dum!=par->nside)
	report_error(1,"Read wrong nside\n");
      he_ring2nest_inplace(map_dum,nside_dum);
      for(ip=0;ip<par->n_pix;ip++)
	par->maps_data[ii+par->n_nu*(jj+par->n_pol*ip)]=map_dum[ip];
      free(map_dum);
    }
  }

  if(NodeThis==0)
    printf("Reading noise from %sXXX.fits\n",par->input_noise_prefix);
  flouble amin2_per_pix=4*M_PI*pow(180*60/M_PI,2)/par->n_pix;
  for(ii=0;ii<par->n_nu;ii++) {
    int jj;
    sprintf(fname_in,"%s%03d.fits",par->input_noise_prefix,ii+1);
    for(jj=0;jj<par->n_pol;jj++) {
      int ip;
      map_dum=he_read_healpix_map(fname_in,&nside_dum,jj);
      if(nside_dum!=par->nside)
	report_error(1,"Read wrong nside\n");
      he_ring2nest_inplace(map_dum,nside_dum);
      for(ip=0;ip<par->n_pix;ip++) {
	par->maps_noise_weight[ii+par->n_nu*(jj+par->n_pol*ip)]=
	  amin2_per_pix/(map_dum[ip]*map_dum[ip]);
      }
      free(map_dum);
    }
  }

  int npix_leftover,npix_pernode;

  npix_leftover=par->n_pix_spec%NNodes;
  npix_pernode=(par->n_pix_spec-npix_leftover)/NNodes;
  if(NodeThis<npix_leftover) {
    par->ipix_0=NodeThis*npix_pernode+NodeThis;
    par->ipix_f=par->ipix_0+npix_pernode+1;
  }
  else {
    par->ipix_0=NodeThis*npix_pernode+npix_leftover;
    par->ipix_f=par->ipix_0+npix_pernode;
  }
#ifdef _DEBUG
  printf("Node %d will compute pixels %d to %d\n",NodeThis,par->ipix_0,par->ipix_f-1);
#endif //_DEBUG

  //Print parameters
  if(NodeThis==0)
    param_bfore_print(par);

  return par;
}

#ifdef _DEBUG
static void write_debug_info(ParamBFoRe *par)
{
  FILE *fo;
  int size_float=sizeof(flouble);
  char fname[256];
  sprintf(fname,"%s_node%d_pix%d.dbg",par->output_prefix,NodeThis,par->dbg_ipix);
  fo=my_fopen(fname,"wb");

  my_fwrite(&size_float,sizeof(size_float),1,fo);
  my_fwrite(&(par->dbg_ipix),sizeof(par->dbg_ipix),1,fo);
  my_fwrite(&(par->nside),sizeof(par->nside),1,fo);
  my_fwrite(&(par->nside_spec),sizeof(par->nside_spec),1,fo);
  my_fwrite(&(par->n_sub),sizeof(par->n_sub),1,fo);
  my_fwrite(&(par->n_nu),sizeof(par->n_nu),1,fo);
  my_fwrite(&(par->n_pol),sizeof(par->n_pol),1,fo);
  my_fwrite(&(par->n_comp),sizeof(par->n_comp),1,fo);
  my_fwrite(&(par->n_spec_vary),sizeof(par->n_spec_vary),1,fo);
  my_fwrite(&(par->n_samples),sizeof(par->n_samples),1,fo);
  my_fwrite(&(par->n_samples_burn),sizeof(par->n_samples_burn),1,fo);
  
  //Write spectral chains
  my_fwrite(par->dbg_extra,sizeof(flouble),
	    (par->n_samples+par->n_samples_burn)*par->n_spec_vary,fo);
  
  //Write spectral mean
  my_fwrite(&(par->map_indices_mean[par->dbg_ipix*par->n_spec_vary]),sizeof(flouble),
	    par->n_spec_vary,fo);

  //Write spectral covar
  my_fwrite(&(par->map_indices_covar[par->dbg_ipix*par->n_spec_vary*par->n_spec_vary]),
	    sizeof(flouble),par->n_spec_vary*par->n_spec_vary,fo);

  //Write amplitudes mean
  my_fwrite(&(par->map_components_mean[par->dbg_ipix*par->n_comp*par->n_pol*par->n_sub]),
	    sizeof(flouble),par->n_comp*par->n_pol*par->n_sub,fo);
  
  //Write amplitudes covar
  my_fwrite(&(par->map_components_covar[par->dbg_ipix*par->n_comp*par->n_comp*
					par->n_pol*par->n_sub]),sizeof(flouble),
	    par->n_comp*par->n_comp*par->n_pol*par->n_sub,fo);

  //Write input data
  my_fwrite(&(par->maps_data[par->dbg_ipix*par->n_nu*par->n_pol*par->n_sub]),sizeof(flouble),
	    par->n_nu*par->n_pol*par->n_sub,fo);
  
  //Write input noise weights
  my_fwrite(&(par->maps_noise_weight[par->dbg_ipix*par->n_nu*par->n_pol*par->n_sub]),
	    sizeof(flouble),par->n_nu*par->n_pol*par->n_sub,fo);

  fclose(fo);
}
#endif //_DEBUG

void write_output(ParamBFoRe *par)
{
#ifdef _DEBUG
  write_debug_info(par);
#else //_DEBUG
  int ic1,ic2,ipol,ipix,ncorr,is1,is2,ispec;
  char fname[256];
  flouble **map_out;

#ifdef _WITH_MPI
  if(NodeThis==0) {
    MPI_Reduce(MPI_IN_PLACE,par->map_components_mean,
	       par->n_comp*par->n_pol+par->n_pix,
	       FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE,par->map_components_covar,
	       par->n_comp*par->n_comp*par->n_pol+par->n_pix,
	       FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE,par->map_indices_mean,
	       par->n_spec_vary*par->n_pix_spec,
	       FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE,par->map_indices_covar,
	       par->n_spec_vary*par->n_spec_vary*par->n_pix_spec,
	       FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  }
  else {
    MPI_Reduce(par->map_components_mean,NULL,
	       par->n_comp*par->n_pol+par->n_pix,
	       FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(par->map_components_covar,NULL,
	       par->n_comp*par->n_comp*par->n_pol+par->n_pix,
	       FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(par->map_indices_mean,NULL,
	       par->n_spec_vary*par->n_pix_spec,
	       FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(par->map_indices_covar,NULL,
	       par->n_spec_vary*par->n_spec_vary*par->n_pix_spec,
	       FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  }
#endif //_WITH_MPI

  if(NodeThis==0) {
    //Write output amplitude means
    sprintf(fname,"!%s_components_mean.fits",par->output_prefix);
    map_out=my_malloc(par->n_pol*par->n_comp*sizeof(flouble *));
    for(ipol=0;ipol<par->n_pol;ipol++) {
      for(ic1=0;ic1<par->n_comp;ic1++) {
	int imap=ic1+par->n_comp*ipol;
	map_out[imap]=my_malloc(par->n_pix*sizeof(flouble));
	for(ipix=0;ipix<par->n_pix;ipix++)
	  map_out[imap][ipix]=par->map_components_mean[ic1+par->n_comp*(ipol+par->n_pol*ipix)];
	he_nest2ring_inplace(map_out[imap],par->nside);
      }
    }
    he_write_healpix_map(map_out,par->n_pol*par->n_comp,par->nside,fname);
    for(ipol=0;ipol<par->n_pol;ipol++) {
      for(ic1=0;ic1<par->n_comp;ic1++) {
	int imap=ic1+par->n_comp*ipol;
	free(map_out[imap]);
      }
    }
    free(map_out);
    
    //Write output amplitude covariances
    sprintf(fname,"!%s_components_covar.fits",par->output_prefix);
    ncorr=par->n_comp*(par->n_comp+1)/2;
    map_out=my_malloc(par->n_pol*ncorr*sizeof(flouble *));
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int icov=0;
      for(ic1=0;ic1<par->n_comp;ic1++) {
	for(ic2=ic1;ic2<par->n_comp;ic2++) {
	  int imap=icov+ncorr*ipol;
	  map_out[imap]=my_malloc(par->n_pix*sizeof(flouble));
	  for(ipix=0;ipix<par->n_pix;ipix++) {
	    map_out[imap][ipix]=
	      par->map_components_mean[ic2+par->n_comp*(ic1+par->n_comp*
							(ipol+par->n_pol*ipix))];
	  }
	  he_nest2ring_inplace(map_out[imap],par->nside);
	  icov++;
	}
      }
    }
    he_write_healpix_map(map_out,par->n_pol*ncorr,par->nside,fname);
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int icov=0;
      for(ic1=0;ic1<par->n_comp;ic1++) {
	for(ic2=ic1;ic2<par->n_comp;ic2++) {
	  int imap=icov+ncorr*ipol;
	  free(map_out[imap]);
	  icov++;
	}
      }
    }
    free(map_out);
    
    //Write output spectral means
    sprintf(fname,"!%s_spec_mean.fits",par->output_prefix);
    map_out=my_malloc(par->n_spec_vary*sizeof(flouble *));
    for(is1=0;is1<par->n_spec_vary;is1++) {
      map_out[is1]=my_malloc(par->n_pix_spec*sizeof(flouble));
      for(ipix=0;ipix<par->n_pix_spec;ipix++)
	map_out[is1][ipix]=par->map_indices_mean[is1+par->n_spec_vary*ipix];
      he_nest2ring_inplace(map_out[is1],par->nside_spec);
    }
    he_write_healpix_map(map_out,par->n_spec_vary,par->nside_spec,fname);
    for(is1=0;is1<par->n_spec_vary;is1++)
      free(map_out[is1]);
    free(map_out);
    
    //Write output spectral covariances
    sprintf(fname,"!%s_spec_covar.fits",par->output_prefix);
    ncorr=par->n_spec_vary*(par->n_spec_vary+1)/2;
    map_out=my_malloc(ncorr*sizeof(flouble *));
    ispec=0;
    for(is1=0;is1<par->n_spec_vary;is1++) {
      for(is2=is1;is2<par->n_spec_vary;is2++) {
	map_out[ispec]=my_malloc(par->n_pix_spec*sizeof(flouble));
	for(ipix=0;ipix<par->n_pix_spec;ipix++)
	  map_out[ispec][ipix]=par->map_indices_covar[is2+par->n_spec_vary*
						      (is1+par->n_spec_vary*ipix)];
	he_nest2ring_inplace(map_out[ispec],par->nside_spec);
	ispec++;
      }
    }
    he_write_healpix_map(map_out,par->n_spec_vary,par->nside_spec,fname);
    ispec=0;
    for(is1=0;is1<par->n_spec_vary;is1++) {
      for(is2=is1;is2<par->n_spec_vary;is2++) {
	free(map_out[ispec]);
	ispec++;
      }
    }
    free(map_out);
  }
#endif //_DEBUG
}
