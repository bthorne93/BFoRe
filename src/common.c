#include "common.h"

int my_linecount(FILE *f)
{
  int i0=0;
  char ch[1000];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
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

static ParamFGRM *param_fgrm_new(void)
{
  ParamFGRM *par=(ParamFGRM *)my_malloc(sizeof(ParamFGRM));

  par->nside=256;
  par->nside_spec=32;
  par->n_side_sub=8; //
  par->n_sub=64; //
  par->n_pix=786432; //
  
  sprintf(par->input_data_prefix,"default");
  par->maps_data=NULL; //

  sprintf(par->input_noise_prefix,"default");
  par->maps_noise_weight=NULL; //

  sprintf(par->output_prefix,"default");
  par->map_components_mean=NULL; //
  par->map_components_covar=NULL; //
  par->map_indices_mean=NULL; //
  par->map_indices_covar=NULL; //
  par->map_chi2=NULL; //

  sprintf(par->fname_nulist,"default");
  par->n_nu=-1; //
  par->freqs=NULL; //
  
  par->flag_include_polarization=0;
  par->n_pol=1; //

  par->flag_include_cmb=0;
  par->flag_include_synchrotron=0;
  par->flag_include_dust=0;
  par->n_comp=0; //
  par->index_cmb=-1; //
  par->index_synchrotron=-1; //
  par->index_dust=-1; //

  par->flag_independent_polarization=0;
  par->flag_beta_s_free=0;
  par->flag_beta_d_free=0;
  par->flag_temp_d_free=0;
  par->n_param_max=0; //
  par->n_spec_vary=0; //
  par->n_dof_pix=0; //
  par->index_beta_s_t=-1; //
  par->index_beta_s_p=-1; //
  par->index_beta_d_t=-1; //
  par->index_beta_d_p=-1; //
  par->index_temp_d_t=-1; //
  par->index_temp_d_p=-1; //
  
  par->beta_s_0=-1.;
  par->beta_d_0=1.54;
  par->temp_d_0=20.;
  par->sigma_beta_s=0.1;
  par->sigma_beta_d=0.1;
  par->sigma_temp_d=0.1;
  par->nu0_s=23.;
  par->nu0_d=353.;

  return par;
}

static void param_fgrm_print(ParamFGRM *par)
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
    printf(" - Will only include\n");
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
  printf(" D(beta_s) = %.3lf\n",par->sigma_beta_s);
  printf(" - beta_d_0 = %.3lf, ",par->beta_d_0);
  printf(" D(beta_d) = %.3lf\n",par->sigma_beta_d);
  printf(" - temp_d_0 = %.3lf, ",par->temp_d_0);
  printf(" D(temp_d) = %.3lf\n",par->sigma_temp_d);
}

void param_fgrm_free(ParamFGRM *par)
{
  free(par->maps_data);
  free(par->maps_noise_weight);
  free(par->map_components_mean);
  free(par->map_components_covar);
  free(par->map_indices_mean);
  free(par->map_indices_covar);
  free(par->map_chi2);
  free(par->freqs);
  free(par);
}

ParamFGRM *read_params(char *fname)
{
  FILE *fi;
  char fname_in[256];
  int n_lin,ii;
  long nside_dum;
  ParamFGRM *par=param_fgrm_new();
  flouble *map_dum;

  //Read parameters from file
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
    else if(!strcmp(s1,"sigma_beta_s="))
      par->sigma_beta_s=atof(s2);
    else if(!strcmp(s1,"sigma_beta_d="))
      par->sigma_beta_d=atof(s2);
    else if(!strcmp(s1,"sigma_temp_d="))
      par->sigma_temp_d=atof(s2);
    else if(!strcmp(s1,"nside="))
      par->nside=atoi(s2); 
    else if(!strcmp(s1,"nside_spec="))
      par->nside_spec=atoi(s2);
   else
      fprintf(stderr,"FGRM: Unknown parameter %s\n",s1);
  }
  fclose(fi);

  par->n_side_sub=par->nside/par->nside_spec;
  par->n_sub=par->n_side_sub*par->n_side_sub;
  par->n_pix=12*par->nside*par->nside;
  par->n_pix_spec=12*par->nside_spec*par->nside_spec;

  if(par->flag_include_polarization==0)
    par->flag_independent_polarization=0;

  fi=my_fopen(par->fname_nulist,"r");
  par->n_nu=my_linecount(fi); rewind(fi);
  par->freqs=my_malloc(par->n_nu*sizeof(flouble));
  for(ii=0;ii<par->n_nu;ii++) {
    flouble dum;
    int stat=fscanf(fi,"%lf %lf %lf %lf %lf\n",&dum,&(par->freqs[ii]),&dum,&dum,&dum);
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

  //Allocate maps
  printf("Allocating memory for maps\n");
  par->maps_data=my_malloc(par->n_pix*par->n_pol*par->n_nu*sizeof(flouble));
  par->maps_noise_weight=my_malloc(par->n_pix*par->n_pol*par->n_nu*sizeof(flouble));
  par->map_components_mean=my_calloc(par->n_pix*par->n_pol*par->n_comp,sizeof(flouble));
  par->map_components_covar=my_calloc(par->n_pix*par->n_pol*par->n_comp*par->n_comp,sizeof(flouble));
  par->map_indices_mean=my_calloc(par->n_pix_spec*par->n_spec_vary,sizeof(flouble));
  par->map_indices_covar=my_calloc(par->n_pix_spec*par->n_spec_vary*par->n_spec_vary,sizeof(flouble));
  par->map_chi2=my_calloc(par->n_pix,sizeof(flouble));

  //Read input maps
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
      for(ip=0;ip<par->n_pix;ip++)
	par->maps_noise_weight[ii+par->n_nu*(jj+par->n_pol*ip)]=amin2_per_pix/(map_dum[ip]*map_dum[ip]);
      free(map_dum);
    }
  }

  //Print parameters
  param_fgrm_print(par);

  return par;
}
