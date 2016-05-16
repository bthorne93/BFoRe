#ifndef _COMMON_BFORE
#define _COMMON_BFORE

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#ifdef _WITH_OMP
#include <omp.h>
#endif //_WITH_OMP
#ifdef _WITH_MPI
#include <mpi.h>
#endif //_WITH_MPI
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>
#include <fitsio.h>
#include <chealpix.h>
#ifdef _WITH_SHT
#include <sharp_almhelpers.h>
#include <sharp_geomhelpers.h>
#include <sharp.h>
#ifdef _WITH_NEEDLET
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#endif //_WITH_NEEDLET
#endif //_WITH_SHT

#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
typedef int lint;
#endif //_LONGIDS

#ifdef _SPREC
typedef float flouble;
typedef float complex fcomplex;
#ifdef _WITH_MPI
#define FLOUBLE_MPI MPI_FLOAT
#endif //_WITH_MPI
#ifdef _WITH_SHT
#define SHT_TYPE 0
#endif //_WITH_SHT
#else //_SPREC
typedef double flouble;
typedef double complex fcomplex;
#ifdef _WITH_MPI
#define FLOUBLE_MPI MPI_DOUBLE
#endif //_WITH_MPI
#ifdef _WITH_SHT
#define SHT_TYPE SHARP_DP
#endif //_WITH_SHT
#endif //_SPREC

#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a,b)  (((a)<(b)) ? (a) : (b)) // minimum


extern int NNodes;
extern int NodeThis;
extern int IThread0;

typedef struct {
  int nside;
  int n_sub;
  int n_pix;

  char input_data_prefix[256];
  flouble *maps_data;

  char input_noise_prefix[256];
  flouble *maps_noise_weight;

  char output_prefix[256];

  char fname_nulist[256];
  int n_nu;
  flouble *freqs;

  int flag_include_polarization;
  int n_pol;

  int ipix_0;
  int ipix_f;

  int dbg_ipix;
  flouble *dbg_extra;

  //The following are only needed for Bayesian CS
  int do_bayes;
  int nside_spec;
  int n_side_sub;
  int n_pix_spec;

  char input_beta_s_t_prior[256];
  char input_beta_s_p_prior[256];
  char input_beta_d_t_prior[256];
  char input_beta_d_p_prior[256];
  char input_temp_d_t_prior[256];
  char input_temp_d_p_prior[256];
  flouble *map_prior_centres;
  flouble *map_prior_widths;

  int flag_write_samples;
  flouble *map_components_mean;
  flouble *map_components_covar;
  flouble *map_indices_mean;
  flouble *map_indices_covar;
  flouble *map_chi2;

  int flag_include_cmb;
  int flag_include_synchrotron;
  int flag_include_dust;
  int flag_include_volume_prior;
  int flag_use_marginal;
  int n_comp;
  int index_cmb;
  int index_synchrotron;
  int index_dust;

  int flag_independent_polarization;
  int flag_beta_s_free;
  int flag_beta_d_free;
  int flag_temp_d_free;
  int n_param_max;
  int n_spec_vary;
  int n_dof_pix;
  int index_beta_s_t;
  int index_beta_s_p;
  int index_beta_d_t;
  int index_beta_d_p;
  int index_temp_d_t;
  int index_temp_d_p;
  flouble beta_s_step;
  flouble beta_d_step;
  flouble temp_d_step;
  flouble nu0_s;
  flouble nu0_d;

  unsigned long seed;
  int n_samples;
  int n_output_rate;
  flouble frac_samples_burn;
  int n_update_covar;
  int n_samples_burn;
  int n_spec_resample;

} ParamBFoRe;

//Defined in common.c
int my_linecount(FILE *f);
void report_error(int level,char *fmt,...);
void *my_malloc(size_t size);
void *my_calloc(size_t nmemb,size_t size);
FILE *my_fopen(const char *path,const char *mode);
size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream);
size_t my_fread(void *ptr,size_t size,size_t count,FILE *stream);
void param_bfore_free(ParamBFoRe *par);
ParamBFoRe *read_params(char *fname);
void write_output(ParamBFoRe *par);
void dbg_printf(int do_print,char *fmt,...);

//Defined in healpix_extra.c
//HE_IO
void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname);
flouble *he_read_healpix_map(char *fname,long *nside,int nfield);
//HE_PIX
int he_ring_num(long nside,double z);
long *he_query_strip(long nside,double theta1,double theta2,long *npix_strip);
void he_ring2nest_inplace(flouble *map_in,long nside);
void he_nest2ring_inplace(flouble *map_in,long nside);
void he_udgrade(flouble *map_in,long nside_in,flouble *map_out,long nside_out,int nest);
//HE_SHT
#ifdef _WITH_SHT
#define HE_MAX_SHT 32
#define HE_FWHM2SIGMA 0.00012352884853326381 //Transforms FWHM in arcmin to sigma_G in rad:
long he_nalms(int lmax);
long he_indexlm(int l,int m,int lmax);
void he_alm2map(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
void he_map2alm(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
void he_alm2cl(fcomplex **alms_1,fcomplex **alms_2,
	       int nmaps_1,int nmaps_2,int pol_1,int pol_2,flouble **cls,int lmax);
void he_anafast(flouble **maps_1,flouble **maps_2,
		int nmaps_1,int nmaps_2,int pol_1,int pol_2,
		flouble **cls,int nside,int lmax);
flouble *he_generate_beam_window(int lmax,flouble fwhm_amin);
void he_alter_alm(int lmax,flouble fwhm_amin,fcomplex *alm_in,
		  fcomplex *alm_out,flouble *window);
flouble *he_synfast(flouble *cl,int nside,int lmax,unsigned int seed);
//HE_NT
#ifdef _WITH_NEEDLET
#define HE_NBAND_NX 512
#define HE_NORM_FT 2.2522836206907617
#define HE_NL_INTPREC 1E-6
#define HE_NT_NSIDE_MIN 32
typedef struct {
  double b;
  double inv_b;
  gsl_spline *b_spline;
  gsl_interp_accel *b_intacc;
  int nside0;
  int nj;
  int *nside_arr;
  int *lmax_arr;
  flouble **b_arr;
} HE_nt_param;
void he_nt_end(HE_nt_param *par);
HE_nt_param *he_nt_init(flouble b_nt,int nside0);
flouble ***he_alloc_needlet(HE_nt_param *par,int pol);
void he_nt_get_window(HE_nt_param *par,int j,flouble *b);
fcomplex **he_map2needlet(HE_nt_param *par,flouble **map,int pol,flouble ***nt,int return_alm);
fcomplex **he_needlet2map(HE_nt_param *par,flouble **map,int pol,flouble ***nt,int return_alm);
#endif //_WITH_NEEDLET
#endif //_WITH_SHT

//Defined in rng.c
#define RNG_NRAN 624
#define RNG_MRAN 397
#define RNG_MATRIX_A 0x9908b0df
#define RNG_UPPER_MASK 0x80000000UL
#define RNG_LOWER_MASK 0x7fffffffUL
typedef struct {
  unsigned long mt[RNG_NRAN];
  int mti;
  int calc_gauss;
  double u;
  double phi;
} Rng;
Rng *init_rng(unsigned long seed);
void end_rng(Rng *rng);
unsigned long rand_ulong(Rng *rng);
double rand_real01(Rng *rng);
double rand_gauss(Rng *rng);

//Defined in powell.c
#define TINY 1.0E-25
#define ZEPS 1.0E-10
#define GLIMIT 100.0
#define GOLD 1.618034
#define CGOLD 0.3819660
#define LIN_TOL 2.0E-4
#define LIN_ITMAX 100
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);  

typedef struct {
  int n;
  double *p;
  double **xi;
  int iter;
  int max_iter;
  double fret;
  double ftol;
  double *xdum;
  double *xdir;
  double (*fun)(double *,void *);
  void *params;
} PowellParams;
void free_powell_params(PowellParams *par);
PowellParams *powell_params_new(int n,flouble *p,flouble (*fun)(flouble *,void *),
				void *params,int max_iter,flouble ftol);
void powell(PowellParams *par);

//Defined in bfore.c
typedef struct {
  flouble *f_matrix;
  gsl_matrix **cov_inv;
  gsl_vector **vec_mean;
  flouble *prior_mean;
  flouble *prior_isigma;
  flouble *rand_spec;
  gsl_vector *vaux;
  double chi2;
} PixelState;
PixelState *pixel_state_new(ParamBFoRe *par);
void pixel_state_free(PixelState *pst,ParamBFoRe *par);
void clean_pixel(ParamBFoRe *par,Rng *rng,PixelState *pst,int ipix_big);
void clean_pixel_from_marginal(ParamBFoRe *par,Rng *rng,PixelState *pst_old,
			       PixelState *pst_new,int ipix_big);

#endif //_COMMON_BFORE
