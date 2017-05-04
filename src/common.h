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
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_deriv.h>
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

//Defined in healpix_extra.c
//HE_IO
void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname);
flouble *he_read_healpix_map(char *fname,long *nside,int nfield);
void he_get_file_params(char *fname,long *nside,int *nfields,int *isnest);
//HE_PIX
int he_ring_num(long nside,double z);
void he_query_strip(long nside,double theta1,double theta2,int *pixlist,long *npix_strip);
void he_ring2nest_inplace(flouble *map_in,long nside);
void he_nest2ring_inplace(flouble *map_in,long nside);
void he_udgrade(flouble *map_in,long nside_in,flouble *map_out,long nside_out,int nest);
long he_nside2npix(long nside);
void he_pix2vec_ring(long nside, long ipix, double *vec);
long he_ang2pix(long nside,double cth,double phi);
void he_in_ring(int nside,int iz,flouble phi0,flouble dphi,int *listir,int *nir);
void he_query_disc(int nside,double cth0,double phi,flouble radius,int *listtot,int *nlist,int inclusive);
//HE_SHT
#ifdef _WITH_SHT
#define HE_MAX_SHT 32
#define HE_FWHM2SIGMA 0.00012352884853326381 //Transforms FWHM in arcmin to sigma_G in rad:
long he_nalms(int lmax);
long he_indexlm(int l,int m,int lmax);
void he_alm2map(int nside,int lmax,int ntrans,int spin,flouble **maps,fcomplex **alms);
void he_map2alm(int nside,int lmax,int ntrans,int spin,flouble **maps,fcomplex **alms,int niter);
void he_alm2cl(fcomplex **alms_1,fcomplex **alms_2,int pol_1,int pol_2,flouble **cls,int lmax);
void he_anafast(flouble **maps_1,flouble **maps_2,int pol_1,int pol_2,flouble **cls,
		int nside,int lmax,int iter);
flouble *he_generate_beam_window(int lmax,double fwhm_amin);
void he_zero_alm(int lmax,fcomplex *alm);
void he_alter_alm(int lmax,double fwhm_amin,fcomplex *alm_in,fcomplex *alm_out,flouble *window,int add_to_out);
void he_map_product(int nside,flouble *mp1,flouble *mp2,flouble *mp_out);
flouble he_map_dot(int nside,flouble *mp1,flouble *mp2);
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
  int niter;
  int nside0;
  int jmax_min;
  int nj;
  int *nside_arr;
  int *lmax_arr;
  flouble **b_arr;
} he_needlet_params;
void he_nt_end(he_needlet_params *par);
he_needlet_params *he_nt_init(flouble b_nt,int nside0,int niter);
void he_free_needlet(he_needlet_params *par,int pol,flouble ***nt);
flouble ***he_alloc_needlet(he_needlet_params *par,int pol);
void he_nt_get_window(he_needlet_params *par,int j,flouble *b);
fcomplex **he_needlet2map(he_needlet_params *par,flouble **map,flouble ***nt,
			  int return_alm,int pol,int input_TEB,int output_TEB);
fcomplex **he_map2needlet(he_needlet_params *par,flouble **map,flouble ***nt,
			  int return_alm,int pol,int input_TEB,int output_TEB);
#endif //_WITH_NEEDLET
#endif //_WITH_SHT

typedef struct {
  int nside; //Map nside
  int n_pix; //Number of pixels per map

  char input_data_prefix[256]; //Prefix pointing to the data maps
  flouble *maps_data; //Data frequency maps

  char output_prefix[256]; //Output prefix

  char fname_nulist[256]; //File containing frequencies
  int n_nu; //Number of frequencies
  flouble *freqs; //Frequencies

  int flag_include_polarization; //Do we have polarization?
  int n_pol; //Number of polarization channels (1-> T, 3-> T,Q,U)

  int ipix_0; //First pixel corresponding to this node
  int ipix_f; //Last pixel for this node (this one actually corresponds to the next node)

  flouble dec_dbg_ipix;
  flouble ra_dbg_ipix;
  long dbg_ipix; //Pixel index used for debugging
  flouble *dbg_extra; //Debugging data

  //The following are only needed for NILC
  int do_nilc;
  flouble b_nt;
  int nilc_ndim;
  int nilc_niter;
  flouble ***map_nilc_in;
  flouble **map_nilc_out;

  //The following are only needed for Bayesian CS
  int do_bayes;
  int nside_spec; //Nside for spectral indices
  int n_side_sub; //Number of sub_pixels per side in each spectral index pixel
  int n_sub; //Number of sub_pixels in each spectral index pixel
  int n_pix_spec; //Number of spectral index pixels
  int n_pix_spec_unmasked; //Number of unmasked spectral index pixels

  char input_noise_prefix[256]; //Prefix pointing to the noixe variance maps
  flouble *maps_noise_weight; //Noise weight maps (basically 1/noise_variance_per_pixel)

  char input_mask_fname[256]; //Mask filename
  int *ipix_unmasked; //Indices of unmasked pixels

  char input_beta_s_t_prior[256]; //Path to prior map on beta_s (temperature)
  char input_beta_s_p_prior[256]; //Path to prior map on beta_s (polarization)
  char input_beta_d_t_prior[256]; //Path to prior map on beta_d (temperature)
  char input_beta_d_p_prior[256]; //Path to prior map on beta_d (polarization)
  char input_temp_d_t_prior[256]; //Path to prior map on temp_d (temperature)
  char input_temp_d_p_prior[256]; //Path to prior map on temp_d (polarization)
  flouble *map_prior_centres; //Maps of the spectral index prior mean
  flouble *map_prior_widths; //Maps of the spectral index prior width

  int flag_write_samples; //Do we want to output samples?
  flouble *map_components_mean; //Mean of amplitudes
  flouble *map_components_covar; //Covariance of amplitudes
  flouble *map_indices_mean; //Mean of spectral indices
  flouble *map_indices_covar; //Covariance of spectral indices
  flouble *map_chi2; //Chi^2 map

  int flag_include_cmb; //Include CMB in sky model?
  int flag_include_synchrotron; //Incude synchrotron in sky model?
  int flag_include_dust; //Include dust in sky model?
  int flag_include_volume_prior; //Use volume (Jeffeys) prior?
  int flag_use_marginal; //Sample spectral indices from marginal distribution?
  int n_comp; //Number of components (up to 3)
  int index_cmb; //Index for CMB component
  int index_synchrotron; //Index for synchrotron component
  int index_dust; //Index for dust component

  int flag_independent_polarization; //Assume independent spectral indices in polarization?
  int flag_beta_s_free; //Is beta_s free?
  int flag_beta_d_free; //Is beta_d free?
  int flag_temp_d_free; //Is temp_d free?
  int n_param_max; //Maximum number of parameters to sample
  int n_spec_vary; //Number of free spectral indices
  int n_dof_pix; //Number of degrees of freedom per spectral index pixel
  int index_beta_s_t; //Index for beta_s (temperature)
  int index_beta_s_p; //Index for beta_s (polarization)
  int index_beta_d_t; //Index for beta_d (temperature)
  int index_beta_d_p; //Index for beta_d (polarization)
  int index_temp_d_t; //Index for temp_d (temperature)
  int index_temp_d_p; //Index for temp_d (polarization)
  flouble beta_s_step; //Initial step size in beta_s
  flouble beta_d_step; //Initial step size in beta_d
  flouble temp_d_step; //Initial step size in temp_d
  flouble nu0_s; //Reference frequency for synchrotron
  flouble nu0_d; //Reference frequency for dust

  unsigned long seed; //Seed
  int n_samples; //Total number of samples per spectral index pixel
  int n_output_rate; //Output sample every so many samples
  flouble frac_samples_burn; //Fraction of samples used for burning
  int n_samples_burn; //Number of samples used for burning
  int n_update_covar; //Number of samples used to compute the initial covariance matrix
  int n_spec_resample; //Number of spectral index samples takeng for each amplitudes sample

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
flouble freq_evolve(int spec_type,double nu_0,double beta,double temp,double nu);

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
  flouble *f_matrix; //Frequency evolution matrix
  gsl_matrix **cov_inv; //Set of covariance matrices for component amplitudes (one per pixel)
  gsl_vector **vec_mean; //Set of vectors of mean component amplitudes (one per pixel)
  flouble *prior_mean; //Prior mean (one per spectral param)
  flouble *prior_isigma; //1/sigma of prior (one per spectral param)
  flouble *rand_spec; //Dummy array to fill with random numbers (one element per spectral param)
  gsl_vector *vaux; //Dummy vector (one element per component)
  double chi2; //Pixel chi2
} PixelState;
PixelState *pixel_state_new(ParamBFoRe *par);
void pixel_state_free(PixelState *pst,ParamBFoRe *par);
void clean_pixel(ParamBFoRe *par,Rng *rng,PixelState *pst,int ipix_big);
void clean_pixel_from_marginal(ParamBFoRe *par,Rng *rng,PixelState *pst_old,
			       PixelState *pst_new,int ipix_big);

//Defined in nilc.c
void do_nilc(ParamBFoRe *par);

#endif //_COMMON_BFORE
