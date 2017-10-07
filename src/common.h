/** @file common.h
 *  @brief Includes for all the required libraries and declarations of all the
 *  functions used in BFoRe.
 *
 *  This file includes all the external libraries and declares all the
 *  functions used by BFoRe.
 *
 *  @author David Alonso (Primary author)
 *  @author Ben Thorne (commenter)
 *  @bug No known bugs.
 */

#ifndef _COMMON_BFORE
  #define _COMMON_BFORE

  #include <stdio.h>
  #include <stdlib.h>
  #include <stdarg.h>
  #include <string.h>
  #include <math.h>
  #include <time.h>
  #include <complex.h>

  #include <gsl/gsl_vector.h>
  #include <gsl/gsl_matrix.h>
  #include <gsl/gsl_linalg.h>
  #include <gsl/gsl_blas.h>
  #include <gsl/gsl_multifit.h>
  #include <gsl/gsl_multimin.h>

  #ifdef _WITH_OMP
    #include <omp.h>
  #endif //_WITH_OMP

  #ifdef _WITH_MPI
    #include <mpi.h>
  #endif //_WITH_MPI

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
  #else //_SPREC
    typedef double flouble;
    typedef double complex fcomplex;
    #ifdef _WITH_MPI
      #define FLOUBLE_MPI MPI_DOUBLE
    #endif //_WITH_MPI
  #endif //_SPREC

  extern int NNodes;
  extern int NodeThis;
  extern int IThread0;

  typedef struct {
    // nside on which amplitudes vary
    int nside;
    // nside on which spectral parameter vary
    int nside_spec;
    //
    int n_side_sub;
    //  square of the ratio between nsides
    int n_sub;
    // number of pixels at base resolution
    int n_pix;
    // number of pixels at spectral parameter resolution
    int n_pix_spec;

    // prefix for the input data. This should be a path to a file
    // with the ending %03d.fits, but not including that ending.
    char input_data_prefix[256];
    // data vector containing the input maps.
    flouble *maps_data;

    // prefix for the input noise maps.
    char input_noise_prefix[256];
    // data vector containing the noise level maps.
    flouble *maps_noise_weight;

    char input_beta_s_t_prior[256];
    char input_beta_s_p_prior[256];
    char input_curv_s_t_prior[256];
    char input_curv_s_p_prior[256];
    char input_beta_d_t_prior[256];
    char input_beta_d_p_prior[256];
    char input_temp_d_t_prior[256];
    char input_temp_d_p_prior[256];
    flouble *map_prior_centres;
    flouble *map_prior_widths;

    char output_prefix[256];
    int flag_write_samples;
    flouble *map_components_mean;
    flouble *map_components_covar;
    flouble *map_indices_mean;
    flouble *map_indices_covar;
    flouble *map_chi2;

    char fname_nulist[256];
    int n_nu;
    flouble *freqs;

    int flag_include_polarization;
    int n_pol;

    int flag_include_cmb;
    int flag_include_synchrotron;
    int flag_include_curvature;
    int flag_include_dust;
    int flag_include_volume_prior;
    int flag_use_marginal;
    int n_comp;
    int index_cmb;
    int index_synchrotron;
    int index_dust;

    int flag_independent_polarization;
    int flag_beta_s_free;
    int flag_curv_s_free;
    int flag_beta_d_free;
    int flag_temp_d_free;
    int n_param_max;
    // number of spectral parameters to be fitted.
    int n_spec_vary;
    int n_dof_pix;
    int index_beta_s_t;
    int index_beta_s_p;
    int index_curv_s_t;
    int index_curv_s_p;
    int index_beta_d_t;
    int index_beta_d_p;
    int index_temp_d_t;
    int index_temp_d_p;
    flouble beta_s_step;
    flouble curv_s_step;
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

    // first and last large pixel that an individual process is responsible
    // for cleaning.
    int ipix_0;
    int ipix_f;

    int dbg_ipix;
    flouble *dbg_extra;
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
  #ifdef _WITH_SHT
  long he_nalms(int lmax);
  long he_indexlm(int l,int m,int lmax);
  void he_alm2map(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
  void he_map2alm(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
  void he_anafast(flouble **maps_1,flouble **maps_2,
  		int nmaps_1,int nmaps_2,
  		int pol_1,int pol_2,
  		flouble **cls,int nside,int lmax);
  double *he_generate_beam_window(int lmax,double fwhm_amin);
  void he_alter_alm(int lmax,double fwhm_amin,fcomplex *alms,double *window);
  #endif //_WITH_SHT
  void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname);
  flouble *he_read_healpix_map(char *fname,long *nside,int nfield);
  int he_ring_num(long nside,double z);
  void he_ring2nest_inplace(flouble *map_in,long nside);
  void he_nest2ring_inplace(flouble *map_in,long nside);
  void he_udgrade(flouble *map_in,long nside_in,
  		flouble *map_out,long nside_out,
  		int nest);

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
  #define MAX(a, b)  (((a) > (b)) ? (a) : (b))
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

  /** @brief function to allocate memory for a PixelState struct.
   *  @param par an instance of the ParamBFoRe struct with all the information
   *         for the run
   *  @return PixelState
   */
  PixelState *pixel_state_new(ParamBFoRe *par);

  /** @brief function to free memory allocated to a PixelState struct.
   *  @param pst instance of PixelState to be freed.
   *  @param par instance of ParamBFoRe that was used to initialize pst.
   *  @return void
   */
  void pixel_state_free(PixelState *pst,ParamBFoRe *par);

  /** @brief clean pixel using old method.
   *  @param par
   *  @param rng
   *  @param pst
   *  @ipix_big
   *  @param par instance of ParamBFoRe that was used to initialize pst.
   *  @return void
   */
  void clean_pixel(ParamBFoRe *par,Rng *rng,PixelState *pst,int ipix_big);

  /** @brief clean a large pixel (ipix_big) having marginalized over amplitudes.
   *  @param rng seeded random number generator.
   *  @param pst initial instance of PixelState which contains the priors.
   *  @return void
   */
  void clean_pixel_from_marginal(ParamBFoRe *par,Rng *rng,PixelState *pst_old,
  			       PixelState *pst_new,int ipix_big);

#endif //_COMMON_BFORE
