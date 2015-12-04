#include "common.h"

int main(int argc,char **argv)
{ 
  char fname_init[256];
  if(argc!=2) {
    printf("Usage: fg_rm.x param_file\n");
    exit(0);
  }
  sprintf(fname_init,"%s",argv[1]);

  gsl_set_error_handler_off();
  ParamBFoRe *par=read_params(fname_init);

#ifndef _DEBUG
#pragma omp parallel default(none) shared(par)
#endif //_DEBUG
  {
    int ipix_big;
#ifdef _DEBUG
    unsigned long seed_thr=par->seed;
#else //_DEBUG
    unsigned long seed_thr=par->seed+omp_get_thread_num();
#endif //_DEBUG
    PixelState *pst=pixel_state_new(par,seed_thr);

#ifdef _DEBUG
    ipix_big=par->dbg_ipix;
#else //_DEBUG
#pragma omp for
    for(ipix_big=0;ipix_big<par->n_pix_spec;ipix_big++)
#endif //_DEBUG
      {
	clean_pixel(par,pst,ipix_big);
      }//end omp for
    pixel_state_free(pst,par);
  }//end omp parallel

  printf("Writing output\n");
#ifdef _DEBUG
  write_debug_info(par);
#else //_DEBUG
  write_output(par);
#endif //_DEBUG
  param_bfore_free(par);

  return 0;
}
