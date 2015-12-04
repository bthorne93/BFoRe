#include "common.h"

void mpi_init(int *p_argc,char ***p_argv)
{
#ifdef _WITH_MPI
  int ii,nthreads_this;
  int *nthreads_all;
  MPI_Init(p_argc,p_argv);

  MPI_Comm_size(MPI_COMM_WORLD,&NNodes);
  MPI_Comm_rank(MPI_COMM_WORLD,&NodeThis);

  nthreads_all=my_malloc(NNodes*sizeof(int));
#ifdef _WITH_OMP
  nthreads_this=omp_get_num_threads();
#else //_WITH_OMP
  nthreads_this=1;
#endif //_WITH_OMP
  MPI_Allgather(&nthreads_this,1,MPI_INT,nthreads_all,1,MPI_INT,MPI_COMM_WORLD);
  if(NodeThis==0) {
    for(ii=0;ii<NNodes;ii++)
      printf("Node %d has %d threads\n",ii,nthreads_all[ii]);
  }
  IThread0=0;
  for(ii=0;ii<NodeThis;ii++)
    IThread0+=nthreads_all[ii];
#ifdef _DEBUG
  printf("Node %d, thread count starts at %d\n",NodeThis,IThread0);
#endif //_DEBUG
  free(nthreads_all);

#else //_WITH_MPI

  NNodes=1;
  NodeThis=0;
  IThread0=0;
#endif //_WITH_MPI
}

int main(int argc,char **argv)
{ 
  char fname_init[256];
  if(argc!=2) {
    printf("Usage: fg_rm.x param_file\n");
    exit(0);
  }
  sprintf(fname_init,"%s",argv[1]);

  mpi_init(&argc,&argv);
  gsl_set_error_handler_off();
  ParamBFoRe *par=read_params(fname_init);

#ifdef _WITH_OMP
#ifndef _DEBUG
#pragma omp parallel default(none) shared(par,IThread0)
#endif //_DEBUG
#endif //_WITH_OMP
  {
    int ipix_big,ithr;
    unsigned long seed_thr;
    PixelState *pst;

    ithr=0;
#ifdef _WITH_OMP
#ifndef _DEBUG
    ithr=omp_get_thread_num();
#endif //_DEBUG
#endif //_WITH_OMP
    seed_thr=par->seed+IThread0+ithr;
    pst=pixel_state_new(par,seed_thr);

#ifdef _DEBUG
    ipix_big=par->dbg_ipix;
#else //_DEBUG
#ifdef _WITH_OMP
#pragma omp for
#endif //_WITH_OMP
    for(ipix_big=par->ipix_0;ipix_big<par->ipix_f;ipix_big++)
#endif //_DEBUG
      {
	clean_pixel(par,pst,ipix_big);
      }//end omp for
    pixel_state_free(pst,par);
  }//end omp parallel

  printf("Writing output\n");
  write_output(par);
  param_bfore_free(par);

#ifdef _WITH_MPI
  MPI_Finalize();
#endif //_WITH_MPI

  return 0;
}
