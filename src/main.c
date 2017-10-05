/** @file main.c
 *  @brief Main file of BFoRe.
 *
 *  This contains the main function of BFoRe controlling the parallelization
 *  and implementation of the code.
 *
 *  @author David Alonso (Primary author)
 *  @author Ben Thorne (commenter)
 *  @bug No known bugs.
 */


#include "common.h"

void mpi_init(int *p_argc, char ***p_argv)
{
#ifdef _WITH_MPI
  int ii, nthreads_this;
  int *nthreads_all;

  /* Initialize MPI processes and get number of nodes, NNodes, and individual
  node number, NodeThis.
  */
  MPI_Init(p_argc, p_argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NNodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &NodeThis);

  nthreads_all=my_malloc(NNodes*sizeof(int));

  /* If compiling with openmp each process checks for number of available
  threads. Indent is due to this being nested within _WITH_MPI
  */
  #ifdef _WITH_OMP
  nthreads_this=omp_get_max_threads();
  #else //_WITH_OMP
  nthreads_this=1;
  #endif //_WITH_OMP

  /* Gather information about the number of threads available to each node.
  Sneds each node's nthreads_this to all other nodes. Therefore every node
  knows how many threads each other node has, these are stored in
  nthreads_all, which is NNodes integers long.

  Master node then prints summary of this information.
  */
  MPI_Allgather(&nthreads_this, 1, MPI_INT, nthreads_all, 1, MPI_INT, MPI_COMM_WORLD);
  if(NodeThis == 0)
  {
    for(ii = 0; ii < NNodes; ii++)
    {
      printf("Node %d has %d threads\n", ii, nthreads_all[ii]);
    }
  }

  /* IThread0 is built up by summing number of threads for all nodes below
  the current node.
  */
  IThread0 = 0;
  for(ii = 0; ii < NodeThis; ii++)
  {
    IThread0 += nthreads_all[ii];
  }

  #ifdef _DEBUG
  printf("Node %d, thread count starts at %d \n", NodeThis, IThread0);
  #endif //_DEBUG
  free(nthreads_all);
  #else //_WITH_MPI

  NNodes = 1;
  NodeThis = 0;
  IThread0 = 0;
#endif //_WITH_MPI
}

int main(int argc, char **argv)
{
  /* Main function of BFoRe. BFoRe expects only one argument: pat
  */
  char fname_init[256];

  if(argc! = 2)
  {
    // Check command line inputs.
    printf("Usage: bfore.x param_file\n");
    exit(0);
  }

  // Store parameter file path into fname_init
  sprintf(fname_init, "%s", argv[1]);

  // Initialize MPI environment.
  mpi_init(&argc, &argv);
  gsl_set_error_handler_off();

  // Read in parameter file to ParamBFoRe struct defined in common.h
  ParamBFoRe *par = read_params(fname_init);

  // If compiled with openmp then each process checks the number of available
  // threads.
  int ii, n_threads;
#ifdef _WITH_OMP
  n_threads = omp_get_max_threads();
#else //_WITH_OMP
  n_threads=1;
#endif //_WITH_OMP

  // Pixel state is initialized. This is defined in common.h
  PixelState **pst_old, **pst_new=NULL;

  pst_old = my_malloc(n_threads * sizeof(PixelState *));
  for (ii = 0; ii < n_threads; ii++)
  {
    pst_old[ii] = pixel_state_new(par);
  }
  // This code seems to be a copy of the previous few lines. Why not just
  // copy pst_old?
  if (par -> flag_use_marginal)
  {
    pst_new = my_malloc(n_threads * sizeof(PixelState *));
    for (ii = 0; ii < n_threads; ii++)
    {
      pst_new[ii] = pixel_state_new(par);
    }
  }

// If compiling with openmp and not debugging a single pixel loop through
// using openmp. This distributes jobs across threads within a particular node.
#ifdef _WITH_OMP
  #ifndef _DEBUG_SINGLEPIX
    #pragma omp parallel default(none) shared(par, IThread0, NodeThis, pst_old, pst_new)
  #endif //_DEBUG_SINGLEPIX
#endif //_WITH_OMP
  {
    int ipix_big, ithr;
    unsigned long seed_thr;
    Rng *rng;
    ithr = 0;

#ifdef _WITH_OMP
  #ifndef _DEBUG_SINGLEPIX
    // Get the number of the thread we are in. This is a thread within a node.
    ithr = omp_get_thread_num();
  #endif //_DEBUG_SINGLEPIX
#endif //_WITH_OMP
    // seed the thread we are on. IThread0 counts all the threads on the nodes
    // below us. Then ithr takes us to the current thread. Therefore each
    // thread in the entire program gets a different seed and all seeds are
    // sequential. Then seed the rng.
    seed_thr = par -> seed + IThread0 + ithr;
    rng = init_rng(seed_thr);

// We now loop over the large pixels over which spectral parameters vary.
// If we are debugging a single pixel, we do not want to loop across all the
// pixels to be sampled. If we are using openmp then distribute pixels across
// threads. If we are not using openmp, and still want to consider all pixels,
// then use a normal for loop.
#ifdef _DEBUG_SINGLEPIX
    ipix_big = par -> dbg_ipix;
#else //_DEBUG_SINGLEPIX
  #ifdef _WITH_OMP
    #pragma omp for
  #endif //_WITH_OMP
    for(ipix_big = par -> ipix_0; ipix_big < par -> ipix_f; ipix_big++)
#endif //_DEBUG_SINGLEPIX
    {
	     printf("Node %d, thread %d, pixel %d\n", NodeThis, ithr, ipix_big);
       if (par -> flag_use_marginal)
       {
         clean_pixel_from_marginal(par, rng, pst_old[ithr], pst_new[ithr], ipix_big);
       }
       else
       {
         clean_pixel(par, rng, pst_old[ithr], ipix_big);
       }
     } //end omp for
     end_rng(rng);
  }//end omp parallel

  // Now do the cleanup.
  for (ii = 0; ii < n_threads; ii++)
  {
    pixel_state_free(pst_old[ii], par);
  }
  free(pst_old);

  if (par -> flag_use_marginal)
  {
    for (ii = 0; ii < n_threads; ii++)
    {
      pixel_state_free(pst_new[ii], par);
    }
    free(pst_new);
  }
  if(NodeThis == 0)
  {
    printf("Writing output\n");
  }

  write_output(par);
  param_bfore_free(par);

#ifdef _WITH_MPI
  MPI_Finalize();
#endif //_WITH_MPI

  return 0;
}
