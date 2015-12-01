#include "common.h"

int main(int argc,char **argv)
{ 
  char fname_init[256];
  if(argc!=2) {
    printf("Usage: fg_rm.x param_file\n");
    exit(0);
  }
  sprintf(fname_init,"%s",argv[1]);

  ParamFGRM *par=read_params(fname_init); //Done - 2 CHECKs
  //  clean_maps(par);
  //  write_maps(par);
  param_fgrm_free(par);

  return 0;
}
