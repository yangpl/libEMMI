/* EM modeling and inversion using FDTD method
 *------------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------------*/
#include <mpi.h>

#include "cstd.h"
#include "mpi_info.h"
#include "emf.h"
#include "acqui.h"


int iproc, nproc, ierr;

void emf_init(emf_t *emf);
void emf_close(emf_t *emf);

void acqui_init(acqui_t *acqui, emf_t * emf);
void acqui_close(acqui_t *acqui);

void do_modeling(acqui_t *acqui, emf_t *emf);
void do_fwi(acqui_t *acqui, emf_t *emf);


int main(int argc, char* argv[])
{
  emf_t *emf;
  acqui_t *acqui;
  
  char cur_time[128];
  time_t      t;
  struct tm*  ptm;

  MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  initargs(argc,argv);
  
  emf = (emf_t *)malloc(sizeof(emf_t));
  acqui = (acqui_t *)malloc(sizeof(acqui_t));


  if(!getparint("verb", &emf->verb)) emf->verb = (iproc==0)?1:0; 
  /* emf->verbose=1 only on master process */
  if(!getparint("mode", &emf->mode)) emf->mode = 0;/* 0=modeling; 1=inversion*/
  if(emf->verb){
    t = time(NULL);
    ptm = localtime(&t);
    strftime(cur_time, 128, "%d-%b-%Y %H:%M:%S", ptm);
    printf("  Current date and time: %s\n", cur_time);

    printf("=====================================================\n");
    printf("    Welcome to EM modeling and inversion code        \n");
    printf("            Author: Pengliang Yang                   \n");
    printf("            E-mail: ypl.2100@gmail.com               \n");
    printf("  Copyright (c) 2020. All rights reserved. \n");
    printf("=====================================================\n");
    printf("\n");
    if(emf->mode==0)
      printf("Task: Electromagnetic forward modeling \n");
    else if(emf->mode==1)
      printf("Task: Inversion for medium parameters \n");
    else if(emf->mode==3)
      printf("Task: Dot product test for finite difference and adjoint gradient \n");
    else
      printf("Task: Output inversion gradient \n");
    printf("\n");
    printf("=====================================================\n");
  }
  emf_init(emf);
  acqui_init(acqui, emf);
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  
  if(emf->mode==0)  do_modeling(acqui, emf); /* emf->mode=0 */
  else do_fwi(acqui, emf);


  acqui_close(acqui);
  emf_close(emf);
  free(emf);
  free(acqui);
    
  MPI_Finalize();
  return 0;
}
