/* EM modeling using FDTD method
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include <mpi.h>

#include "cstd.h"

#include "acqui.h"
#include "emf.h"
#include "interp.h"
#include "mpi_info.h"
#include "constants.h"


void emf_init(emf_t *emf);
void emf_close(emf_t *emf);

void extend_model_init(emf_t *emf);
void extend_model_close(emf_t *emf);
void sanity_check(emf_t *emf);

void acqui_init(acqui_t *acqui, emf_t * emf);
void acqui_close(acqui_t *acqui);

void fdtd_init(emf_t *emf);
void fdtd_null(emf_t *emf);
void fdtd_close(emf_t *emf);
void fdtd_curlH(emf_t *emf, int it, int adj);
void fdtd_update_E(emf_t *emf, int it, int adj);
void fdtd_curlE(emf_t *emf, int it, int adj);
void fdtd_update_H(emf_t *emf, int it, int adj);

void interpolation_init(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg);
void interpolation_close(emf_t *emf, interp_t *interp_rg, interp_t *interp_sg);
void interpolation_weights(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg);

void inject_electric_src_fwd(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg, int it);

void airwave_bc_init(emf_t *emf);
void airwave_bc_close(emf_t *emf);
void airwave_bc_update_H(emf_t *emf, float ***H3, float ***H1, float ***H2);
void airwave_bc_update_E(emf_t *emf, float ***E1, float ***E2);

void dtft_emf_init(emf_t *emf, int adj);
void dtft_emf_close(emf_t *emf, int adj);
void dtft_emf(emf_t *emf, int it, int adj);

void compute_green_function(emf_t *emf, int adj);
void extract_emf(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg);

void write_data(acqui_t *acqui, emf_t *emf, char *fname, float _Complex ***dcal_fd);

int check_convergence(emf_t *emf, int adj);


void computing_box_init(acqui_t *acqui, emf_t *emf, int adj);
void computing_box_close(emf_t *emf, int adj);

#ifdef GPU
void cuda_modeling(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg, int adj);
#endif

/*----------------------------------------------------*/
void do_modeling(acqui_t *acqui, emf_t *emf)
{
  int it, adj;
  interp_t *interp_rg, *interp_sg;

  interp_rg = (interp_t *)malloc(sizeof(interp_t));
  interp_sg = (interp_t *)malloc(sizeof(interp_t));
  interpolation_init(acqui, emf, interp_rg, interp_sg);
  interpolation_weights(acqui, emf, interp_rg, interp_sg);

  sanity_check(emf);

  emf->stf = alloc1float(emf->nt);
  memset(emf->stf, 0, emf->nt*sizeof(float));
  emf->stf[0] = 1.;
  emf->dcal_fd = alloc3complexf(acqui->nrec, emf->nfreq, emf->nchrec);
  memset(&emf->dcal_fd[0][0][0], 0, acqui->nrec*emf->nchrec*emf->nfreq*sizeof(float _Complex));

  computing_box_init(acqui, emf, 0);
  emf->expfactor = alloc2complexf(emf->nfreq, emf->nt);
  for(int ifreq=0; ifreq<emf->nfreq; ifreq++){
    float _Complex omegap = (1.0+I)*sqrt(emf->omega0*emf->omegas[ifreq]);
    for(it=0; it<emf->nt; it++) emf->expfactor[it][ifreq] = cexp(I*omegap*(it+0.5)*emf->dt);
  }


  adj = 0;
  extend_model_init(emf);
  fdtd_init(emf);
  fdtd_null(emf);
  if(emf->airwave) airwave_bc_init(emf);
  dtft_emf_init(emf, adj);


#ifdef GPU    
  cuda_modeling(acqui, emf, interp_rg, interp_sg, adj); /* mode=0 */
#else
  for(it=0; it<emf->nt; it++){
    if(it%50==0 && emf->verb) printf("it-----%d\n", it);

    fdtd_curlH(emf, it, adj);
    inject_electric_src_fwd(acqui, emf, interp_rg, interp_sg, it);
    fdtd_update_E(emf, it, adj);
    if(emf->airwave) airwave_bc_update_E(emf, emf->E1, emf->E2);
    
    fdtd_curlE(emf, it, adj); 
    fdtd_update_H(emf, it, adj); 
    if(emf->airwave) airwave_bc_update_H(emf, emf->H3, emf->H1, emf->H2);
    
    dtft_emf(emf, it, adj);
    
    if(it%100==0){/* convergence check */
      emf->ncorner = check_convergence(emf, adj);
      if(iproc==0) printf("%d corners of the cube converged!\n", emf->ncorner);
      if(emf->ncorner==8) break;/* all 8 corners converged, exit now */
    }    
  }
#endif
  compute_green_function(emf, adj);
  extract_emf(acqui, emf, interp_rg, interp_sg);

  char fname[sizeof("emf_0000.txt")];
  sprintf(fname, "emf_%04d.txt", acqui->shot_idx[iproc]);
  write_data(acqui, emf, fname, emf->dcal_fd);
  printf("modelling for emf_%04d.txt completed!\n", acqui->shot_idx[iproc]);  

  extend_model_close(emf);
  fdtd_close(emf);
  if(emf->airwave) airwave_bc_close(emf);
  dtft_emf_close(emf, 0);
  interpolation_close(emf, interp_rg, interp_sg);
  free(interp_rg);
  free(interp_sg);

  free(emf->stf);

  free2complexf(emf->expfactor);
  computing_box_close(emf, 0);

}
