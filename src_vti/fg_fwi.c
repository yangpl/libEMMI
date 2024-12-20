/* function and gradient evaluation for EM full waveform inversion
 *--------------------------------------------------------------------
 *
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 *--------------------------------------------------------------------*/
#include <mpi.h>

#include "cstd.h"
#include "acqui.h"
#include "emf.h"
#include "fwi.h"
#include "interp.h"
#include "constants.h"
#include "mpi_info.h"

acqui_t *acqui;
emf_t *emf;
fwi_t *fwi;

/*---------------- modelling part ------------*/
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
void inject_electric_src_adj(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg, int it);

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


/*----------------- inversion part -------------------*/
void write_data(acqui_t *acqui, emf_t *emf, char *fname, float _Complex ***dcal_fd);
void read_data(acqui_t *acqui, emf_t *emf);
void write_misfit(acqui_t *acqui, emf_t *emf, char *fname, float ***obj);
void cal_uncertainty_noise(acqui_t *acqui, emf_t *emf, fwi_t *fwi);

void homogenization(emf_t *emf);/* homogenization method by Davedycheva 2003 */
void build_gradient(emf_t *emf, fwi_t *fwi);
void log_exp_scaling(int n, float *x, float *y, int flag);

void fcost_adjsrc_init();
float fcost_adjoint_source(acqui_t *acqui, emf_t *emf);


void fg_mod_reg(acqui_t *acqui, emf_t *emf, fwi_t *fwi, float *x, float *g);
void Hv_mod_reg(emf_t *emf, fwi_t *fwi, float *r, float *Hv);

void computing_box_init(acqui_t *acqui, emf_t *emf, int adj);
void computing_box_close(emf_t *emf, int adj);

#ifdef GPU
void cuda_modeling(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg, int adj);
#endif

/*--------------------------------------------------------*/
void fg_fwi_init(acqui_t *acqui_, emf_t *emf_, fwi_t *fwi_)
{
  char *fibathy, *fmask;
  int irec, ntrace, k, i1, i2, istat;
  FILE *fp;
  
  acqui = acqui_;
  emf = emf_;
  fwi = fwi_;
  
  if(!getparint("addnoise", &emf->addnoise)) emf->addnoise = 1;//1=addnoise;0=not
  if(!getparint("invmask", &emf->invmask)) emf->invmask = 0;
  if(!getparfloat("amp_perc", &emf->amp_perc)) emf->amp_perc = 0.03;
  /* uncertainty, default=1% noise/deviation from the true value*/
  if(!getparfloat("delta_phi", &emf->delta_phi)) emf->delta_phi = 1.5*PI/180.;
  if(!getparfloat("noisefloorE", &emf->noisefloorE)) emf->noisefloorE = 1e-16;
  if(!getparfloat("noisefloorH", &emf->noisefloorH)) emf->noisefloorH = 1e-14;
  if(!getparfloat("offset_start", &emf->offset_start)) emf->offset_start = 1000.;
  if(!getparfloat("gamma1", &fwi->gamma1)) fwi->gamma1 = 1.e3;//Tikhonov
  if(!getparfloat("gamma2", &fwi->gamma2)) fwi->gamma2 = 0.;//TV
  if(!getparint("preco", &fwi->preco)) fwi->preco = 1;/* 1 = precondition on, 0 = off */
  if(!getparint("r1", &fwi->r1)) fwi->r1 = 1;//no triangle smoothing if r1=1
  if(!getparint("r2", &fwi->r2)) fwi->r2 = 1;//no triangle smoothing if r2=1
  if(!getparint("r3", &fwi->r3)) fwi->r3 = 1;//no triangle smoothing if r3=1
  if(!getparint("repeat", &fwi->repeat)) fwi->repeat = 1;//repeat times of smoothing
  if(!getparint("nmutebathy", &emf->nmutebathy)) emf->nmutebathy = 5;//mute 5 points below bathy
  if(emf->verb) {
    printf("preco=%d\n", fwi->preco);
    printf("offset >=%g \n", emf->offset_start);
    printf("Amplitude uncertainty: %g\%%\n", emf->amp_perc*100.);
    printf("noisefloorE=%g noisefloorH=%g\n", emf->noisefloorE, emf->noisefloorH);
    printf("Tikhonov regularization: gamma1=%g\n", fwi->gamma1);
    printf("TV regularization:       gamma2=%g\n", fwi->gamma2);
  }


  emf->ibathy = alloc2int(emf->n1, emf->n2);
  emf->offset_ok = alloc1int(acqui->nrec);
  emf->mask = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->dobs_fd = alloc3complexf(acqui->nrec, emf->nfreq, emf->nchrec);
  emf->dcal_fd = alloc3complexf(acqui->nrec, emf->nfreq, emf->nchrec);
  emf->dres_fd = alloc3complexf(acqui->nrec, emf->nfreq, emf->nchrec);
  emf->delta_emf = alloc3float(acqui->nrec, emf->nfreq, emf->nchrec);
  emf->obj = alloc3float(acqui->nrec, emf->nfreq, emf->nchrec);//significant misfit
  fwi->xref = alloc1float(fwi->n);
  fwi->g_h = alloc3float(emf->n1, emf->n2, emf->n3);//grad_h
  fwi->g_v = alloc3float(emf->n1, emf->n2, emf->n3);//grad_v
  fwi->fcost_list = alloc1float(nproc);
  fwi->ndp_list = alloc1int(nproc);

  if(!(getparstring("fibathy", &fibathy))){
    if(!getparfloat("waterdepth", &emf->waterdepth)) err("Must have waterdepth= ");
    for(i2=0; i2<emf->n2; i2++)
      for(i1=0; i1<emf->n1; i1++)
	emf->ibathy[i2][i1] = ceil(emf->waterdepth/emf->d3+0.5);//the index of bathymetry along z
  }else{
    fp = fopen(fibathy, "rb");
    if(fp==NULL) err("canot open file");
    istat = fread(&emf->ibathy[0][0], sizeof(float), emf->n1*emf->n2, fp);
    if(istat != emf->n1*emf->n2) err("size parameter does not match the file - fibathy!");
    fclose(fp);
    for(i2=0; i2<emf->n2; i2++)
      for(i1=0; i1<emf->n1; i1++)
	emf->ibathy[i2][i1] += emf->nmutebathy;//the index of bathymetry along z    
    if(emf->verb) printf("bathymetry specified!\n");
  }

  if(emf->invmask){
    if(!(getparstring("fmask", &fmask))) err("Must have fmask= ");
    else{
      fp = fopen(fmask, "rb");
      if(fp==NULL) err("canot open file");
      istat = fread(&emf->mask[0][0][0], sizeof(int), emf->n123, fp);
      if(istat != emf->n123) err("size parameter does not match the file - fmask!");
      fclose(fp);
    }
    if(emf->verb) printf("inversion mask set up!\n");
  }
  
  fwi->firstgrad = 1;
  fwi->alpha = 1;
  fwi->iter = 0;

  ntrace = 0;
  for(irec = 0; irec<acqui->nrec; irec++){
    float d1 = acqui->src_x1[0]-acqui->rec_x1[irec];
    float d2 = acqui->src_x2[0]-acqui->rec_x2[irec];
    float tmp = sqrtf(d1*d1 + d2*d2);//offset between source and receiver
    if(tmp>emf->offset_start) {
      emf->offset_ok[irec] = 1;
      ntrace++;
    }else emf->offset_ok[irec] = 0;
  }
  if(emf->verb) printf("%d of %d receivers in offset of interest\n", ntrace, acqui->nrec);
  ntrace *= emf->nchsrc*emf->nchrec*emf->nfreq;
  ierr = MPI_Gather(&ntrace, 1, MPI_INT, fwi->ndp_list, 1, MPI_INT, 0, MPI_COMM_WORLD);
  ierr = MPI_Allreduce(&ntrace, &fwi->ndp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);/*sum over process*/
  if(emf->verb){
    printf("total number of data points: ndp=%d\n", fwi->ndp);
    for(k=0; k<nproc; k++) printf("isrc=%d ndp=%d\n", acqui->shot_idx[k], fwi->ndp_list[k]);
  }
  
  read_data(acqui, emf);/* initialize dobs_fd */
  if(emf->verb) printf("read observed data completed!\n");
  cal_uncertainty_noise(acqui, emf, fwi);
  if(emf->verb) printf("data uncertainty computed!\n");

  fcost_adjsrc_init();
  memcpy(emf->rho_h[0][0], emf->rho11[0][0], emf->n123*sizeof(float));//set it to be rho11
  memcpy(emf->rho_v[0][0], emf->rho33[0][0], emf->n123*sizeof(float));//set it to be rho33
  if(emf->verb) printf("----------- fwi init completed ------------\n");
}

void fg_fwi_close()
{
  free2int(emf->ibathy);
  free1int(emf->offset_ok);
  free3float(emf->mask);
  free3complexf(emf->dobs_fd);
  free3complexf(emf->dcal_fd);
  free3complexf(emf->dres_fd);
  free3float(emf->delta_emf);
  free3float(emf->obj);
  free1float(fwi->xref);
  free3float(fwi->g_h);
  free3float(fwi->g_v);
  free1float(fwi->fcost_list);
  free1int(fwi->ndp_list);
}


/*----------------------------------------------------------*/
float fg_fwi(float *x, float *g)
/*< misfit function and gradient evaluation of FWI >*/
{
  interp_t *interp_rg;/* interpolation for regular grid */
  interp_t *interp_sg;/* interpolation for staggered grid */
  int adj, it, k, offset;
  FILE *fp;
  float tmp, s1, s2;

  fwi->fcost_dat = 0.;
  fwi->fcost_mod = 0.;
  emf->gnopt = 0;
  memset(fwi->g_h[0][0], 0, emf->n123*sizeof(float));
  memset(fwi->g_v[0][0], 0, emf->n123*sizeof(float));
  offset = 0;
  for(k=0; k<fwi->npar; k++){
    if(fwi->idxpar[k]==1){//Rv with index 1
      log_exp_scaling(emf->n123, x+offset, &emf->rho_v[0][0][0], 2);
      offset += emf->n123;
    }else if(fwi->idxpar[k]==2){// Rh with index 2
      log_exp_scaling(emf->n123, x+offset, &emf->rho_h[0][0][0], 2);
      offset += emf->n123;
    }
  }

  /*-----------------------------------------------------------------------*/
  /* step 0: convert inversion parameter x into physical resistiviity model */
  /*-----------------------------------------------------------------------*/
  homogenization(emf);//get rho11, rho22 and rho33 from isotropic rho
  sanity_check(emf);  
  
  emf->stf = alloc1float(emf->nt);
  memset(emf->stf, 0, emf->nt*sizeof(float));
  emf->stf[0] = 1.;  
  emf->dres_td = alloc3float(acqui->nrec, emf->nchrec, emf->nt);
  memset(&emf->dcal_fd[0][0][0], 0, emf->nfreq*acqui->nrec*emf->nchrec*sizeof(float _Complex));
  memset(&emf->dres_fd[0][0][0], 0, emf->nfreq*acqui->nrec*emf->nchrec*sizeof(float _Complex));

  computing_box_init(acqui, emf, 0);
  computing_box_init(acqui, emf, 1);
  emf->expfactor = alloc2complexf(emf->nfreq, emf->nt);
  for(int ifreq=0; ifreq<emf->nfreq; ifreq++){
    float _Complex omegap = (1.0+I)*sqrt(emf->omega0*emf->omegas[ifreq]);
    for(it=0; it<emf->nt; it++) emf->expfactor[it][ifreq] = cexp(I*omegap*(it+0.5)*emf->dt);
  }

  interp_rg = (interp_t *)malloc(sizeof(interp_t));
  interp_sg = (interp_t *)malloc(sizeof(interp_t));
  interpolation_init(acqui, emf, interp_rg, interp_sg);
  interpolation_weights(acqui, emf, interp_rg, interp_sg);
  extend_model_init(emf);
  fdtd_init(emf);
  if(emf->airwave) airwave_bc_init(emf);
  dtft_emf_init(emf, 0);
  dtft_emf_init(emf, 1);
  
  /*=================================================================*/
  if(emf->verb) printf("------stage 1: forward modelling ---------\n");
  adj = 0;
#ifdef GPU
  cuda_modeling(acqui, emf, interp_rg, interp_sg, adj); /* mode=0 */
#else
  fdtd_null(emf);
  for(it=0; it<emf->nt; it++){
    if(it%50==0 && emf->verb) printf("it-----%d\n", it);

    fdtd_curlH(emf, it, adj);
    inject_electric_src_fwd(acqui, emf, interp_rg, interp_sg, it);
    fdtd_update_E(emf, it, adj);
    if(emf->airwave) airwave_bc_update_E(emf, emf->E1, emf->E2);
    
    dtft_emf(emf, it, adj);
    fdtd_curlE(emf, it, adj); 
    fdtd_update_H(emf, it, adj); 
    if(emf->airwave) airwave_bc_update_H(emf, emf->H3, emf->H1, emf->H2);
    
    if(it%100==0){/* convergence check */
      emf->ncorner = check_convergence(emf, adj);
      if(iproc==0) printf("%d corners of the cube converged!\n", emf->ncorner);
      if(emf->ncorner==8) { break; }/* all 8 corners converged, exit now */
    }    
  }
#endif
  compute_green_function(emf, adj);
  extract_emf(acqui, emf, interp_rg, interp_sg);

  char fname[sizeof("syn_0000.txt")];
  sprintf(fname, "syn_%04d.txt", acqui->shot_idx[iproc]);
  write_data(acqui, emf, fname, emf->dcal_fd);
  printf("modelling for syn_%04d.txt completed!\n", acqui->shot_idx[iproc]);

  /*=================================================================*/
  if(emf->verb) printf("------stage 2: adjoint source estimation ---------\n");
  tmp = fcost_adjoint_source(acqui, emf);/* data misfit on each processor */
  ierr = MPI_Allreduce(&tmp, &fwi->fcost_dat, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);/*sum over all process*/
  ierr = MPI_Gather(&tmp, 1, MPI_FLOAT, fwi->fcost_list, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  if(emf->verb){
    printf("fcost_dat=%g, rmse=%g\n", fwi->fcost_dat, sqrt(2.*fwi->fcost_dat/fwi->ndp));
    for(k=0; k<nproc; k++) {
      tmp = 2.*fwi->fcost_list[k]/fwi->ndp_list[k];
      printf("isrc=%d rmse=%g\n", acqui->shot_idx[k], sqrt(tmp));
    }
  }

  /*=================================================================*/
  if(emf->verb) printf("-----------stage 3: adjoint modelling ----------\n");
  adj = 1;
#ifdef GPU
  cuda_modeling(acqui, emf, interp_rg, interp_sg, adj); /* mode=0 */
#else
  memset(emf->stf, 0, emf->nt*sizeof(float));
  emf->stf[0] = 1.;
  fdtd_null(emf);
  for(it=0; it<emf->nt; it++){
    if(it%50==0 && emf->verb) printf("it-----%d\n", it);

    fdtd_curlH(emf, it, adj);
    inject_electric_src_adj(acqui, emf, interp_rg, interp_sg, it);
    fdtd_update_E(emf, it, adj);
    if(emf->airwave) airwave_bc_update_E(emf, emf->E1, emf->E2);
    
    dtft_emf(emf, it, adj);
    fdtd_curlE(emf, it, adj); 
    fdtd_update_H(emf, it, adj); 
    if(emf->airwave) airwave_bc_update_H(emf, emf->H3, emf->H1, emf->H2);

    if(it%100==0){/* convergence check */
      emf->ncorner = check_convergence(emf, adj);
      if(iproc==0) printf("%d corners of the cube converged!\n", emf->ncorner);
      if(emf->ncorner==8) { break; }/* all 8 corners converged, exit now */
    }    
  }
#endif
  compute_green_function(emf, adj);

  /*=================================================================*/
  if(emf->verb) printf("-------- stage 4: build the gradient ----------\n");
  build_gradient(emf, fwi);//dJ/d(log(rho))
  offset = 0;
  for(k=0; k<fwi->npar; k++){
    if(fwi->idxpar[k]==1){//Rv with index 1
      ierr = MPI_Allreduce(&fwi->g_v[0][0][0], g+offset, emf->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
      memcpy(&fwi->g_v[0][0][0], g+offset, emf->n123*sizeof(float));
      offset += emf->n123;
    }else if(fwi->idxpar[k]==2){// Rh with index 2
      ierr = MPI_Allreduce(&fwi->g_h[0][0][0], g+offset, emf->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
      memcpy(&fwi->g_h[0][0][0], g+offset, emf->n123*sizeof(float));
      offset += emf->n123;
    }
  }
  if(emf->verb){
    fp = fopen("gradient_fwi_logRv", "wb");
    fwrite(&fwi->g_v[0][0][0], emf->n123*sizeof(float), 1, fp);
    fclose(fp);

    fp = fopen("gradient_fwi_logRh", "wb");
    fwrite(&fwi->g_h[0][0][0], emf->n123*sizeof(float), 1, fp);
    fclose(fp);
  }
  ierr = MPI_Barrier(MPI_COMM_WORLD);

  /*=================================================================*/
  extend_model_close(emf);
  fdtd_close(emf);
  if(emf->airwave) airwave_bc_close(emf);
  dtft_emf_close(emf, 0);
  dtft_emf_close(emf, 1);
  interpolation_close(emf, interp_rg, interp_sg);
  free(interp_rg);
  free(interp_sg);
  free(emf->stf);
  free3float(emf->dres_td);

  free2complexf(emf->expfactor);
  computing_box_close(emf, 0);
  computing_box_close(emf, 1);

  
  /*===================================================================*/
  if(emf->verb) printf("-------- stage 5: model regularization ----------\n");
  fg_mod_reg(acqui, emf, fwi, x, g);
  fwi->fcost = fwi->fcost_dat + fwi->fcost_mod;
  if(emf->verb){
    printf("fcost_dat=%g\n", fwi->fcost_dat);
    printf("fcost_mod=%g\n", fwi->fcost_mod);
    printf("fcost=%g\n", fwi->fcost);
  }

  /* l-BFGS algorithm is not scale invariant, we need to estimate a scaling factor
     alpha at 1st iteration and use it for FWI objective function in later iterations */
  if(emf->mode == 1 && fwi->firstgrad){
    s1 = 0;
    s2 = 0;
    for(k = 0; k<fwi->n; k++){
      s1 += fabs(x[k]);
      s2 += fabs(g[k]);
    }
    fwi->alpha = 2e-3*s1/s2;
    if(emf->verb) printf("|x|_1 = %g, |g|_1 = %g scaling = %g\n", s1, s2, fwi->alpha);
  }
  for(k=0; k<fwi->n; k++) g[k] *= fwi->alpha;
  if(fwi->firstgrad) fwi->firstgrad = 0;

  return fwi->fcost*fwi->alpha; 
}


/*----------------------------------------------------------*/
void Hv_fwi(float *x, float *v, float *Hv)
/*< compute Gauss-Newton Hessian vector product for FWI >*/
{

}
