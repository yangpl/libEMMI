/* CSEM FDTD modeling using Nvidia GPU
 *
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 */
#include <mpi.h>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>



#ifdef __cplusplus
extern "C" {
#endif
  
#include "cstd.h"
#include "acqui.h"
#include "emf.h"
#include "interp.h"
#include "constants.h"
#include "mpi_info.h"
  
void compute_green_function(emf_t *emf);
void extract_emf(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg);
void write_data(acqui_t *acqui, emf_t *emf, char *fname, float _Complex ***dcal_fd);
  
#ifdef __cplusplus
}
#endif

#include "cuda_fdtd.cuh"


float c11, c21, c12, c22, c13, c23;

cudaError_t status;
dim3 dimBlock, dimGrid, dimGrid_dft;

float *d_inveps11, *d_inveps22, *d_inveps33;
float *d_E1, *d_E2, *d_E3, *d_H1, *d_H2, *d_H3;
float *d_curlE1, *d_curlE2, *d_curlE3, *d_curlH1, *d_curlH2, *d_curlH3;
float *d_memD2E1, *d_memD3E1, *d_memD1E2, *d_memD3E2, *d_memD1E3, *d_memD2E3;
float *d_memD2H1, *d_memD3H1, *d_memD1H2, *d_memD3H2, *d_memD1H3, *d_memD2H3;
float *d_a1, *d_b1, *d_a2, *d_b2, *d_a3, *d_b3;
float *d_omegas;

cuFloatComplex *d_fwd_E1, *d_fwd_E2, *d_fwd_E3;
cuFloatComplex *d_backup, *d_expfactor;
int *d_corner_id, *h_ncorner;

cufftHandle fftPlan;
cuFloatComplex *d_sH1kxky, *d_sH2kxky, *d_emfft, *d_emfft0;
float *d_sE12kxky;


int *d_rg_src_i1, *d_rg_src_i2, *d_rg_src_i3;
float *d_rg_src_w1, *d_rg_src_w2, *d_rg_src_w3;
int *d_sg_src_i1, *d_sg_src_i2, *d_sg_src_i3;
float *d_sg_src_w1, *d_sg_src_w2, *d_sg_src_w3;

int *d_rg_rec_i1, *d_rg_rec_i2, *d_rg_rec_i3;
float *d_rg_rec_w1, *d_rg_rec_w2, *d_rg_rec_w3;
int *d_sg_rec_i1, *d_sg_rec_i2, *d_sg_rec_i3;
float *d_sg_rec_w1, *d_sg_rec_w2, *d_sg_rec_w3;

int *d_chrec, *d_chsrc;
float *d_dres_td;//adjoint source


void cuda_fdtd_init(emf_t *emf)
{
  int ic, *h_chrec, *h_chsrc;
  int i1, i2, i3;

  int nchsrc = emf->nchsrc;
  int nchrec = emf->nchrec;
  int corner_id[8];

  //4-th order staggered FD, backward difference using shared memory:
  //c1*(D[0]-D[-1])+c2*(D[1]-D[-2])
  c11 = fd_c1/emf->d1;
  c21 = fd_c2/emf->d1;
  c12 = fd_c1/emf->d2;
  c22 = fd_c2/emf->d2;
  c13 = fd_c1/emf->d3;
  c23 = fd_c2/emf->d3;

  /* allocate memory on device */
  cudaMalloc(&d_inveps11, emf->n123pad*sizeof(float));
  cudaMalloc(&d_inveps22, emf->n123pad*sizeof(float));
  cudaMalloc(&d_inveps33, emf->n123pad*sizeof(float));
  cudaMalloc(&d_E1, emf->n123pad*sizeof(float));
  cudaMalloc(&d_E2, emf->n123pad*sizeof(float));
  cudaMalloc(&d_E3, emf->n123pad*sizeof(float));
  cudaMalloc(&d_H1, emf->n123pad*sizeof(float));
  cudaMalloc(&d_H2, emf->n123pad*sizeof(float));
  cudaMalloc(&d_H3, emf->n123pad*sizeof(float));
  cudaMalloc(&d_curlE1, emf->n123pad*sizeof(float));
  cudaMalloc(&d_curlE2, emf->n123pad*sizeof(float));
  cudaMalloc(&d_curlE3, emf->n123pad*sizeof(float));
  cudaMalloc(&d_curlH1, emf->n123pad*sizeof(float));
  cudaMalloc(&d_curlH2, emf->n123pad*sizeof(float));
  cudaMalloc(&d_curlH3, emf->n123pad*sizeof(float));
  cudaMalloc(&d_memD2E1, emf->n1pad*2*emf->nb*emf->n3pad*sizeof(float));
  cudaMalloc(&d_memD3E1, emf->n1pad*emf->n2pad*2*emf->nb*sizeof(float));
  cudaMalloc(&d_memD1E2, 2*emf->nb*emf->n2pad*emf->n3pad*sizeof(float));
  cudaMalloc(&d_memD3E2, emf->n1pad*emf->n2pad*2*emf->nb*sizeof(float));
  cudaMalloc(&d_memD1E3, 2*emf->nb*emf->n2pad*emf->n3pad*sizeof(float));
  cudaMalloc(&d_memD2E3, emf->n1pad*2*emf->nb*emf->n3pad*sizeof(float));
  cudaMalloc(&d_memD2H1, emf->n1pad*2*emf->nb*emf->n3pad*sizeof(float));
  cudaMalloc(&d_memD3H1, emf->n1pad*emf->n2pad*2*emf->nb*sizeof(float));
  cudaMalloc(&d_memD1H2, 2*emf->nb*emf->n2pad*emf->n3pad*sizeof(float));
  cudaMalloc(&d_memD3H2, emf->n1pad*emf->n2pad*2*emf->nb*sizeof(float));
  cudaMalloc(&d_memD1H3, 2*emf->nb*emf->n2pad*emf->n3pad*sizeof(float));
  cudaMalloc(&d_memD2H3, emf->n1pad*2*emf->nb*emf->n3pad*sizeof(float));
  cudaMalloc(&d_fwd_E1, emf->n123pad*emf->nfreq*sizeof(cuFloatComplex));
  cudaMalloc(&d_fwd_E2, emf->n123pad*emf->nfreq*sizeof(cuFloatComplex));
  cudaMalloc(&d_fwd_E3, emf->n123pad*emf->nfreq*sizeof(cuFloatComplex));
  cudaMalloc(&d_a1, emf->nb*sizeof(float));
  cudaMalloc(&d_b1, emf->nb*sizeof(float));
  cudaMalloc(&d_a2, emf->nb*sizeof(float));
  cudaMalloc(&d_b2, emf->nb*sizeof(float));
  cudaMalloc(&d_a3, emf->nb*sizeof(float));
  cudaMalloc(&d_b3, emf->nb*sizeof(float));

  cudaMalloc(&d_corner_id, 8*sizeof(int));
  cudaHostAlloc(&h_ncorner, sizeof(int), cudaHostAllocMapped);	
  cudaMalloc(&d_backup, 8*sizeof(cuFloatComplex));
  cudaMalloc(&d_expfactor, emf->nfreq*emf->nt*sizeof(cuFloatComplex));
  status = cudaGetLastError();
  if (cudaSuccess!=status) { printf("Failed to allocate memory on device - fdtd !\n"); exit(0); }
  
  //initialize memory on device
  cudaMemcpy(d_inveps11, emf->inveps11[0][0], emf->n123pad*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_inveps22, emf->inveps22[0][0], emf->n123pad*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_inveps33, emf->inveps33[0][0], emf->n123pad*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemset(d_E1, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_E2, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_E3, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_H1, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_H2, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_H3, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_curlE1, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_curlE2, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_curlE3, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_curlH1, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_curlH2, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_curlH3, 0, emf->n123pad*sizeof(float));
  cudaMemset(d_memD2E1, 0, 2*emf->nb*emf->n1pad*emf->n3pad*sizeof(float));
  cudaMemset(d_memD3E1, 0, 2*emf->nb*emf->n1pad*emf->n2pad*sizeof(float));
  cudaMemset(d_memD1E2, 0, 2*emf->nb*emf->n2pad*emf->n3pad*sizeof(float));
  cudaMemset(d_memD3E2, 0, 2*emf->nb*emf->n1pad*emf->n2pad*sizeof(float));
  cudaMemset(d_memD1E3, 0, 2*emf->nb*emf->n2pad*emf->n3pad*sizeof(float));
  cudaMemset(d_memD2E3, 0, 2*emf->nb*emf->n1pad*emf->n3pad*sizeof(float));
  cudaMemset(d_memD2H1, 0, 2*emf->nb*emf->n1pad*emf->n3pad*sizeof(float));
  cudaMemset(d_memD3H1, 0, 2*emf->nb*emf->n1pad*emf->n2pad*sizeof(float));
  cudaMemset(d_memD1H2, 0, 2*emf->nb*emf->n2pad*emf->n3pad*sizeof(float));
  cudaMemset(d_memD3H2, 0, 2*emf->nb*emf->n1pad*emf->n2pad*sizeof(float));
  cudaMemset(d_memD1H3, 0, 2*emf->nb*emf->n2pad*emf->n3pad*sizeof(float));
  cudaMemset(d_memD2H3, 0, 2*emf->nb*emf->n1pad*emf->n3pad*sizeof(float));
  cudaMemset(d_fwd_E1, 0, emf->n123pad*emf->nfreq*sizeof(cuFloatComplex));
  cudaMemset(d_fwd_E2, 0, emf->n123pad*emf->nfreq*sizeof(cuFloatComplex));
  cudaMemset(d_fwd_E3, 0, emf->n123pad*emf->nfreq*sizeof(cuFloatComplex));
  cudaMemcpy(d_a1, emf->a1, emf->nb*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b1, emf->b1, emf->nb*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_a2, emf->a2, emf->nb*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b2, emf->b2, emf->nb*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_a3, emf->a3, emf->nb*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b3, emf->b3, emf->nb*sizeof(float), cudaMemcpyHostToDevice);

  i1 = emf->nbe;
  i2 = emf->nbe;
  i3 = emf->nbe;
  corner_id[0] = i1 + emf->n1pad*(i2 + emf->n2pad*i3);
  i1 = emf->nbe + emf->n1 - 1;
  i2 = emf->nbe;
  i3 = emf->nbe;
  corner_id[1] = i1 + emf->n1pad*(i2 + emf->n2pad*i3);
  i1 = emf->nbe;
  i2 = emf->nbe + emf->n2 - 1;
  i3 = emf->nbe;
  corner_id[2] = i1 + emf->n1pad*(i2 + emf->n2pad*i3);
  i1 = emf->nbe;
  i2 = emf->nbe;
  i3 = emf->nbe + emf->n3 - 1;
  corner_id[3] = i1 + emf->n1pad*(i2 + emf->n2pad*i3);
  i1 = emf->nbe + emf->n1 -1;
  i2 = emf->nbe + emf->n2 -1;
  i3 = emf->nbe;
  corner_id[4] = i1 + emf->n1pad*(i2 + emf->n2pad*i3);
  i1 = emf->nbe + emf->n1 - 1;
  i2 = emf->nbe;
  i3 = emf->nbe + emf->n3 - 1;
  corner_id[5] = i1 + emf->n1pad*(i2 + emf->n2pad*i3);
  i1 = emf->nbe;
  i2 = emf->nbe + emf->n2 - 1;
  i3 = emf->nbe + emf->n3 - 1;
  corner_id[6] = i1 + emf->n1pad*(i2 + emf->n2pad*i3);
  i1 = emf->nbe + emf->n1 - 1;
  i2 = emf->nbe + emf->n2 - 1;
  i3 = emf->nbe + emf->n3 - 1;
  corner_id[7] = i1 + emf->n1pad*(i2 + emf->n2pad*i3);
  cudaMemcpy(d_corner_id, corner_id, 8*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemset(d_backup, 0, 8*sizeof(cuFloatComplex));
  cudaMemcpy(d_expfactor, &emf->expfactor[0][0], emf->nfreq*emf->nt*sizeof(cuFloatComplex), cudaMemcpyHostToDevice);
    
  h_chsrc = alloc1int(nchsrc);
  h_chrec = alloc1int(nchrec);
  for(ic=0; ic<nchsrc; ++ic) {
    if     (strcmp(emf->chsrc[ic],"Ex")==0) h_chsrc[ic] = 1;
    else if(strcmp(emf->chsrc[ic],"Ey")==0) h_chsrc[ic] = 2;
    else if(strcmp(emf->chsrc[ic],"Ez")==0) h_chsrc[ic] = 3;
    else if(strcmp(emf->chsrc[ic],"Hx")==0) h_chsrc[ic] = 4;
    else if(strcmp(emf->chsrc[ic],"Hy")==0) h_chsrc[ic] = 5;
    else if(strcmp(emf->chsrc[ic],"Hz")==0) h_chsrc[ic] = 6;
  }
  for(ic=0; ic<nchrec; ++ic) {
    if     (strcmp(emf->chrec[ic],"Ex")==0) h_chrec[ic] = 1;
    else if(strcmp(emf->chrec[ic],"Ey")==0) h_chrec[ic] = 2;
    else if(strcmp(emf->chrec[ic],"Ez")==0) h_chrec[ic] = 3;
    else if(strcmp(emf->chrec[ic],"Hx")==0) h_chrec[ic] = 4;
    else if(strcmp(emf->chrec[ic],"Hy")==0) h_chrec[ic] = 5;
    else if(strcmp(emf->chrec[ic],"Hz")==0) h_chrec[ic] = 6;
  }
  cudaMalloc(&d_chsrc, nchsrc*sizeof(int));
  cudaMalloc(&d_chrec, nchrec*sizeof(int));
  cudaMemcpy(d_chsrc, h_chsrc, nchsrc*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_chrec, h_chrec, nchrec*sizeof(int), cudaMemcpyHostToDevice);
  free(h_chsrc);
  free(h_chrec);

  status = cudaGetLastError();
  if (cudaSuccess!=status) { printf("Failed to initialize memory on device - fdtd !\n"); exit(0); }

}

void cuda_fdtd_close()
{
  /* free memory on device */
  cudaFree(d_inveps11);
  cudaFree(d_inveps22);
  cudaFree(d_inveps33);
  cudaFree(d_E1);
  cudaFree(d_E2);
  cudaFree(d_E3);
  cudaFree(d_H1);
  cudaFree(d_H2);
  cudaFree(d_H3);
  cudaFree(d_curlE1);
  cudaFree(d_curlE2);
  cudaFree(d_curlE3);
  cudaFree(d_curlH1);
  cudaFree(d_curlH2);
  cudaFree(d_curlH3);
  cudaFree(d_memD2E1);
  cudaFree(d_memD3E1);
  cudaFree(d_memD1E2);
  cudaFree(d_memD3E2);
  cudaFree(d_memD1E3);
  cudaFree(d_memD2E3);
  cudaFree(d_memD2H1);
  cudaFree(d_memD3H1);
  cudaFree(d_memD1H2);
  cudaFree(d_memD3H2);
  cudaFree(d_memD1H3);
  cudaFree(d_memD2H3);
  cudaFree(d_fwd_E1);
  cudaFree(d_fwd_E2);
  cudaFree(d_fwd_E3);
  cudaFree(d_a1);
  cudaFree(d_b1);
  cudaFree(d_a2);
  cudaFree(d_b2);
  cudaFree(d_a3);
  cudaFree(d_b3);

  cudaFree(d_corner_id);
  cudaFreeHost(h_ncorner);
  cudaFree(d_backup);
  cudaFree(d_expfactor);
  
  cudaFree(d_chsrc);
  cudaFree(d_chrec);

}


void cuda_airwave_bc_init(emf_t *emf)
{
  int n1fft = emf->n1fft;
  int n2fft = emf->n2fft;
  
  // create FFT plan
  cufftPlan2d(&fftPlan, n1fft, n2fft, CUFFT_C2C);
  cudaMalloc(&d_sH1kxky, n1fft*n2fft*emf->rd*sizeof(cuFloatComplex));
  cudaMalloc(&d_sH2kxky, n1fft*n2fft*emf->rd*sizeof(cuFloatComplex));
  cudaMalloc(&d_sE12kxky, n1fft*n2fft*(emf->rd-1)*sizeof(float));
  cudaMalloc(&d_emfft, n1fft*n2fft*sizeof(cuFloatComplex));
  cudaMalloc(&d_emfft0, n1fft*n2fft*sizeof(cuFloatComplex));
  status = cudaGetLastError();
  if (cudaSuccess!=status) { printf("Failed to allocate memory on device - airwave !\n"); exit(0); }

  
  cudaMemcpy(d_sH1kxky, &emf->sH1kxky[0][0][0], n1fft*n2fft*emf->rd*sizeof(float _Complex), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sH2kxky, &emf->sH2kxky[0][0][0], n1fft*n2fft*emf->rd*sizeof(float _Complex), cudaMemcpyHostToDevice);
  if(emf->rd>1) cudaMemcpy(d_sE12kxky, &emf->sE12kxky[0][0][0], n1fft*n2fft*(emf->rd-1)*sizeof(float), cudaMemcpyHostToDevice);
  status = cudaGetLastError();
  if (cudaSuccess!=status) { printf("Failed to initialize memory on device - airwave!\n"); exit(0); }
  
}


void cuda_airwave_bc_close()
{
  cufftDestroy(fftPlan);
  cudaFree(d_sH1kxky);
  cudaFree(d_sH2kxky);
  cudaFree(d_sE12kxky);
  cudaFree(d_emfft);
  cudaFree(d_emfft0);
  
}



void cuda_interpolation_init(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg)
{
  int rd = emf->rd;

  cudaMalloc(&d_rg_src_i1, acqui->nsrc*acqui->nsubsrc*sizeof(int));
  cudaMalloc(&d_rg_src_i2, acqui->nsrc*acqui->nsubsrc*sizeof(int));
  cudaMalloc(&d_rg_src_i3, acqui->nsrc*acqui->nsubsrc*sizeof(int));
  cudaMalloc(&d_rg_src_w1, 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float));
  cudaMalloc(&d_rg_src_w2, 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float));
  cudaMalloc(&d_rg_src_w3, 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float));

  cudaMalloc(&d_sg_src_i1, acqui->nsrc*acqui->nsubsrc*sizeof(int));
  cudaMalloc(&d_sg_src_i2, acqui->nsrc*acqui->nsubsrc*sizeof(int));
  cudaMalloc(&d_sg_src_i3, acqui->nsrc*acqui->nsubsrc*sizeof(int));
  cudaMalloc(&d_sg_src_w1, 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float));
  cudaMalloc(&d_sg_src_w2, 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float));
  cudaMalloc(&d_sg_src_w3, 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float));
  
  cudaMalloc(&d_rg_rec_i1, acqui->nrec*acqui->nsubrec*sizeof(int));
  cudaMalloc(&d_rg_rec_i2, acqui->nrec*acqui->nsubrec*sizeof(int));
  cudaMalloc(&d_rg_rec_i3, acqui->nrec*acqui->nsubrec*sizeof(int));
  cudaMalloc(&d_rg_rec_w1, 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float));
  cudaMalloc(&d_rg_rec_w2, 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float));
  cudaMalloc(&d_rg_rec_w3, 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float));

  cudaMalloc(&d_sg_rec_i1, acqui->nrec*acqui->nsubrec*sizeof(int));
  cudaMalloc(&d_sg_rec_i2, acqui->nrec*acqui->nsubrec*sizeof(int));
  cudaMalloc(&d_sg_rec_i3, acqui->nrec*acqui->nsubrec*sizeof(int));
  cudaMalloc(&d_sg_rec_w1, 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float));
  cudaMalloc(&d_sg_rec_w2, 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float));
  cudaMalloc(&d_sg_rec_w3, 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float));
  status = cudaGetLastError();
  if (cudaSuccess!=status) { printf("Failed to allocate memory on device - interpolation!\n"); exit(0);  }

  //-------------------------------------------------------------------------
  cudaMemcpy(d_rg_src_i1, interp_rg->src_i1[0], acqui->nsrc*acqui->nsubsrc*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rg_src_i2, interp_rg->src_i2[0], acqui->nsrc*acqui->nsubsrc*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rg_src_i3, interp_rg->src_i3[0], acqui->nsrc*acqui->nsubsrc*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rg_src_w1, interp_rg->src_w1[0][0], 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rg_src_w2, interp_rg->src_w2[0][0], 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rg_src_w3, interp_rg->src_w3[0][0], 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_sg_src_i1, interp_sg->src_i1[0], acqui->nsrc*acqui->nsubsrc*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sg_src_i2, interp_sg->src_i2[0], acqui->nsrc*acqui->nsubsrc*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sg_src_i3, interp_sg->src_i3[0], acqui->nsrc*acqui->nsubsrc*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sg_src_w1, interp_sg->src_w1[0][0], 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sg_src_w2, interp_sg->src_w2[0][0], 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sg_src_w3, interp_sg->src_w3[0][0], 2*rd*acqui->nsrc*acqui->nsubsrc*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_rg_rec_i1, interp_rg->rec_i1[0], acqui->nrec*acqui->nsubrec*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rg_rec_i2, interp_rg->rec_i2[0], acqui->nrec*acqui->nsubrec*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rg_rec_i3, interp_rg->rec_i3[0], acqui->nrec*acqui->nsubrec*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rg_rec_w1, interp_rg->rec_w1[0][0], 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rg_rec_w2, interp_rg->rec_w2[0][0], 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rg_rec_w3, interp_rg->rec_w3[0][0], 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_sg_rec_i1, interp_sg->rec_i1[0], acqui->nrec*acqui->nsubrec*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sg_rec_i2, interp_sg->rec_i2[0], acqui->nrec*acqui->nsubrec*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sg_rec_i3, interp_sg->rec_i3[0], acqui->nrec*acqui->nsubrec*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sg_rec_w1, interp_sg->rec_w1[0][0], 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sg_rec_w2, interp_sg->rec_w2[0][0], 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sg_rec_w3, interp_sg->rec_w3[0][0], 2*rd*acqui->nrec*acqui->nsubrec*sizeof(float), cudaMemcpyHostToDevice);

  status = cudaGetLastError();
  if(cudaSuccess!=status) { printf("Failed to initialize memory on device - interpolation !\n"); exit(0); }
}

void cuda_interpolation_close()
{
  cudaFree(d_rg_src_i1);
  cudaFree(d_rg_src_i2);
  cudaFree(d_rg_src_i3);
  cudaFree(d_rg_src_w1);
  cudaFree(d_rg_src_w2);
  cudaFree(d_rg_src_w3);

  cudaFree(d_sg_src_i1);
  cudaFree(d_sg_src_i2);
  cudaFree(d_sg_src_i3);
  cudaFree(d_sg_src_w1);
  cudaFree(d_sg_src_w2);
  cudaFree(d_sg_src_w3);

  cudaFree(d_rg_rec_i1);
  cudaFree(d_rg_rec_i2);
  cudaFree(d_rg_rec_i3);
  cudaFree(d_rg_rec_w1);
  cudaFree(d_rg_rec_w2);
  cudaFree(d_rg_rec_w3);

  cudaFree(d_sg_rec_i1);
  cudaFree(d_sg_rec_i2);
  cudaFree(d_sg_rec_i3);
  cudaFree(d_sg_rec_w1);
  cudaFree(d_sg_rec_w2);
  cudaFree(d_sg_rec_w3);
}


extern "C"
void cuda_modeling(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg, int adj)
{
  double t_update_E,t_update_H,t_inject_E,t_inject_H,t_curlE, t_curlH,t_dft_emf,t0,t_convergence;
  static int doneinit = 0;

  if(!doneinit){
    cudaSetDevice(0);// initialize device, default device=0;
    status = cudaGetLastError();
    if (cudaSuccess!=status) { printf("Failed to initialize device!\n"); exit(0); }
    doneinit = 1;
  }

  /*==========================================================*/
  cuda_interpolation_init(acqui, emf, interp_rg, interp_sg);
  cuda_fdtd_init(emf);
  if(emf->airwave) cuda_airwave_bc_init(emf);

  if(adj){
    cudaMalloc(&d_dres_td, emf->nt*acqui->nrec*emf->nchrec*sizeof(float));
    cudaMemcpy(d_dres_td, &emf->dres_td[0][0][0], emf->nt*acqui->nrec*emf->nchrec*sizeof(float), cudaMemcpyHostToDevice);
    status = cudaGetLastError();
    if (cudaSuccess!=status) { printf("Failed to allocate memory on device - fdtd !\n"); exit(0); }
    
  }
    
  dimBlock.x = BlockSize1; 
  dimBlock.y = BlockSize2;
  dimGrid.x = (emf->n1pad+BlockSize1-1)/BlockSize1;
  dimGrid.y = (emf->n2pad+BlockSize2-1)/BlockSize2;
  dimGrid_dft.x = (emf->n1fft+BlockSize1-1)/BlockSize1;
  dimGrid_dft.y = (emf->n2fft+BlockSize2-1)/BlockSize2;
  if(emf->verb){
    printf("dimBlock.x=%d \n", dimBlock.x);
    printf("dimBlock.y=%d \n", dimBlock.y);
    printf("dimGrid.x=%d \n", dimGrid.x);
    printf("dimGrid.y=%d \n", dimGrid.y);
    printf("dimGrid_dft.x=%d \n", dimGrid_dft.x);
    printf("dimGrid_dft.y=%d \n", dimGrid_dft.y);

    t0 = 0;
    t_curlE = 0.;
    t_inject_H = 0.;
    t_update_H = 0.;
    t_curlH = 0.;
    t_inject_E = 0.;
    t_update_E = 0.;
    t_dft_emf= 0.;
    t_convergence=0.;
  }

  int it;
  float mstimer;
  cudaEvent_t start, stop;

  cudaEventCreate(&start);	
  cudaEventCreate(&stop);
  cudaEventRecord(start);
  for(it=0; it<emf->nt; it++){
    if(emf->verb && it%50==0) printf("it---- %d\n", it);

    /*--------------------------------------------------------------*/
    if(emf->verb) t0 = MPI_Wtime();
    cuda_fdtd_curlE<<<dimGrid,dimBlock>>>
      (d_E1, d_E2, d_E3, d_curlE1, d_curlE2, d_curlE3, d_a1, d_b1, d_a2, d_b2, d_a3, d_b3,
       d_memD1E2, d_memD1E3, d_memD2E1, d_memD2E3, d_memD3E1, d_memD3E2, 
       c11, c21, c12, c22, c13, c23,
       emf->n1pad, emf->n2pad, emf->n3pad, emf->nb, emf->nbe, emf->airwave, 
       adj?emf->i3min_adj[it]:emf->i3min_fwd[it], adj?emf->i3max_adj[it]:emf->i3max_fwd[it]);
    status = cudaGetLastError();
    if (cudaSuccess!=status) { printf("Failed to compute curlE on device!\n"); exit(0); }    
    if(emf->verb) t_curlE += MPI_Wtime()-t0;

    
    if(emf->verb) t0 = MPI_Wtime();
    cuda_inject_magnetic_source<<<(emf->nchsrc*acqui->nsrc+BlockSize-1)/BlockSize,BlockSize>>>
      (d_rg_src_i1, d_rg_src_i2, d_rg_src_i3, d_rg_src_w1, d_rg_src_w2, d_rg_src_w3,
       d_sg_src_i1, d_sg_src_i2, d_sg_src_i3, d_sg_src_w1, d_sg_src_w2, d_sg_src_w3,
       d_curlE1, d_curlE2, d_curlE3, d_chsrc,
       emf->stf[it], emf->d1, emf->d2, emf->d3, emf->nchsrc, acqui->nsrc, acqui->nsubsrc,
       emf->n1pad, emf->n2pad, emf->n3pad, emf->nbe, emf->rd);
    status = cudaGetLastError();
    if (cudaSuccess!=status) { printf("Failed to inject magnetic source on device!\n"); exit(0); }
    if(emf->verb) t_inject_H += MPI_Wtime()-t0;

    
    if(emf->verb) t0 = MPI_Wtime();
    cuda_fdtd_update_H<<<dimGrid,dimBlock>>>(d_H1, d_H2, d_H3, d_curlE1, d_curlE2, d_curlE3,
					     emf->dt, emf->n1pad, emf->n2pad, emf->n3pad);
    status = cudaGetLastError();
    if (cudaSuccess!=status) { printf("Failed to update H on device!\n"); exit(0); }    
    if(emf->airwave){
      cuda_airwave_bc_copy<<<dimGrid_dft,dimBlock>>>
    	(d_emfft, &d_H3[emf->n1pad*emf->n2pad*emf->nbe], emf->n1pad, emf->n2pad, emf->n1fft, emf->n2fft);
      cufftExecC2C(fftPlan, d_emfft, d_emfft0, CUFFT_FORWARD);//FFT into wavenumber domain

      cudaMemcpy(d_emfft, d_emfft0, emf->n1fft*emf->n2fft*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice);
      cuda_airwave_bc_scale_FH<<<dimGrid_dft,dimBlock>>>(d_emfft, d_sH1kxky, emf->n1fft, emf->n2fft);
      cufftExecC2C(fftPlan, d_emfft, d_emfft, CUFFT_INVERSE);//IFFT back to space domain
      cuda_airwave_bc_back2emf<<<dimGrid_dft,dimBlock>>>
	(&d_H1[emf->n1pad*emf->n2pad*(emf->nbe-1)], d_emfft, emf->n1pad, emf->n2pad, emf->n1fft, emf->n2fft);

      cudaMemcpy(d_emfft, d_emfft0, emf->n1fft*emf->n2fft*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice);
      cuda_airwave_bc_scale_FH<<<dimGrid_dft,dimBlock>>>(d_emfft, &d_sH1kxky[emf->n1fft*emf->n2fft], emf->n1fft, emf->n2fft);
      cufftExecC2C(fftPlan, d_emfft, d_emfft, CUFFT_INVERSE);//IFFT back to space domain
      cuda_airwave_bc_back2emf<<<dimGrid_dft,dimBlock>>>
	(&d_H1[emf->n1pad*emf->n2pad*(emf->nbe-2)], d_emfft, emf->n1pad, emf->n2pad, emf->n1fft, emf->n2fft);

      cudaMemcpy(d_emfft, d_emfft0, emf->n1fft*emf->n2fft*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice);
      cuda_airwave_bc_scale_FH<<<dimGrid_dft,dimBlock>>>(d_emfft, d_sH2kxky, emf->n1fft, emf->n2fft);
      cufftExecC2C(fftPlan, d_emfft, d_emfft, CUFFT_INVERSE);//IFFT back to space domain
      cuda_airwave_bc_back2emf<<<dimGrid_dft,dimBlock>>>
	(&d_H2[emf->n1pad*emf->n2pad*(emf->nbe-1)], d_emfft, emf->n1pad, emf->n2pad, emf->n1fft, emf->n2fft);

      cudaMemcpy(d_emfft, d_emfft0, emf->n1fft*emf->n2fft*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice);
      cuda_airwave_bc_scale_FH<<<dimGrid_dft,dimBlock>>>(d_emfft, &d_sH2kxky[emf->n1fft*emf->n2fft], emf->n1fft, emf->n2fft);
      cufftExecC2C(fftPlan, d_emfft, d_emfft, CUFFT_INVERSE);//IFFT back to space domain
      cuda_airwave_bc_back2emf<<<dimGrid_dft,dimBlock>>>
	(&d_H2[emf->n1pad*emf->n2pad*(emf->nbe-2)], d_emfft, emf->n1pad, emf->n2pad, emf->n1fft, emf->n2fft);

      status = cudaGetLastError();
      if (cudaSuccess!=status) { printf("Failed to handle air-water interface!\n"); exit(0); }
    }
    if(emf->verb) t_update_H += MPI_Wtime()-t0;


    /*--------------------------------------------------------------*/
    if(emf->verb) t0 = MPI_Wtime();
    cuda_fdtd_curlH<<<dimGrid,dimBlock>>>(d_H1, d_H2, d_H3, d_curlH1, d_curlH2, d_curlH3,
    					  d_a1, d_b1, d_a2, d_b2, d_a3, d_b3,
    					  d_memD1H2, d_memD1H3, d_memD2H1,
    					  d_memD2H3, d_memD3H1, d_memD3H2, 
    					  d_inveps11, d_inveps22, d_inveps33,
    					  c11, c21, c12, c22, c13, c23,
    					  emf->n1pad, emf->n2pad, emf->n3pad, emf->nb, emf->nbe, emf->airwave,
					  adj?emf->i3min_adj[it]:emf->i3min_fwd[it], adj?emf->i3max_adj[it]:emf->i3max_fwd[it]);
    status = cudaGetLastError();
    if (cudaSuccess!=status) { printf("Failed to compute curlH on device!\n"); exit(0); }    
    if(emf->verb) t_curlH += MPI_Wtime()-t0;

    if(emf->verb) t0 = MPI_Wtime();
    if(adj){
      cuda_inject_electric_src_adj<<<(emf->nchrec*acqui->nrec+BlockSize-1)/BlockSize,BlockSize>>>
	(d_rg_rec_i1, d_rg_rec_i2, d_rg_rec_i3, d_rg_rec_w1, d_rg_rec_w2, d_rg_rec_w3,
	 d_sg_rec_i1, d_sg_rec_i2, d_sg_rec_i3, d_sg_rec_w1, d_sg_rec_w2, d_sg_rec_w3,
	 d_inveps11, d_inveps22, d_inveps33, d_curlH1, d_curlH2, d_curlH3, d_chrec, 
	 d_dres_td, it, emf->nt, emf->nchrec, acqui->nrec, acqui->nsubrec,
	 emf->n1pad, emf->n2pad, emf->n3pad, emf->nbe, emf->rd);
    }else{
      cuda_inject_electric_src_fwd<<<(emf->nchsrc*acqui->nsrc+BlockSize-1)/BlockSize,BlockSize>>>
	(d_rg_src_i1, d_rg_src_i2, d_rg_src_i3, d_rg_src_w1, d_rg_src_w2, d_rg_src_w3,
	 d_sg_src_i1, d_sg_src_i2, d_sg_src_i3, d_sg_src_w1, d_sg_src_w2, d_sg_src_w3,
	 d_inveps11, d_inveps22, d_inveps33, d_curlH1, d_curlH2, d_curlH3, d_chsrc, 
	 emf->stf[it], emf->d1, emf->d2, emf->d3, emf->nchsrc, acqui->nsrc, acqui->nsubsrc,
	 emf->n1pad, emf->n2pad, emf->n3pad, emf->nbe, emf->rd);
    }
    status = cudaGetLastError();
    if (cudaSuccess!=status) { printf("Failed to inject electric source on device!\n"); exit(0); }
    if(emf->verb) t_inject_E += MPI_Wtime()-t0;

    
    if(emf->verb) t0 = MPI_Wtime();
    cuda_fdtd_update_E<<<dimGrid,dimBlock>>>
      (d_E1, d_E2, d_E3, d_curlH1, d_curlH2, d_curlH3, d_inveps11, d_inveps22, d_inveps33,
       emf->n1pad, emf->n2pad, emf->n3pad, emf->dt);
    status = cudaGetLastError();
    if (cudaSuccess!=status) { printf("Failed to update E on device!\n"); exit(0); }    
    if(emf->airwave){
      cuda_airwave_bc_copy<<<dimGrid_dft,dimBlock>>>
	(d_emfft, &d_E1[emf->n1pad*emf->n2pad*emf->nbe], emf->n1pad, emf->n2pad, emf->n1fft, emf->n2fft);
      cufftExecC2C(fftPlan, d_emfft, d_emfft, CUFFT_FORWARD);//FFT into wavenumber domain
      cuda_airwave_bc_scale_FE<<<dimGrid_dft,dimBlock>>>(d_emfft, d_sE12kxky, emf->n1fft, emf->n2fft);
      cufftExecC2C(fftPlan, d_emfft, d_emfft, CUFFT_INVERSE);//IFFT back to space domain
      cuda_airwave_bc_back2emf<<<dimGrid,dimBlock>>>
    	(&d_E1[emf->n1pad*emf->n2pad*(emf->nbe-1)], d_emfft, emf->n1pad, emf->n2pad, emf->n1fft, emf->n2fft);

      cuda_airwave_bc_copy<<<dimGrid_dft,dimBlock>>>
	(d_emfft, &d_E2[emf->n1pad*emf->n2pad*emf->nbe], emf->n1pad, emf->n2pad, emf->n1fft, emf->n2fft);
      cufftExecC2C(fftPlan, d_emfft, d_emfft, CUFFT_FORWARD);//FFT into wavenumber domain
      cuda_airwave_bc_scale_FE<<<dimGrid_dft,dimBlock>>>(d_emfft, d_sE12kxky, emf->n1fft, emf->n2fft);
      cufftExecC2C(fftPlan, d_emfft, d_emfft, CUFFT_INVERSE);//IFFT back to space domain
      cuda_airwave_bc_back2emf<<<dimGrid,dimBlock>>>
    	(&d_E2[emf->n1pad*emf->n2pad*(emf->nbe-1)], d_emfft, emf->n1pad, emf->n2pad, emf->n1fft, emf->n2fft);

      status = cudaGetLastError();
      if (cudaSuccess!=status) { printf("Failed to handle air-water interface!\n"); exit(0); }
    }
    if(emf->verb) t_update_E += MPI_Wtime()-t0;
    
    /*--------------------------------------------------------------*/
    if(emf->verb) t0 = MPI_Wtime();
    cuda_dtft_emf<<<dimGrid,dimBlock>>>(d_fwd_E1, d_fwd_E2, d_fwd_E3, &d_expfactor[it*emf->nfreq], d_E1, d_E2, d_E3, 
					emf->nb, emf->n123pad, emf->n1pad, emf->n2pad, emf->n3pad, emf->nfreq,
					adj?emf->i3min_adj[it]:emf->i3min_fwd[it], adj?emf->i3max_adj[it]:emf->i3max_fwd[it]);
    status = cudaGetLastError();
    if (cudaSuccess!=status) { printf("Failed to compute DFT of E + H on device!\n"); exit(0); }    
    if(emf->verb) t_dft_emf += MPI_Wtime()-t0;


    /*--------------------------------------------------------------*/
    if(emf->verb) t0 = MPI_Wtime();
    if(it%100==0){/* convergence check */
      cuda_check_convergence1<<<1,8>>>(d_corner_id, d_fwd_E1, d_backup, h_ncorner);
      if(emf->verb) printf("%d corners of the cube converged!\n", h_ncorner[0]);
      if(h_ncorner[0]==8) { emf->nt = it; printf("converge after %d steps\n", it); break; }/* all 8 corners converged, exit now */
    }
    if(emf->verb) t_convergence += MPI_Wtime()-t0;
  }

  
  if(emf->verb) {
    t0 = t_curlH + t_inject_E + t_update_E + t_curlE + t_inject_H + t_update_H
      + t_dft_emf + t_convergence;
    printf("-------------- elapsed time --------------------\n");
    printf("    compute curlE: %e s\n", t_curlE);
    printf("    inject magnetic source: %e s\n", t_inject_H);
    printf("    update magnetic field: %e s\n", t_update_H);

    printf("    compute curlH: %e s\n", t_curlH);
    printf("    inject electric source: %e s\n", t_inject_E);
    printf("    update electric field: %e s\n", t_update_E);
    
    printf("    DFT EM field: %e s\n", t_dft_emf);
    printf("    convergence check: %e s\n", t_convergence);
    printf("    Total modeling time: %e s\n", t0);
    printf("------------------------------------------------\n");
  }
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&mstimer, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  if(emf->verb) printf("Elapsed time: %g (s)\n", mstimer*1.e-3);

  if(adj){
    cudaMemcpy(&emf->adj_E1[0][0][0][0], &d_fwd_E1[0], emf->nfreq*emf->n123pad*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&emf->adj_E2[0][0][0][0], &d_fwd_E2[0], emf->nfreq*emf->n123pad*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&emf->adj_E3[0][0][0][0], &d_fwd_E3[0], emf->nfreq*emf->n123pad*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
  }else{
    cudaMemcpy(&emf->fwd_E1[0][0][0][0], &d_fwd_E1[0], emf->nfreq*emf->n123pad*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&emf->fwd_E2[0][0][0][0], &d_fwd_E2[0], emf->nfreq*emf->n123pad*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&emf->fwd_E3[0][0][0][0], &d_fwd_E3[0], emf->nfreq*emf->n123pad*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
  }

  cuda_interpolation_close();
  cuda_fdtd_close(); 
  if(emf->airwave) cuda_airwave_bc_close();

  if(adj) cudaFree(d_dres_td);
}
