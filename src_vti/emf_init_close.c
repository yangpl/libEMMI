/* initialize important parameter for modelling electromagnetic fields (emf)
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"
#include "mpi_info.h"
#include "constants.h"


int cmpfunc(const void *a, const void *b) { return ( *(int*)a - *(int*)b ); }

void nugrid_init(emf_t *emf);
void nugrid_close(emf_t *emf);

void homogenization(emf_t *emf);

void emf_init(emf_t *emf)
{
  char *frho_h, *frho_v;
  FILE *fp=NULL;
  int ifreq, ic, istat;

  if(!getparint("reciprocity", &emf->reciprocity)) emf->reciprocity = 0; 

  if(!getparint("n1", &emf->n1)) emf->n1 = 101; 
  /* number of cells in axis-1, nx */
  if(!getparint("n2", &emf->n2)) emf->n2 = 101; 
  /* number of cells in axis-2, ny */
  if(!getparint("n3", &emf->n3)) emf->n3 = 51; 
  /* number of cells in axis-3, nz */
  if(!getparint("nb", &emf->nb)) emf->nb = 12; 
  /* number of PML layers on each side */
  if(!getparfloat("d1", &emf->d1)) emf->d1 = 100.;
  /* grid spacing in 1st dimension, dx */
  if(!getparfloat("d2", &emf->d2)) emf->d2 = 100.;
  /* grid spacing in 2nd dimension, dy */
  if(!getparfloat("d3", &emf->d3)) emf->d3 = 100.;
  /* grid spacing in 3rd dimension, dz */
  if(!getparint("rd1", &emf->rd1)) emf->rd1 = 2; 
  /* half length of FD stencil */
  if(!getparint("rd2", &emf->rd2)) emf->rd2 = 2; 
  /* half length of FD stencil */
  if(!getparint("rd3", &emf->rd3)) emf->rd3 = 1;
#ifdef GPU
  emf->rd3 = 2;//GPU code does not support rd3=1 yet
#endif
  /* half length of FD stencil */
  if(!getparint("ne", &emf->ne)) emf->ne = MAX(MAX(emf->rd1, emf->rd2), emf->rd3)+5;
  /* number of buffer layers on each side */
  if(!getparint("airwave", &emf->airwave)) emf->airwave = 1; 
  /* simulate airwave on top boundary */

  emf->n123 = emf->n1*emf->n2*emf->n3;
  emf->nbe = emf->nb + emf->ne;/* number of PML layers + extra 2 points due to 4-th order FD */
  emf->n1pad = emf->n1+2*emf->nbe;/* total number of grid points after padding PML+extra */
  emf->n2pad = emf->n2+2*emf->nbe;/* total number of grid points after padding PML+extra */
  emf->n3pad = emf->n3+2*emf->nbe;/* total number of grid points after padding PML+extra */
  emf->n123pad = emf->n1pad*emf->n2pad*emf->n3pad;
  if(iproc==0){
    printf("reciprocity=%d\n", emf->reciprocity);
    printf("PML layers on each side: nb=%d\n", emf->nb);
    printf("Buffer layers on each side: ne=%d\n", emf->ne);
    printf("Number of layers extended outside model: nbe=%d\n", emf->nbe);
    printf("[d1, d2, d3]=[%g, %g, %g]\n", emf->d1, emf->d2, emf->d3);
    printf("[n1, n2, n3]=[%d, %d, %d]\n", emf->n1, emf->n2, emf->n3);
    printf("[n1pad, n2pad, n3pad]=[%d, %d, %d]\n", emf->n1pad, emf->n2pad, emf->n3pad);
    printf("[rd1, rd2, rd3]=[%d, %d, %d]\n", emf->rd1, emf->rd2, emf->rd3);
  }
  
  if(!getparfloat("f0", &emf->f0)) emf->f0 = 1.;
  emf->omega0 = 2.*PI*emf->f0;
  /* reference frequency */
  if(!getparfloat("Glim", &emf->Glim)) emf->Glim = 5;/* 5 points/wavelength */
  if(!(getparstring("frho_h", &frho_h))) err("Need frho_h= ");
  if(!(getparstring("frho_v", &frho_v))) err("Need frho_v= ");

  
  if(!(emf->nfreq=countparval("freqs"))) err("Need freqs= vector");
  /* number of frequencies for electromagnetic emf->modeling */
  emf->freqs=alloc1float(emf->nfreq);
  emf->omegas=alloc1float(emf->nfreq);
  getparfloat("freqs", emf->freqs);/* a list of frequencies separated by comma */
  qsort(emf->freqs, emf->nfreq, sizeof(float), cmpfunc);/*sort frequencies in ascending order*/
  for(ifreq=0; ifreq<emf->nfreq; ++ifreq) {
    emf->omegas[ifreq]=2.*PI*emf->freqs[ifreq];
    if(iproc==0) printf("freq[%d]=%g ", ifreq+1, emf->freqs[ifreq]);
  }
  if(iproc==0) printf("\n");
  
  /* read active source channels */
  if((emf->nchsrc=countparval("chsrc"))!=0) {
    emf->chsrc=(char**)alloc1(emf->nchsrc, sizeof(void*));
    getparstringarray("chsrc", emf->chsrc);
    /* active source channels: Ex, Ey, Ez, Hx, Hy, Hz or their combinations */
  }else{
    emf->nchsrc=1;
    emf->chsrc=(char**)alloc1(emf->nchsrc, sizeof(void*));
    emf->chsrc[0]="Ex";
  }
  /* read active receiver channels */
  if((emf->nchrec=countparval("chrec"))!=0) {
    emf->chrec=(char**)alloc1(emf->nchrec, sizeof(void*));
    getparstringarray("chrec", emf->chrec);
    /* active receiver channels: Ex, Ey, Ez, Hx, Hy, Hz or their combinations */
  }else{
    emf->nchrec=2;
    emf->chrec=(char**)alloc1(emf->nchrec, sizeof(void*));
    emf->chrec[0] = "Ex";
    emf->chrec[1] = "Ey";
  }
  if(iproc==0){
    printf("Active source channels:");
    for(ic=0; ic<emf->nchsrc; ++ic) printf(" %s", emf->chsrc[ic]);
    printf("\n");
    printf("Active recever channels:");
    for(ic=0; ic<emf->nchrec; ++ic) printf(" %s", emf->chrec[ic]);
    printf("\n");
  }

  /*-------------------------------------------------------*/
  emf->rho_h = alloc3float(emf->n1, emf->n2, emf->n3);//effective rho used in inversion
  emf->rho_v = alloc3float(emf->n1, emf->n2, emf->n3);//effective rho used in inversion
  emf->rho11 = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->rho22 = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->rho33 = alloc3float(emf->n1, emf->n2, emf->n3);

  
  /* read  resistivity */
  fp = fopen(frho_h, "rb");
  if(fp==NULL) err("cannot open file");
  istat = fread(emf->rho_h[0][0], sizeof(float), emf->n123, fp);
  if(istat != emf->n123) err("size parameter does not match the file!");
  fclose(fp);
   
  fp = fopen(frho_v, "rb");
  if(fp==NULL) err("cannot open file");
  istat = fread(emf->rho_v[0][0], sizeof(float), emf->n123, fp);
  if(istat != emf->n123) err("size parameter does not match the file!");
  fclose(fp);

  homogenization(emf);/*create rho11, rho22, rho33 by homogenization along x, y, z */

  if(!getparint("nugrid", &emf->nugrid)) emf->nugrid=0;/* 1=nonuniform grid; 0=uniform grid */
  if(iproc==0) printf("nugrid=%d\n", emf->nugrid);
  if(emf->nugrid) nugrid_init(emf);
}



void emf_close(emf_t *emf)
{
  free1float(emf->freqs);
  free1float(emf->omegas);
  free3float(emf->rho_h);
  free3float(emf->rho_v);
  free3float(emf->rho11); 
  free3float(emf->rho22);
  free3float(emf->rho33);

  if(emf->nugrid) nugrid_close(emf);
}

