/* read and write CSEM data in binary
 *
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 */
#include <unistd.h>
#include "cstd.h"
#include "mpi_info.h"

#include "acqui.h"
#include "emf.h"

void write_data(acqui_t *acqui, emf_t *emf, char *fname, float _Complex ***dcal_fd)
/*< write synthetic data according to shot/process index >*/
{
  FILE *fp;
  int isrc, irec, ichrec, ifreq;
  float dp_re, dp_im;

  fp=fopen(fname,"w");
  if(fp==NULL) err("error opening file for writing");
  fprintf(fp, "iTx 	 iRx    chrec  ifreq 	 emf_real 	 emf_imag\n");
  isrc = acqui->shot_idx[iproc];//index starts from 1
  for(ichrec=0; ichrec<emf->nchrec; ichrec++){
    for(ifreq=0; ifreq<emf->nfreq; ifreq++){
      for(irec=0; irec<acqui->nrec; irec++){
	dp_re = creal(dcal_fd[ichrec][ifreq][irec]);
	dp_im = cimag(dcal_fd[ichrec][ifreq][irec]);
	fprintf(fp, "%d \t %d \t %s \t %d \t %e \t %e\n",
		isrc, irec+1, emf->chrec[ichrec], ifreq+1, dp_re, dp_im);
      }
    }
  }
  fclose(fp);
}


void read_data(acqui_t *acqui, emf_t *emf)
/*< read observed data according to shot/process index >*/
{
  char fname[sizeof("emf_0000.txt")];
  int isrc, irec, ifreq, iseof;
  float dp_re, dp_im;
  char chrec[sizeof("Ex")];
  FILE *fp;

  memset(&emf->Ea[0][0], 0, emf->nfreq*acqui->nrec*sizeof(float _Complex));
  memset(&emf->Eb[0][0], 0, emf->nfreq*acqui->nrec*sizeof(float _Complex));
  memset(&emf->Ha[0][0], 0, emf->nfreq*acqui->nrec*sizeof(float _Complex));
  memset(&emf->Hb[0][0], 0, emf->nfreq*acqui->nrec*sizeof(float _Complex));

  if(emf->reciprocity){
    
  }else{
    sprintf(fname, "emf_%04d.txt", acqui->shot_idx[iproc]);
    if(access(fname, F_OK)!=0) err("file %s does not exist", fname);
    fp=fopen(fname,"r");
    fscanf(fp, "%*[^\n]\n");//skip a line at the beginning of the file
    while(1){
      iseof=fscanf(fp, "%d \t %d \t %s \t %d \t %e \t %e\n",
		   &isrc, &irec, chrec, &ifreq, &dp_re, &dp_im);
      if(iseof==EOF){
	break;
      }else{
	if(acqui->shot_idx[iproc]==isrc && ifreq<=emf->nfreq && irec<=acqui->nrec){
	  if(strcmp(chrec,"Ex")==0){
	    emf->Ea[ifreq-1][irec-1] = dp_re + I* dp_im;
	  }
	  if(strcmp(chrec,"Ey")==0){
	    emf->Eb[ifreq-1][irec-1] = dp_re + I* dp_im;
	  }
	  if(strcmp(chrec,"Hx")==0){
	    emf->Ha[ifreq-1][irec-1] = dp_re + I* dp_im;
	  }
	  if(strcmp(chrec,"Hy")==0){
	    emf->Hb[ifreq-1][irec-1] = dp_re + I* dp_im;
	  }
	}//end if
      }
    }
    fclose(fp);
  }
    
}

void write_misfit(acqui_t *acqui, emf_t *emf, char *fname, float ***obj)
/*< write synthetic data according to shot/process index >*/
{
  FILE *fp;
  int isrc, irec, ichrec, ifreq;

  fp=fopen(fname,"w");
  if(fp==NULL) err("error opening file for writing");
  fprintf(fp, "iTx 	 iRx    ichrec  ifreq 	 sig_misfit   x       y      z\n");
  isrc = acqui->shot_idx[iproc];//index starts from 1
  for(ichrec=0; ichrec<emf->nchrec; ichrec++){
    for(ifreq=0; ifreq<emf->nfreq; ifreq++){
      for(irec=0; irec<acqui->nrec; irec++){
	fprintf(fp, "%d \t %d \t %s \t %d \t %e \t %e \t %e \t %e\n",
		isrc, irec+1, emf->chrec[ichrec], ifreq+1, obj[ichrec][ifreq][irec],
		acqui->rec_x1[irec], acqui->rec_x2[irec], acqui->rec_x3[irec]);
      }
    }

  }
  fclose(fp);

}

