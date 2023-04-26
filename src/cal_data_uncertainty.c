/* calculate EM data uncertainty using all data components
 *
 *  Reference: Morten 2009 SEG, CSEM data uncertainty analysis for 3D inversion
 *
 * Copyright (c) 2020 Pengliang Yang. All rights reserved.
 * Email: ypl.2100@gmail.com
 */
#include <mpi.h>

#include "cstd.h"
#include "acqui.h"
#include "emf.h"
#include "fwi.h"
#include "mpi_info.h"

float frannor(void);/* random noise generator*/
void write_data(acqui_t *acqui, emf_t *emf, char *fname, float _Complex ***dcal_fd);

void write_misfit(acqui_t *acqui, emf_t *emf, char *fname, float ***obj);

void cal_data_uncertainty(acqui_t *acqui, emf_t *emf, fwi_t *fwi)
{
  int ifreq, isrc, irec, ichrec, ic;
  float d1, d2, tmp,tmp1,tmp2, phi;
  float _Complex cgwn;

  //step 1: compute Fx, Fy using Fa, Fb, where F=E,H
  memset(&emf->dobs_fd[0][0][0], 0, emf->nchrec*emf->nfreq*acqui->nrec*sizeof(float _Complex));
  isrc = 0; //current shot global index=acqui->shot_idx[iproc]-1, local index=0 (always only 1 shot on each processor)
  for(ichrec=0; ichrec<emf->nchrec; ichrec++){
    if(strcmp(emf->chrec[ichrec],"Ex")==0){
      for(irec=0; irec<acqui->nrec; irec++){
	phi = acqui->rec_azimuth[irec]-acqui->src_azimuth[isrc];
	for(ifreq=0; ifreq<emf->nfreq; ifreq++){
	  emf->dobs_fd[ichrec][ifreq][irec] = emf->Ea[ifreq][irec]*cos(phi) + emf->Eb[ifreq][irec]*sin(phi);
	}//end for ifreq
      }//end for irec
	
    }else if(strcmp(emf->chrec[ichrec],"Ey")==0){
      for(irec=0; irec<acqui->nrec; irec++){
	phi = acqui->rec_azimuth[irec]-acqui->src_azimuth[isrc];
	for(ifreq=0; ifreq<emf->nfreq; ifreq++){
	  emf->dobs_fd[ichrec][ifreq][irec] = -emf->Ea[ifreq][irec]*sin(phi) + emf->Eb[ifreq][irec]*cos(phi);
	}//end for ireq
      }//end for irec
	
    }else if(strcmp(emf->chrec[ichrec],"Hx")==0){
      for(irec=0; irec<acqui->nrec; irec++){
	phi = acqui->rec_azimuth[irec]-acqui->src_azimuth[isrc];
	for(ifreq=0; ifreq<emf->nfreq; ifreq++){
	  emf->dobs_fd[ichrec][ifreq][irec] = emf->Ha[ifreq][irec]*cos(phi) + emf->Hb[ifreq][irec]*sin(phi);
	}//end for ifreq
      }//end for irec
	
    }else if(strcmp(emf->chrec[ichrec],"Hy")==0){
      for(irec=0; irec<acqui->nrec; irec++){
	phi = acqui->rec_azimuth[irec]-acqui->src_azimuth[isrc];
	for(ifreq=0; ifreq<emf->nfreq; ifreq++){
	  emf->dobs_fd[ichrec][ifreq][irec] = -emf->Ha[ifreq][irec]*sin(phi) + emf->Hb[ifreq][irec]*cos(phi);
	}//end for ifreq
      }//end for irec

    }//end if
      
  }

  //step 2: compute uncertainties associated with different components:
  //delta_Fx, delta_Fy depends on Fx and Fy (involving Fa and Fb)
  fwi->mse = 0;  
  for(ichrec=0; ichrec<emf->nchrec; ichrec++){
    for(ifreq=0; ifreq<emf->nfreq; ifreq++){
      for(irec=0; irec<acqui->nrec; irec++){

	tmp1 = 0.;
	phi = acqui->rec_azimuth[irec]-acqui->src_azimuth[isrc];
	if(strcmp(emf->chrec[ichrec],"Ex")==0){
	  d1 = cabs(emf->Ea[ifreq][irec])*cos(phi);
	  d2 = cabs(emf->Eb[ifreq][irec])*sin(phi);
	  tmp1 = emf->amp_perc*sqrt(d1*d1+d2*d2);
	}else	  if(strcmp(emf->chrec[ichrec],"Ey")==0){
	  d1 = cabs(emf->Ea[ifreq][irec])*sin(phi);
	  d2 = cabs(emf->Eb[ifreq][irec])*cos(phi);
	  tmp1 = emf->amp_perc*sqrt(d1*d1+d2*d2);
	}else	  if(strcmp(emf->chrec[ichrec],"Hx")==0){
	  d1 = cabs(emf->Ha[ifreq][irec])*cos(phi);
	  d2 = cabs(emf->Hb[ifreq][irec])*sin(phi);
	  tmp1 = emf->amp_perc*sqrt(d1*d1+d2*d2);
	}else	  if(strcmp(emf->chrec[ichrec],"Hy")==0){
	  d1 = cabs(emf->Ha[ifreq][irec])*sin(phi);
	  d2 = cabs(emf->Hb[ifreq][irec])*cos(phi);
	  tmp1 = emf->amp_perc*sqrt(d1*d1+d2*d2);
	}
	  
	tmp2 = 0.;
	for(ic=0; ic<emf->nchrec; ic++){
	  //reset tmp2 to non zero value if the other components exist, otherwise tmp2=0
	  if(strcmp(emf->chrec[ichrec],"Ex")==0 && strcmp(emf->chrec[ic],"Ey")==0){
	    tmp2 = emf->delta_phi*cabs(emf->dobs_fd[ic][ifreq][irec]);
	  }
	  if(strcmp(emf->chrec[ichrec],"Ey")==0 && strcmp(emf->chrec[ic],"Ex")==0){
	    tmp2 = emf->delta_phi*cabs(emf->dobs_fd[ic][ifreq][irec]);
	  }
	  if(strcmp(emf->chrec[ichrec],"Hx")==0 && strcmp(emf->chrec[ic],"Hy")==0){
	    tmp2 = emf->delta_phi*cabs(emf->dobs_fd[ic][ifreq][irec]);
	  }
	  if(strcmp(emf->chrec[ichrec],"Hy")==0 && strcmp(emf->chrec[ic],"Hx")==0){
	    tmp2 = emf->delta_phi*cabs(emf->dobs_fd[ic][ifreq][irec]);
	  }
	}
	/* if(strcmp(emf->chrec[ichrec],"Ex")==0||strcmp(emf->chrec[ichrec],"Ey")==0) */
	/*   emf->delta_emf[ichrec][ifreq][irec] = sqrt(tmp1*tmp1 + tmp2*tmp2 + emf->noisefloorE*emf->noisefloorE); */
	/* if(strcmp(emf->chrec[ichrec],"Hx")==0||strcmp(emf->chrec[ichrec],"Hy")==0) */
	/*   emf->delta_emf[ichrec][ifreq][irec] = sqrt(tmp1*tmp1 + tmp2*tmp2 + emf->noisefloorH*emf->noisefloorH); */
	
	if(strcmp(emf->chrec[ichrec],"Ex")==0||strcmp(emf->chrec[ichrec],"Ey")==0)
	  emf->delta_emf[ichrec][ifreq][irec] = sqrt(tmp1*tmp1 + tmp2*tmp2);
	if(strcmp(emf->chrec[ichrec],"Hx")==0||strcmp(emf->chrec[ichrec],"Hy")==0)
	  emf->delta_emf[ichrec][ifreq][irec] = sqrt(tmp1*tmp1 + tmp2*tmp2);

	if(emf->addnoise && emf->offset_ok[irec]){
	  if(strcmp(emf->chrec[ichrec],"Ex")==0||strcmp(emf->chrec[ichrec],"Ey")==0)
	    emf->delta_emf[ichrec][ifreq][irec] += emf->noisefloorE;
	  if(strcmp(emf->chrec[ichrec],"Hx")==0||strcmp(emf->chrec[ichrec],"Hy")==0)
	    emf->delta_emf[ichrec][ifreq][irec] += emf->noisefloorH;

	  /* frannor() compute Gaussian white noise (GWN) following N(0,1) distribution  */
	  cgwn = (frannor() + I*frannor())/sqrt(2.);//complex-valued GWN
	  emf->dobs_fd[ichrec][ifreq][irec] += emf->delta_emf[ichrec][ifreq][irec]*cgwn;
	    
	  fwi->mse += conj(cgwn)*cgwn;//accumulate mse value
	}
      }//end for irec
    }//end for ifreq      
      
  }//end for ichrec

  MPI_Allreduce(&fwi->mse, &tmp, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);/*sum over all process*/
  fwi->mse = tmp;
  if(emf->verb) printf("addnoise=%d, target rmse=%e\n", emf->addnoise, sqrt(fwi->mse/fwi->ndp));
  
  /* char fname[sizeof("uncertainty_0000.txt")]; */
  /* sprintf(fname, "uncertainty_%04d.txt", acqui->shot_idx[iproc]); */
  /* write_misfit(acqui, emf, fname, emf->delta_emf); */
  /* sprintf(fname, "obs_%04d.txt", acqui->shot_idx[iproc]); */
  /* write_data(acqui, emf, fname, emf->dobs_fd);/\* print out noise corrupted data *\/ */

}
