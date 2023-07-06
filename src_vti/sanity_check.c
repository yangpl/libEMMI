/* - Sanity check for stability condition, dispersion requirement
 * - determine an optimal temporal sampling dt and number of time steps nt
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"
#include "constants.h"
#include "mpi_info.h"


void sanity_check(emf_t *emf)
{
  int i, i1, i2, i3;
  float tmp,tmp1,tmp2,kappa,Rmax,cfl;
  float D1, D2, D3, s3, t3;
  float ***inveps11, ***inveps22, ***inveps33;

  
  inveps11 = alloc3float(emf->n1, emf->n2, emf->n3);
  inveps22 = alloc3float(emf->n1, emf->n2, emf->n3);
  inveps33 = alloc3float(emf->n1, emf->n2, emf->n3);


  /*----------------------------------------------------------------------------------*/
  /* step 1: find rhomax and rhomin, sigma=2*omega0*eps --> inveps=2*omega0*rho */
  /*----------------------------------------------------------------------------*/
  emf->rhomax = MAX( MAX(emf->rho11[0][0][0], emf->rho22[0][0][0]), emf->rho33[0][0][0]);
  emf->rhomin = MIN( MIN(emf->rho11[0][0][0], emf->rho22[0][0][0]), emf->rho33[0][0][0]);
  for(i3=0; i3<emf->n3; i3++){
    for(i2=0; i2<emf->n2; i2++){
      for(i1=0; i1<emf->n1; i1++){
	tmp1 = MAX( MAX(emf->rho11[i3][i2][i1], emf->rho22[i3][i2][i1]), emf->rho33[i3][i2][i1]);
	tmp2 = MIN( MIN(emf->rho11[i3][i2][i1], emf->rho22[i3][i2][i1]), emf->rho33[i3][i2][i1]);
	if(emf->airwave==1 && i3==0){//sigma_[11,22]=0.5*sigma_water, rho_interface=2*rho_water
	  tmp1 *= 2;
	  tmp2 *= 2;
	}
	if(tmp1>emf->rhomax) emf->rhomax = tmp1;
	if(tmp2<emf->rhomin) emf->rhomin = tmp2;
	
  	inveps11[i3][i2][i1] = 2.*emf->omega0*emf->rho11[i3][i2][i1];
  	inveps22[i3][i2][i1] = 2.*emf->omega0*emf->rho22[i3][i2][i1];
  	inveps33[i3][i2][i1] = 2.*emf->omega0*emf->rho33[i3][i2][i1];
      }
    }
  }


  /*----------------------------------------------------------------------------------*/
  /* Stage 2: find minimum and maximum velocity for stability conditon and dispersion */
  /*    emf->vmin: important for minimum number of points per wavelength              */
  /*    emf->vmax: important for CFL condition and fdtd computing box                 */
  /*----------------------------------------------------------------------------------*/
  if(emf->rd1==1)      D1 = 1.;
  else if(emf->rd1==2) D1 = (fabs(1.125) + fabs(-0.04167));
  else if(emf->rd1==3) D1 = (fabs(1.17188) + fabs(-0.06510) + fabs(0.00469));
  if(emf->rd2==1)      D2 = 1.;
  else if(emf->rd2==2) D2 = (fabs(1.125) + fabs(-0.04167));
  else if(emf->rd2==3) D2 = (fabs(1.17188) + fabs(-0.06510) + fabs(0.00469));
  if(emf->rd3==1)      D3 = 1.;
  else if(emf->rd3==2) D3 = (fabs(1.125) + fabs(-0.04167));
  else if(emf->rd3==3) D3 = (fabs(1.17188) + fabs(-0.06510) + fabs(0.00469));
  D1 *= 2;
  D2 *= 2;
  D3 *= 2;
  D1 /= emf->d1;
  D2 /= emf->d2;
  D3 /= emf->d3;
  
  kappa = 0;
  tmp1 = MAX( MAX(inveps11[0][0][0], inveps22[0][0][0]), inveps33[0][0][0]);
  tmp2 = MIN( MIN(inveps11[0][0][0], inveps22[0][0][0]), inveps33[0][0][0]);
  emf->vmax = sqrt(tmp1*invmu0);
  emf->vmin = sqrt(tmp2*invmu0);
  for(i3=0; i3<emf->n3; ++i3){
    int i3_ = i3+emf->nbe;
    if(emf->nugrid){
      s3 = 0;
      t3 = 0;
      for(i=0; i<2*emf->rd3; i++) {
	s3 += fabs(emf->v3[i3_][i]);
	t3 += fabs(emf->v3s[i3_][i]);
      }
      D3 = MAX(s3, t3);
    }

    for(i2=0; i2<emf->n2; ++i2){
      for(i1=0; i1<emf->n1; ++i1){

	tmp1 = MAX( MAX(inveps11[i3][i2][i1], inveps22[i3][i2][i1]), inveps33[i3][i2][i1]);
	tmp2 = MIN( MIN(inveps11[i3][i2][i1], inveps22[i3][i2][i1]), inveps33[i3][i2][i1]);
	if(emf->airwave==1 && i3==0){//sigma_[11,22]=0.5*sigma_water, rho_interface=2*rho_water
	  tmp1 *= 2;
	  tmp2 *= 2;
	}
	tmp1 = sqrt(tmp1*invmu0);
	tmp2 = sqrt(tmp2*invmu0);
	if(tmp1>emf->vmax) emf->vmax = tmp1;
	if(tmp2<emf->vmin) emf->vmin = tmp2;
	
	tmp = 0.5*sqrt(D1*D1 + D2*D2 + D3*D3);
	tmp *= emf->vmax;
	if(tmp>kappa) kappa = tmp;
      }
    }
  }

  /*------------------------------------------------------------------------*/
  /* Stage 3: determine the optimal dt and nt automatically                 */
  /*------------------------------------------------------------------------*/
  if(!getparfloat("dt", &emf->dt)) emf->dt = 0.98/kappa;
  /* temporal sampling, determine dt by stability condition if not provided */
  cfl = emf->dt*kappa;
  if(iproc==0) printf("cfl=%g\n", cfl); 
  if(cfl > 1.0) err("CFL condition not satisfied!");
  emf->freqmax = emf->vmin/(emf->Glim*MIN(MIN(emf->d1,emf->d2),emf->d3));

  if(!getparint("nt", &emf->nt)){
    Rmax = MAX((emf->n1-1)*emf->d1, (emf->n2-1)*emf->d2);
    emf->nt = (int)(2.*Rmax/(emf->vmin*emf->dt));
  }/* automatically determine nt using maximum offset if not provided */
  if(iproc==0){
    printf("[vmin, vmax]=[%g, %g] m/s\n", emf->vmin, emf->vmax);
    printf("[rhomin, rhomax]=[%g, %g] Ohm-m\n", emf->rhomin, emf->rhomax);
    printf("freqmax=%g Hz\n", emf->freqmax);
    printf("dt=%g s\n",  emf->dt);
    printf("nt=%d\n",  emf->nt);
  }

  free3float(inveps11);
  free3float(inveps22);
  free3float(inveps33);

}

