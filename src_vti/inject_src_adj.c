/* inject electric and magnetic source for adjoint modelling
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "acqui.h"
#include "emf.h"
#include "interp.h"

void inject_electric_src_adj(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg, int it)
/*< inject a source time function into EM field >*/
{
  int ic,irec,isub,i1,i2,i3,ix1,ix2,ix3,i1_,i2_,i3_;
  float w1,w2,w3,s,t,d3;

  for(ic=0; ic<emf->nchrec; ++ic){
    for(irec = 0; irec<acqui->nrec; irec++) {
      s = emf->dres_td[it][ic][irec];
      s /= acqui->nsubrec;
      for(isub=0; isub<acqui->nsubrec; isub++){
	if(strcmp(emf->chrec[ic],"Ex") == 0){
	  /* staggered grid: E1[i1,i2,i3] = Ex[i1+0.5,i2,i3] */
	  ix1 = interp_sg->rec_i1[irec][isub];
	  ix2 = interp_rg->rec_i2[irec][isub];
	  ix3 = interp_rg->rec_i3[irec][isub];
	  d3 = (emf->nugrid)?(emf->x3s[ix3]-emf->x3s[ix3-1]):emf->d3;
	  t = s/(emf->d1*emf->d2*d3);/* normalized by volume */
	  for(i3=-emf->rd3+1; i3<=emf->rd3; i3++){
	    w3 = interp_rg->rec_w3[irec][isub][i3+emf->rd3-1];
	    i3_ = ix3+i3;
	    for(i2=-emf->rd2+1; i2<=emf->rd2; i2++){
	      w2 = interp_rg->rec_w2[irec][isub][i2+emf->rd2-1];
	      i2_ = ix2+i2;
	      for(i1=-emf->rd1+1; i1<=emf->rd1; i1++){
		w1 = interp_sg->rec_w1[irec][isub][i1+emf->rd1-1];
		i1_ = ix1+i1;
		/*interpolate over continuous field Ex */
		emf->curlH1[i3_][i2_][i1_] -= t*w1*w2*w3;
	      }
	    }
	  }

	}else if(strcmp(emf->chrec[ic],"Ey") == 0){
	  /* staggered grid: E2[i1,i2,i3] = Ey[i1,i2+0.5,i3] */
	  ix1 = interp_rg->rec_i1[irec][isub];
	  ix2 = interp_sg->rec_i2[irec][isub];
	  ix3 = interp_rg->rec_i3[irec][isub];
	  d3 = (emf->nugrid)?(emf->x3s[ix3]-emf->x3s[ix3-1]):emf->d3;
	  t = s/(emf->d1*emf->d2*d3);
	  for(i3=-emf->rd3+1; i3<=emf->rd3; i3++){
	    w3 = interp_rg->rec_w3[irec][isub][i3+emf->rd3-1];
	    i3_ = ix3+i3;
	    for(i2=-emf->rd2+1; i2<=emf->rd2; i2++){
	      w2 = interp_sg->rec_w2[irec][isub][i2+emf->rd2-1];
	      i2_ = ix2+i2;
	      for(i1=-emf->rd1+1; i1<=emf->rd1; i1++){
		w1 = interp_rg->rec_w1[irec][isub][i1+emf->rd1-1];
		i1_ = ix1+i1;
		/*interpolate over continuous field Ey */
		emf->curlH2[i3_][i2_][i1_] -= t*w1*w2*w3;
	      }
	    }
	  }
	}/* end if */

      }/* end for isub */
    }/* end for irec */
  }/* end for ic */

    
}

//==========================================================
void inject_magnetic_src_adj(acqui_t *acqui, emf_t *emf, interp_t *interp_rg, interp_t *interp_sg, int it)
/*< inject a source time function into EM field >*/
{
  
}

