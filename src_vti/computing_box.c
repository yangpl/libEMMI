/* a computing box to avoid unnecessary computation in the domain
 *
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 */
#include "cstd.h"
#include "acqui.h"
#include "emf.h"


void computing_box_init(acqui_t *acqui, emf_t *emf, int adj)
{
  int n1e,n2e,n3e;/* number of grid point extended from the source/receiver */
  int isrc,irec,i1,i2,i3,it;
  float t;
  int shift = 3;/* 1 due to NINT + 2 due to 4th order FD */

  if(adj){
    emf->i1min_adj=alloc1int(emf->nt);
    emf->i1max_adj=alloc1int(emf->nt);
    emf->i2min_adj=alloc1int(emf->nt);
    emf->i2max_adj=alloc1int(emf->nt);
    emf->i3min_adj=alloc1int(emf->nt);
    emf->i3max_adj=alloc1int(emf->nt);
    i1 = (acqui->rec_x1[0]-acqui->x1min)/emf->d1+emf->nbe;
    i2 = (acqui->rec_x2[0]-acqui->x2min)/emf->d2+emf->nbe;
    i3 = (acqui->rec_x3[0]-acqui->x3min)/emf->d3+emf->nbe;
    emf->i1min_adj[emf->nt-1] = emf->i1max_adj[emf->nt-1]  =  i1;
    emf->i2min_adj[emf->nt-1] = emf->i2max_adj[emf->nt-1]  =  i2;
    emf->i3min_adj[emf->nt-1] = emf->i3max_adj[emf->nt-1]  =  i3;
    for(irec = 0; irec<acqui->nrec; irec++){
      i1  =  (acqui->rec_x1[irec]-acqui->x1min)/emf->d1+emf->nbe;
      i2  =  (acqui->rec_x2[irec]-acqui->x2min)/emf->d2+emf->nbe;
      i3  =  (acqui->rec_x3[irec]-acqui->x3min)/emf->d3+emf->nbe;
      emf->i1min_adj[emf->nt-1] = MIN(emf->i1min_adj[emf->nt-1], i1);
      emf->i1max_adj[emf->nt-1] = MAX(emf->i1max_adj[emf->nt-1], i1);
      emf->i2min_adj[emf->nt-1] = MIN(emf->i2min_adj[emf->nt-1], i2);
      emf->i2max_adj[emf->nt-1] = MAX(emf->i2max_adj[emf->nt-1], i2);
      emf->i3min_adj[emf->nt-1] = MIN(emf->i3min_adj[emf->nt-1], i3);
      emf->i3max_adj[emf->nt-1] = MAX(emf->i3max_adj[emf->nt-1], i3);
    }
    for(it = emf->nt-2; it>= 0; it--){
      t = (emf->nt-1-it)*emf->dt;
      n1e = NINT(emf->vmax*t/emf->d1)+shift;
      n2e = NINT(emf->vmax*t/emf->d2)+shift;
      n3e = NINT(emf->vmax*t/emf->d3)+shift;
      emf->i1min_adj[it] = MAX(emf->i1min_adj[emf->nt-1]-n1e, 0);
      emf->i1max_adj[it] = MIN(emf->i1max_adj[emf->nt-1]+n1e, emf->n1pad-1);
      emf->i2min_adj[it] = MAX(emf->i2min_adj[emf->nt-1]-n2e, 0);
      emf->i2max_adj[it] = MIN(emf->i2max_adj[emf->nt-1]+n2e, emf->n2pad-1);
      emf->i3min_adj[it] = MAX(emf->i3min_adj[emf->nt-1]-n3e, 0);
      emf->i3max_adj[it] = MIN(emf->i3max_adj[emf->nt-1]+n3e, emf->n3pad-1);
    }
  }else{
    emf->i1min_fwd = alloc1int(emf->nt);
    emf->i1max_fwd = alloc1int(emf->nt);
    emf->i2min_fwd = alloc1int(emf->nt);
    emf->i2max_fwd = alloc1int(emf->nt);
    emf->i3min_fwd = alloc1int(emf->nt);
    emf->i3max_fwd = alloc1int(emf->nt);
    emf->i1min_fwd[0] = NINT((acqui->src_x1[0]-acqui->x1min)/emf->d1)+emf->nbe;
    emf->i1max_fwd[0] = NINT((acqui->src_x1[0]-acqui->x1min)/emf->d1)+emf->nbe;
    emf->i2min_fwd[0] = NINT((acqui->src_x2[0]-acqui->x2min)/emf->d2)+emf->nbe;
    emf->i2max_fwd[0] = NINT((acqui->src_x2[0]-acqui->x2min)/emf->d2)+emf->nbe;
    emf->i3min_fwd[0] = NINT((acqui->src_x3[0]-acqui->x3min)/emf->d3)+emf->nbe;
    emf->i3max_fwd[0] = NINT((acqui->src_x3[0]-acqui->x3min)/emf->d3)+emf->nbe;
    for(isrc = 0; isrc<acqui->nsrc; ++isrc){
      emf->i1min_fwd[0] = MIN(emf->i1min_fwd[0],NINT((acqui->src_x1[isrc]-acqui->x1min)/emf->d1)+emf->nbe);
      emf->i1max_fwd[0] = MAX(emf->i1max_fwd[0],NINT((acqui->src_x1[isrc]-acqui->x1min)/emf->d1)+emf->nbe);
      emf->i2min_fwd[0] = MIN(emf->i2min_fwd[0],NINT((acqui->src_x2[isrc]-acqui->x2min)/emf->d2)+emf->nbe);
      emf->i2max_fwd[0] = MAX(emf->i2max_fwd[0],NINT((acqui->src_x2[isrc]-acqui->x2min)/emf->d2)+emf->nbe);
      emf->i3min_fwd[0] = MIN(emf->i3min_fwd[0],NINT((acqui->src_x3[isrc]-acqui->x3min)/emf->d3)+emf->nbe);
      emf->i3max_fwd[0] = MAX(emf->i3max_fwd[0],NINT((acqui->src_x3[isrc]-acqui->x3min)/emf->d3)+emf->nbe);
    }
    for(it = 1; it<emf->nt; ++it){
      t = it*emf->dt;
      n1e = NINT(emf->vmax*t/emf->d1)+shift;
      n2e = NINT(emf->vmax*t/emf->d2)+shift;
      n3e = NINT(emf->vmax*t/emf->d3)+shift;
      emf->i1min_fwd[it] = MAX(emf->i1min_fwd[0]-n1e,0);
      emf->i1max_fwd[it] = MIN(emf->i1max_fwd[0]+n1e,emf->n1pad-1);
      emf->i2min_fwd[it] = MAX(emf->i2min_fwd[0]-n2e,0);
      emf->i2max_fwd[it] = MIN(emf->i2max_fwd[0]+n2e,emf->n2pad-1);
      emf->i3min_fwd[it] = MAX(emf->i3min_fwd[0]-n3e,0);
      emf->i3max_fwd[it] = MIN(emf->i3max_fwd[0]+n3e,emf->n3pad-1);

    }
  }
}

void computing_box_close(emf_t *emf, int adj)
{
  if(adj){
    free(emf->i1min_adj);
    free(emf->i1max_adj);
    free(emf->i2min_adj);
    free(emf->i2max_adj);
    free(emf->i3min_adj);
    free(emf->i3max_adj);
  }else{
    free(emf->i1min_fwd);
    free(emf->i1max_fwd);
    free(emf->i2min_fwd);
    free(emf->i2max_fwd);
    free(emf->i3min_fwd);
    free(emf->i3max_fwd);
  }
}

