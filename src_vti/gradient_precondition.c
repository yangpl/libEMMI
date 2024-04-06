#include "cstd.h"
#include "acqui.h"
#include "emf.h"
#include "fwi.h"
#include "constants.h"

/*------------------------------------------------------*/
void gradient_precondition(acqui_t *acqui, emf_t *emf, fwi_t *fwi, float *g)
/*< precondition vector x of size n >*/
{
  int i1, i2, i3, i;
  float skin_depth, z, factor;

  skin_depth = sqrt(2*invmu0/emf->omegas[0]);//consider rho_ref=1 Ohm-m
  for(i3=0; i3<emf->n3; i3++){
    for(i2=0; i2<emf->n2; i2++){
      for(i1=0; i1<emf->n1; i1++){
	i = i1 + emf->n1*(i2 + emf->n2*i3);
	if(i3>emf->ibathy[i2][i1]){
	  z = emf->nugrid?(emf->x3n[i3] - acqui->x3min):(i3*emf->d3+acqui->x3min);
	  z -= acqui->src_x3[0];// z-zb

	  factor = exp(z/skin_depth);
	  g[i] *= factor;
	}else
	  g[i] = 0.;
      }
    }

  }
}
