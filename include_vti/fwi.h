#ifndef fwi_h
#define fwi_h

typedef struct {
  int npar;
  int n; /* total number of unknowns in fwi */

  int preco;//1=precondition; 0=no precodition
  int precoopt;//1=Plessix; 2=Pseudo-Hessian
  int r1, r2, r3, repeat;
  
  int *idxpar;
  float *minpar;
  float *maxpar;
  int iter;
  int niter;

  float *xref; //a reference model to compute regularization term
  float gamma1;//Tikhonov regularization parameter
  float gamma2;//TV regularization parameter

  int firstgrad;//1=first grad; 0= grad at later iterations
  float alpha; //step length
  float fcost_dat, fcost_mod, fcost;
  float *fcost_list;
  float ***g_h, ***g_v, ***Hv;//gradient, pseudo-Hessian and Hessian vector product

  int *ndp_list;
  int ndp; // number of data points=nsrc*nfreq*nchrec*nrec
  float mse; // target mse misfit
} fwi_t;

#endif
