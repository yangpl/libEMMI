/* full waveform inversion of the CSEM data
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
#include "lbfgs.h"
#include "fwi.h"

void log_exp_scaling(int n, float *x, float *y, int flag);
void fg_fwi_init(acqui_t *acqui_, emf_t *emf_, fwi_t *fwi_);
void fg_fwi_close();
float fg_fwi(float *x, float *g);
void Hv_fwi(float *x, float *v, float *Hv);


void do_fwi(acqui_t *acqui, emf_t *emf)
{
  int i, j, k, i1, i2, i3, offset;
  float fcost, *g;
  lbfgs_t *opt; /* pointer for lbfgs_t parameters */
  fwi_t *fwi;
  FILE *fp;

  fwi = (fwi_t*)malloc(sizeof(fwi_t));
  opt = (lbfgs_t*)malloc(sizeof(lbfgs_t));
  
  /*-------------------------------------------------------------------------*/
  if(!getparint("niter", &opt->niter)) opt->niter = 30;/* maximum number of iterations */
  if(!getparint("nls", &opt->nls)) opt->nls = 5;/* maximum number of line searches */
  if(!getparfloat("tol", &opt->tol)) opt->tol = 1e-6;/* convergence tolerance */
  if(!getparint("npair", &opt->npair)) opt->npair = 5; /* l-BFGS memory length */
  if(!getparfloat("c1", &opt->c1)) opt->c1 = 1e-4; /* Nocedal value for Wolfe condition */
  if(!getparfloat("c2", &opt->c2)) opt->c2 = 0.9;  /* Nocedal value for Wolfe condition */
  if(!getparfloat("alpha", &opt->alpha)) opt->alpha = 1.;  /* initial step length */
  if(!getparint("bound", &opt->bound)) opt->bound = 0;/* 1 = bound on, 0 = off */
  if(!getparint("method", &opt->method)) opt->method = 0;//0=lBFGS; 1=Guass-Newton
  if(!getparint("ncg", &opt->ncg)) opt->ncg = 5;//Guass-Newton inversion
  opt->verb = emf->verb; /* output message only on process 0, others remain silent */
  emf->gnopt = opt->method?1:0;
  
  /*-------------------------------------------------------------------------*/
  if(!getparint("npar", &fwi->npar)) fwi->npar = 1;/* number of inversion parameters */
  fwi->idxpar = alloc1int(fwi->npar);
  getparint("idxpar", fwi->idxpar);/* indices of the inversion parameters */
  if(opt->bound) {
    if(!(j=countparval("minpar"))) 
      err("Need lower bound of parameters minpar= vector");
    if(j!=fwi->npar) err("npar=%d, must have length[idxpar]=npar", fwi->npar);
    if(!(j=countparval("maxpar"))) 
      err("Need upper bound of parameters maxpar= vector");
    if(j!=fwi->npar) err("npar=%d, must have length[idxpar]=npar", fwi->npar);

    fwi->minpar=alloc1float(fwi->npar);
    fwi->maxpar=alloc1float(fwi->npar);
    getparfloat("minpar", fwi->minpar);
    getparfloat("maxpar", fwi->maxpar);
  }
  if(opt->verb){
    printf("---------------------------------------\n");
    printf("**** method=%d (0=l-BFGS; 1=Gauss-Newton)\n", opt->method);
    printf("**** niter=%d\n", opt->niter);
    printf("**** nls=%d\n", opt->nls);
    printf("**** bound=%d\n", opt->bound);
    printf("**** ncg=%d\n", opt->ncg);
    if(opt->bound) printf("**** [lb,ub]=[%g, %g]\n", fwi->minpar[0], fwi->maxpar[0]);
  }  

  fwi->niter = opt->niter;
  fwi->iter = 0; /* before start set iter=0 */
  fwi->n = emf->n123*fwi->npar;/* lenth of unknown vector */
  fg_fwi_init(acqui, emf, fwi);


  /*---------------------------------------------------------------------------*/
  opt->x = alloc1float(fwi->n);/* model vector */
  opt->g = alloc1float(fwi->n);/* gradient vector */
  opt->d = alloc1float(fwi->n);/* descent direction */
  opt->sk = alloc2float(fwi->n, opt->npair);
  opt->yk = alloc2float(fwi->n, opt->npair);
  offset = 0;
  for(k=0; k<fwi->npar; k++){
    if(fwi->idxpar[k]==1){//Rv with index 1
      log_exp_scaling(emf->n123, &emf->rho_v[0][0][0], &opt->x[offset], 1);
      offset += emf->n123;
    }else if(fwi->idxpar[k]==2){// Rh with index 2
      log_exp_scaling(emf->n123, &emf->rho_h[0][0][0], &opt->x[offset], 1);
      offset += emf->n123;
    }
  }
  memcpy(fwi->xref, opt->x, fwi->n*sizeof(float));
  g = alloc1float(fwi->n);// a copy of opt->g

  if(opt->bound){
    opt->xmin = alloc1float(fwi->n);
    opt->xmax = alloc1float(fwi->n);
    for(j = 0; j<fwi->npar; j++){
      for(i3=0; i3<emf->n3; i3++){
	for(i2=0; i2<emf->n2; i2++){
	  for(i1=0; i1<emf->n1; i1++){
	    i = i1 + emf->n1*(i2 + emf->n2*i3);
	    if(i3>emf->ibathy[i2][i1]){//below bathymetry, specify xmin,xmax
	      /* note that the unknown x[:] =log(m[:])*/
	      opt->xmin[i+j*emf->n123] = log(fwi->minpar[j]);
	      opt->xmax[i+j*emf->n123] = log(fwi->maxpar[j]);
	    }else{//above the bathymetry, use input parameter, same xmin,xmax
	      //this allows min value set to be formation 1 Ohm-m directly
	      opt->xmin[i+j*emf->n123] = opt->x[i+j*emf->n123];
	      opt->xmax[i+j*emf->n123] = opt->x[i+j*emf->n123];
	    }
	  }
	}
      }
    }
  }


  /*---------------------------------------------------------------------------*/  
  fcost = fg_fwi(opt->x, opt->g);
  /*
  if(emf->mode==3){
    opt->f0 = fcost;
    fwi->gamma1 = 0.;
    fwi->gamma2 = 0.;
    float *r = alloc1float(fwi->n);

    for(i3=0; i3<emf->n3; i3++){
      for(i2=0; i2<emf->n2; i2++){
	for(i1=0; i1<emf->n1; i1++){
	  i = i1 + emf->n1*(i2 + emf->n2*i3);
	  opt->d[i] = opt->g[i]/emf->rho[i3][i2][i1];//convert g_log(rho) to g_rho
	  r[i] = 2.*rand()/(float)RAND_MAX -1.;//random uniform distribution between (-1,1)
	  if(i3<=emf->ibathy[i2][i1]) {
	    r[i] = 0.;
	    opt->d[i] = 0.;
	  }

	}
      }
    }
    float tmp = dotprod(fwi->n, opt->d, r);
    float eps = 1;
    for(int k=0; k<10; k++){
      eps *= 0.1;
      for(i3=0; i3<emf->n3; i3++){
	for(i2=0; i2<emf->n2; i2++){
	  for(i1=0; i1<emf->n1; i1++){
	    i = i1 + emf->n1*(i2 + emf->n2*i3);
	    opt->x[i] = log(exp(fwi->xref[i]) + eps*r[i]);
	  }
	}
      }      
      opt->fk = fg_fwi(opt->x, opt->g);
      if(emf->verb) printf("J(m+eps*r)-J(m)-eps<J',r>=%e ratio=%g\n",
			   opt->fk - opt->f0 - eps*tmp, (opt->fk - opt->f0)/(eps*tmp));
    }
    
    free1float(r);
  }
  */
  
  if(emf->mode==1){/* FWI */
    opt->f0 = fcost;
    opt->fk = fcost;
    opt->igrad = 0;
    opt->kpair = 0;
    opt->ils = 0;
    if(opt->verb){
      opt->gk_norm = l2norm(fwi->n, opt->g);
      fp=fopen("iterate.txt", "w");
      fprintf(fp, "==========================================================\n");
      fprintf(fp, "l-BFGS memory length: %d\n", opt->npair);
      fprintf(fp, "Maximum number of iterations: %d\n", opt->niter);
      fprintf(fp, "Convergence tolerance: %3.2e\n", opt->tol);
      fprintf(fp, "maximum number of line search: %d\n", opt->nls);
      fprintf(fp, "initial step length: alpha=%g\n", opt->alpha);
      fprintf(fp, "==========================================================\n");
      fprintf(fp, "iter    fk       fk/f0      ||gk||    alpha    nls   ngrad\n");
      fclose(fp);
    }
    /* l-BFGS/Newton-CG optimization */
    for(opt->iter=0; opt->iter<opt->niter; opt->iter++){
      fwi->iter = opt->iter;
      if(opt->verb){
	printf("==========================================================\n");
	printf("# iter=%d  fk/f0=%g\n", opt->iter, opt->fk/opt->f0);
	opt->gk_norm=l2norm(fwi->n, opt->g);

	fp=fopen("iterate.txt", "a");
	fprintf(fp, "%3d   %3.2e  %3.2e   %3.2e  %3.2e  %3d  %4d\n", 
		opt->iter, opt->fk, opt->fk/opt->f0, opt->gk_norm, opt->alpha, opt->ils, opt->igrad);
	fclose(fp);

	if(opt->iter==0) fp=fopen("rmse_misfit.txt", "w");
	else             fp=fopen("rmse_misfit.txt", "a");
	fprintf(fp, "%d \t %g\n", fwi->iter, sqrt(2.*fwi->fcost/fwi->ndp));
	fclose(fp);
      }

      if(opt->method==1){
  	cg_solve(fwi->n, opt->x, opt->g, opt->d, Hv_fwi, opt);/* solve Hv = -g */
      }else{
	memcpy(g, opt->g, fwi->n*sizeof(float));
  	if(opt->iter==0){/* first iteration, no stored gradient */
	  flipsign(fwi->n, opt->g, opt->d);/* descent direction=-gradient */
	}else{
	  lbfgs_update(fwi->n, opt->x, opt->g, opt->sk, opt->yk, opt);
	  lbfgs_descent(fwi->n, opt->g, opt->d, opt->sk, opt->yk, opt);
	} 
	lbfgs_save(fwi->n, opt->x, opt->g, opt->sk, opt->yk, opt);
      }

      line_search(fwi->n, opt->x, opt->g, opt->d, fg_fwi, opt);
      if(opt->ls_fail){
	opt->kpair =0;//flush out stored information, restart lbfgs 
	memcpy(opt->g, g, fwi->n*sizeof(float));//copy g back to opt->g
	flipsign(fwi->n, opt->g, opt->d);/* descent direction=-gradient */
	lbfgs_save(fwi->n, opt->x, opt->g, opt->sk, opt->yk, opt);
	line_search(fwi->n, opt->x, opt->g, opt->d, fg_fwi, opt);
      }

      /* not break, then line search succeeds or descent direction accepted */
      if(opt->verb) {/* print out inverted physical parameter at each iteration */
	fp=fopen("param_final_Rv", "wb");
	fwrite(&emf->rho_v[0][0][0], emf->n123*sizeof(float), 1, fp);
	fclose(fp);

	fp=fopen("param_final_Rh", "wb");
	fwrite(&emf->rho_h[0][0][0], emf->n123*sizeof(float), 1, fp);
	fclose(fp);
      }
    } 
    if(opt->verb && opt->iter==opt->niter) {
      fp=fopen("iterate.txt", "a");
      fprintf(fp, "==>Maximum iteration number reached!\n");
      fclose(fp);
    }
  }


  fg_fwi_close();
  free1float(opt->x);
  free1float(opt->g);
  free1float(opt->d);
  free2float(opt->sk);
  free2float(opt->yk);
  if(opt->bound){
    free1float(opt->xmin);
    free1float(opt->xmax);
    free1float(fwi->minpar);
    free1float(fwi->maxpar);
  }
  free1int(fwi->idxpar);
  free(opt);
  free(fwi);
}

