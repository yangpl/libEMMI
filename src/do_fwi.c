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
void precondition(float *x);
float fg_fwi(float *x, float *g);
void Hv_fwi(float *x, float *v, float *Hv);


void do_fwi(acqui_t *acqui, emf_t *emf)
{
  int i, j, i1, i2, i3;
  float fcost;
  lbfgs_t *opt; /* pointer for lbfgs_t parameters */
  fwi_t *fwi;
  FILE *fp;

  fwi = (fwi_t*)malloc(sizeof(fwi_t));
  opt = (lbfgs_t*)malloc(sizeof(lbfgs_t));
  
  /*-------------------------------------------------------------------------*/
  if(!getparint("niter", &opt->niter)) opt->niter = 30;/* maximum number of iterations */
  if(!getparint("nls", &opt->nls)) opt->nls = 10;/* maximum number of line searches */
  if(!getparfloat("tol", &opt->tol)) opt->tol = 1e-8;/* convergence tolerance */
  if(!getparint("npair", &opt->npair)) opt->npair = 5; /* l-BFGS memory length */
  if(!getparfloat("c1", &opt->c1)) opt->c1 = 1e-4; /* Nocedal value for Wolfe condition */
  if(!getparfloat("c2", &opt->c2)) opt->c2 = 0.9;  /* Nocedal value for Wolfe condition */
  if(!getparfloat("alpha", &opt->alpha)) opt->alpha = 1.;  /* initial step length */
  if(!getparint("bound", &opt->bound)) opt->bound = 0;/* 1 = bound on, 0 = off */
  if(!getparint("preco", &opt->preco)) opt->preco = 1;/* 1 = precondition on, 0 = off */
  if(!getparint("method", &opt->method)) opt->method = 0;//0=lBFGS; 1=Guass-Newton
  if(!getparint("ncg", &opt->ncg)) opt->ncg = 5;//Guass-Newton inversion
  opt->verb = emf->verb; /* output message only on process 0, others remain silent */
  emf->gnopt = opt->method?1:0;
  
  /*-------------------------------------------------------------------------*/
  if(!getparint("npar", &fwi->npar)) fwi->npar = 1;/* number of inversion parameters */
  if(!(j=countparval("idxpar"))) err("Need lower bound of parameters idxpar= vector");
  if(j!=fwi->npar) err("must have length[idxpar]=%d", fwi->npar);
  fwi->idxpar = alloc1int(fwi->npar);
  getparint("idxpar", fwi->idxpar);/* indices of the inversion parameters */
  if(opt->bound) {
    if(!(j=countparval("minpar"))) err("Need lower bound of parameters minpar= vector");
    if(j!=fwi->npar) err("must have length[idxpar]=length[minpar]");
    if(!(j=countparval("maxpar"))) err("Need upper bound of parameters maxpar= vector");
    if(j!=fwi->npar) err("must have length[idxpar]=length[maxpar]");

    fwi->minpar=alloc1float(fwi->npar);
    fwi->maxpar=alloc1float(fwi->npar);
    getparfloat("minpar", fwi->minpar);
    getparfloat("maxpar", fwi->maxpar);
  }
  if(opt->verb){
    printf("---------------------------------------\n");
    printf("**** method=%d (0=l-BFGS; 1=Gauss-Newton)\n", opt->method);
    printf("**** niter=%d\n", opt->niter);
    printf("**** preco=%d\n", opt->preco);
    printf("**** nls=%d\n", opt->nls);
    printf("**** bound=%d\n", opt->bound);
    printf("**** ncg=%d\n", opt->ncg);
    if(opt->bound) printf("**** [lb,ub]=[%g, %g]\n", fwi->minpar[0], fwi->maxpar[0]);
  }  

  fwi->niter = opt->niter;
  fwi->iter = 0; /* before start set iter=0 */
  fwi->n = emf->n123*fwi->npar;/* lenth of unknown vector */
  fwi->preco = opt->preco;
  fg_fwi_init(acqui, emf, fwi);

  /*---------------------------------------------------------------------------*/
  opt->x = alloc1float(fwi->n);/* model vector */
  opt->g = alloc1float(fwi->n);/* gradient vector */
  if(opt->preco) opt->pg = alloc1float(fwi->n);
  opt->d = alloc1float(fwi->n);/* descent direction */
  opt->sk = alloc2float(fwi->n, opt->npair);
  opt->yk = alloc2float(fwi->n, opt->npair);
  if(emf->mode==3) memcpy(opt->x, &emf->rho[0][0][0], emf->n123*sizeof(float));
  else {
    log_exp_scaling(emf->n123, &emf->rho[0][0][0], opt->x, 1);
    memcpy(fwi->xref, opt->x, fwi->n*sizeof(float));
  }
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
	      opt->xmin[i+j*emf->n123] = log(emf->rho[i3][i2][i1]);
	      opt->xmax[i+j*emf->n123] = log(emf->rho[i3][i2][i1]);
	    }
	  }
	}
      }
    }
  }


  /*---------------------------------------------------------------------------*/  
  fcost = fg_fwi(opt->x, opt->g);
  if(emf->mode==3){
    opt->f0 = fcost;
    float *r = alloc1float(fwi->n);
    float eps = 0.01;
    for(i3=0; i3<emf->n3; i3++){
      for(i2=0; i2<emf->n2; i2++){
	for(i1=0; i1<emf->n1; i1++){
	  i = i1 + emf->n1*(i2 + emf->n2*i3);
	  r[i] = rand() / (float)RAND_MAX;
	  /* if(i3>emf->ibathy[i2][i1]){//below bathymetry, specify xmin,xmax */
	  /* }else{//above the bathymetry, use input parameter, same xmin,xmax */
	  /*   r[i] = 0.; */
	  /* } */
	  opt->x[i] += eps*r[i];
	  opt->d[i] = opt->g[i];
	}
      }
    }
    
    opt->fk = fg_fwi(opt->x, opt->g);
    float tmp = dotprod(fwi->n, opt->d, r);
    if(emf->verb) printf("(J(m+eps*r)-J(m))/eps=%e \t <J',r>=%e\n", (opt->fk-opt->f0)/eps, tmp);
    free1float(r);
  }

  
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
  	if(opt->iter==0){/* first iteration, no stored gradient */
	  if(opt->preco){
	    memcpy(opt->pg, opt->g, fwi->n*sizeof(float));
	    precondition(opt->pg);
	    flipsign(fwi->n, opt->pg, opt->d);/* descent direction=-preconditioned gradient */
	  }else
	    flipsign(fwi->n, opt->g, opt->d);/* descent direction=-gradient */
	}else{
	  opt->q = alloc1float(fwi->n);
	  opt->rho = alloc1float(opt->kpair);
	  opt->alp = alloc1float(opt->kpair);

	  lbfgs_update(fwi->n, opt->x, opt->g, opt->sk, opt->yk, opt);

	  /* 1st loop of two-loop recursion */
	  opt->loop1=lbfgs_descent1(fwi->n, opt->g, opt->q, opt->rho, opt->alp, opt->sk, opt->yk, opt);
	  /* precondition */
	  if(opt->preco) precondition(opt->q);
	  /* 2nd loop of two-loop recursion if 1st loop was done */
	  if(opt->loop1) lbfgs_descent2(fwi->n, opt->g, opt->q, opt->rho, opt->alp, opt->sk, opt->yk, opt);
	  flipsign(fwi->n, opt->q, opt->d); /* descent direction d=-q where q=H^{-1}g */
	  /* lbfgs_descent(n, opt->g, opt->d, opt->sk, opt->yk, opt); */

	  free(opt->q);
	  free(opt->alp);
	  free(opt->rho);
	} 
	lbfgs_save(fwi->n, opt->x, opt->g, opt->sk, opt->yk, opt);
      }

      line_search(fwi->n, opt->x, opt->g, opt->d, fg_fwi, opt);
      if(opt->ls_fail){
	if(opt->verb) {
	  printf("==>Line search failed, ils=%d\n", opt->ils);
	  fp=fopen("iterate.txt", "a");
	  fprintf(fp, "==>Line search failed!\n");
	  fclose(fp);
	}
	break;
      }

      /* not break, then line search succeeds or descent direction accepted */
      if(opt->verb) {/* print out inverted physical parameter at each iteration */
	log_exp_scaling(emf->n123, opt->x, &emf->rho[0][0][0], 2);
	fp=fopen("param_final", "wb");
	fwrite(&emf->rho[0][0][0], emf->n123*sizeof(float), 1, fp);
	fclose(fp);

	/* if(fwi->iter==0) fp=fopen("param_iter", "wb"); */
	/* else             fp=fopen("param_iter", "ab"); */
	/* fwrite(&emf->rho[0][0][0], emf->n123*sizeof(float), 1, fp); */
	/* fclose(fp); */
      }
      if(emf->addnoise && 2.*fwi->fcost_dat <= fwi->mse){//reach target MSE misfit 
	if(opt->verb){
	  fp=fopen("iterate.txt", "a");
	  fprintf(fp, "==>Convergence reached!\n");
	  fclose(fp);
	}
	break;
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
  if(opt->preco) free1float(opt->pg);
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

