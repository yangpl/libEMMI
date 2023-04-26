/* Pengliang Yang, Nov. 20, 2021 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cstd.h"

const double PI = 3.1415926535897932384626;

double Six(double x)
{
  double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
  double x1, x2, num, den, fx, gx, si;

  if(x==0) return 0;

  x1 = fabs(x);
  x2 = x*x;

  if(x1<2.3){
    t0 = -8.99432372444805e-11;
    t1 = 4.81542491189951e-8 + x2*t0;
    t2 = -9.98260398342126e-6 + x2*t1;
    t3 = 9.86950559158117e-4 + x2*t2;
    t4 = -4.16045660913051e-2 + x2*t3;
    num = 1.0 + x2*t4;

    t0 =1.60360248325839e-12;
    t1 = 1.11222895369314e-9 + x2*t0;
    t2 = 4.07021596174597e-7 + x2*t1;
    t3 = 9.53388627275856e-5 + x2*t2;
    t4 = 1.39509894642504e-2 + x2*t3;
    den = 1.0 + x2*t4;

    si = x1*num/den;
  }else{
    t0 = 2.23396200314283e-4;
    t1 = 9.41205789559801e-3 + x1*t0;
    t2 = 1.92842523449289e-1 + x1*t1;
    t3 = 2.28624888100692e0 + x1*t2;
    t4 = 1.64949384718444e1 + x1*t3;
    t5 = 7.03579911255664e1 + x1*t4;
    t6 = 1.59787502038397e2 + x1*t5;
    t7 = 1.42357662863517e2 + x1*t6;
    t8 = 3.41402258238831e1 + x1*t7;
    num = 8.76577957656212e-2 + x1*t8;

    t0 = 2.23396200314282e-4;
    t1 = 9.41205789560655e-3 + x1*t0;
    t2 = 1.93289315843730e-1 + x1*t1;
    t3 = 2.30507299852978e0 + x1*t2;
    t4 = 1.68761553399714e1 + x1*t3;
    t5 = 7.47422708097091e1 + x1*t4;
    t6 = 1.89060299847209e2 + x1*t5;
    t7 = 2.43359383546479e2 + x1*t6;
    t8 = 1.35358811016980e2 + x1*t7;
    t9 = 2.31300105233788e1 + x1*t8;
    den = x1*t9;
    
    fx = num/den;

    /*---------------------------------------------------*/
    t0 = 8.30826965804889e-5;
    t1 = 3.89884706241243e-3 + x1*t0;
    t2 = 8.93684959082754e-2 + x1*t1;
    t3 = 1.19330218679743e0 + x1*t2;
    t4 = 9.73185708298441e0 + x1*t3;
    t5 = 4.64218699919365e1 + x1*t4;
    t6 = 1.09328851692030e2 + x1*t5;
    t7 = 5.20840820481407e1 + x1*t6;
    t8 = 2.25386141512083e0 + x1*t7;
    num = -8.32280818835294e-2 + x1*t8;

    t0 = 8.30826965804714e-5;
    t1 = 3.89884706247557e-3 + x1*t0;
    t2 = 8.98669920492969e-2 + x1*t1;
    t3 = 1.21669527870366e0 + x1*t2;
    t4 = 1.02610878228560e1 + x1*t3;
    t5 = 5.32542903990403e1 + x1*t4;
    t6 = 1.60523538558852e2 + x1*t5;
    t7 = 2.45539153620874e2 + x1*t6;
    t8 = 1.47667803988686e2 + x1*t7;
    t9 = 2.51189885193833e1 + x1*t8;
    den = x2*t9;

    gx = num/den;

    si = 0.5*PI - fx*cos(x1) - gx*sin(x1);
  }

  if(x>0) return si;
  else return -si;
}

float heaviside(float eta){ return 0.5 + Six(eta*PI)/PI; }


/* cubic spline interpolation using slopes */
void spline_interp(float *x,  float *y,  int n,  float *dy,  float *ddy, float *s)
{
  int j;
  float h0, h1, alpha, beta;

  dy[0] = -0.5;
  h0 = x[1]-x[0];
  s[0] = 3.0*(y[1]-y[0])/(2.0*h0)-ddy[0]*h0/4.0;
  for (j = 1;j<= n-2;j++){
    h1 = x[j+1]-x[j];
    alpha = h0/(h0+h1);
    beta = (1.0-alpha)*(y[j]-y[j-1])/h0;
    beta = 3.0*(beta+alpha*(y[j+1]-y[j])/h1);
    dy[j] = -alpha/(2.0+(1.0-alpha)*dy[j-1]);
    s[j] = (beta-(1.0-alpha)*s[j-1]);
    s[j] = s[j]/(2.0+(1.0-alpha)*dy[j-1]);
    h0 = h1;
  }
  dy[n-1] = (3.0*(y[n-1]-y[n-2])/h1+ddy[n-1]*h1/
	   2.0-s[n-2])/(2.0+dy[n-2]);
  for (j = n-2;j>= 0;j--)
    dy[j] = dy[j]*dy[j+1]+s[j];
  for (j = 0;j<= n-2;j++) s[j] = x[j+1]-x[j];
  for (j = 0;j<= n-2;j++){
    h1 = s[j]*s[j];
    ddy[j] = 6.0*(y[j+1]-y[j])/h1-2.0*(2.0*dy[j]+dy[j+1])/s[j];
  }
  h1 = s[n-2]*s[n-2];
  ddy[n-1] = 6.*(y[n-2]-y[n-1])/h1+2.*(2.*dy[n-1]+dy[n-2])/s[n-2];
}

void getval(float *x,  float *y,  int n,  float *dy,  float *ddy, float *s, float *t,  int m,  float *z,  float *dz,  float *ddz)
{
  int i, j;
  float h0, h1;

  for (j = 0;j<= m-1;j++){
    if (t[j]>= x[n-1]) i = n-2;
    else {
      i = 0;
      while (t[j]>x[i+1]) i = i+1;
    }
    h1 = (x[i+1]-t[j])/s[i];
    h0 = h1*h1;
    z[j] = (3.0*h0-2.0*h0*h1)*y[i];
    z[j] = z[j]+s[i]*(h0-h0*h1)*dy[i];
    dz[j] = 6.0*(h0-h1)*y[i]/s[i];
    dz[j] = dz[j]+(3.0*h0-2.0*h1)*dy[i];
    ddz[j] = (6.0-12.0*h1)*y[i]/(s[i]*s[i]);
    ddz[j] = ddz[j]+(2.0-6.0*h1)*dy[i]/s[i];
    h1 = (t[j]-x[i])/s[i];
    h0 = h1*h1;
    z[j] = z[j]+(3.0*h0-2.0*h0*h1)*y[i+1];
    z[j] = z[j]-s[i]*(h0-h0*h1)*dy[i+1];
    dz[j] = dz[j]-6.0*(h0-h1)*y[i+1]/s[i];
    dz[j] = dz[j]+(3.0*h0-2.0*h1)*dy[i+1];
    ddz[j] = ddz[j]+(6.0-12.0*h1)*y[i+1]/(s[i]*s[i]);
    ddz[j] = ddz[j]-(2.0-6.0*h1)*dy[i+1]/s[i];
  }

}

/* find the index k in x[] such that x[k]<= val <x[k+1] */
int find_index(int n, float *x, float val)
{
  int low=0;
  int high = n-1;
  int i = (low+high)/2;
  
  while(low<high){
    i=(low+high)/2;
    if(x[i]<=val) low = i;
    if(x[i]>val) high=i;
    if(low==high||low==high-1) break;
  }

  return low;
}

/* generate nonuniform grid using geometric progression */
/* determine the optimal ratio r by root finding using fixed point iteration */
float create_nugrid(int n, float len, float dx, float *x)
{
  int i;
  float r, rr;
  float eps = 1e-15;

  if(fabs(n*dx-len)<eps) {
    for(i=0; i<=n; i++) x[i] = i*dx;
    return 1;
  }
  
  //fixed point iterations to find the root-"r" of the equation len = dx *(r^n-1)/(r-1)
  // assume n intervals (n+1 points/nodes)
  r = 1.1;
  rr = 1;
  while(1){
    rr = pow(len*(r-1.)/dx + 1., 1./n);
    if(fabs(rr-r)<eps) break;
    r = rr;
  }

  x[0] = 0;
  for(i=1; i<=n; i++)    x[i] = (pow(r,i) - 1.)*dx/(r-1.);

  return r;
}

/*----------------------------------------------------------------*/
void create_tilted_seafloor()
{
  int nh;
  float *x, *y, *dy_dx, *d2y_dx2, *ss, *tt, *rr;
  float dhx, xi, zi;
  int iz, ix, iy, i;
  float dz_dx, d2z_dx2, xb, zb, xs, zs;
  float tmp, eta;
  float ***rho_h, ***rho_v;
  FILE *fp;

  int len = 4;
  int niter = 7;
  float rho_w = 0.3;
  float rho_f = 1;

  float xmin = -10000;
  float xmax = 10000;
  float ymin = -10000;
  float ymax = 10000;
  float zmin = 0;
  float zmax = 5000;
  
  int nx = 101;
  int ny = 101;
  int nz = 101;

  float dx = (xmax-xmin)/(nx-1);
  float dy = (ymax-ymin)/(ny-1);
  float dz = 20;

  float *x1nu = alloc1float(nx);
  float *x2nu = alloc1float(ny);
  float *x3nu = alloc1float(nz);

  /*--------------------------------------*/
  /* step 1: generate NU grid coordinates */
  /*--------------------------------------*/
  for(ix=0; ix<nx; ix++) x1nu[ix] = ix*dx + xmin;
  for(iy=0; iy<ny; iy++) x2nu[iy] = iy*dy + ymin;
  for(iz=0; iz<nz; iz++) x3nu[iz] = iz*dz + zmin;

  zmin = 1000;
  iz = zmin/dz;
  nh = nz-iz;
  float *x3right = alloc1float(nh);
  eta = create_nugrid(nh-1, zmax-zmin, dz, x3right);
  printf("stretching factor=%g\n", eta);
  for(i=0; i<nh; i++) x3nu[iz+i] = x3nu[iz]+x3right[i];
  free1float(x3right);

  zmin = 0;
  x1nu[0] = xmin;
  x1nu[nx-1] = xmax;
  x2nu[0] = ymin;
  x2nu[ny-1] = ymax;
  x3nu[0] = zmin;
  x3nu[nz-1] = zmax;
  
  
  fp = fopen("x1nu", "wb");
  fwrite(x1nu, nx*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("x2nu", "wb");
  fwrite(x2nu, ny*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("x3nu", "wb");
  fwrite(x3nu, nz*sizeof(float), 1, fp);
  fclose(fp);

  /*--------------------------------------*/
  /* step 2: create seafloor horizon      */
  /*--------------------------------------*/
  dhx = 10;
  nh = (xmax-xmin)/dhx + 1;//grid spacing of horizon

  x = alloc1float(nh);
  y = alloc1float(nh);
  dy_dx = alloc1float(nh);
  d2y_dx2 = alloc1float(nh);
  ss = alloc1float(nh);
  tt = alloc1float(nh);
  rr = alloc1float(nh);

  fp = fopen("topo.txt", "w");
  for(ix=0; ix<nh; ix++) {
    //xi = MIN(ix*dhx+xmin, xmax);
    xi = ix*dhx+xmin;
    x[ix] = xi;
    tmp = xmax-xmin;
    y[ix] = 600. + 200.*sin(2.*PI*xi*3/tmp);//3 periods of sine shape
    eta = 2.*PI*3/tmp;
    tt[ix] = 200.*cos(2.*PI*xi*3/tmp)*eta;//1st derivative, analytically known here
    rr[ix] = -200*sin(2.*PI*xi*3/tmp)*eta*eta;//2nd derivative, analytically known here
    fprintf(fp, "%e \t %e\n", x[ix], y[ix]);
  }
  fclose(fp);

  //in general, 1st + 2nd derivative can be approximated by spline, if only a set of (xi,yi) are given
  dy_dx[0] = tt[0];
  dy_dx[nh-1] = tt[nh-1];
  spline_interp(x, y, nh, dy_dx, d2y_dx2, ss);

  xi = 9*dx;
  xb = xi;
  zi = 400;
  for(i=0; i<niter; i++){
    getval(x, y, nh, dy_dx, d2y_dx2, ss, &xb, 1, &zb, &dz_dx, &d2z_dx2);
    //compute distance to (xi, zi+0.5*dz) to obtain rho_v (z shifted half cell)
    xs = (xb - xi);//normalize distance by grid spacing
    zs = (zb - (zi+0.5*dz));
    //d2z_dx2 = 0.;//set it to be 0 to avoid instability
    xb -= ( xs + zs*dz_dx )/(1. + dz_dx*dz_dx + zs*d2z_dx2);
    float dist = sqrt(xs*xs + zs*zs);
    printf("root finding: dist=%g xb=%g\n", dist, xb);
  }

  
  /*--------------------------------------------*/
  /* step3: asign values to rho_v and rho_h     */
  /*--------------------------------------------*/
  rho_h = alloc3float(nx, ny, nz);
  rho_v = alloc3float(nx, ny, nz);

  iy = 0;
  for(iz=0; iz<nz; iz++){
    zi = iz*dz;
    for(ix=0; ix<nx; ix++){
      xi = ix*dx + xmin;

      xb = xi;
      getval(x, y, nh, dy_dx, d2y_dx2, ss, &xb, 1, &zb, &dz_dx, &d2z_dx2);
      if(zi>zb){//below seafloor
	rho_h[iz][iy][ix] = rho_f;
	rho_v[iz][iy][ix] = rho_f;
      }else{//above seafloor
	rho_h[iz][iy][ix] = rho_w;
	rho_v[iz][iy][ix] = rho_w;
      }
      
      /*----- (ix,iz)=rho_v[ix, iz+0.5]-----*/
      xb = xi;
      for(i=0; i<niter; i++){
	getval(x, y, nh, dy_dx, d2y_dx2, ss, &xb, 1, &zb, &dz_dx, &d2z_dx2);
	//compute distance to (xi, zi+0.5*dz) to obtain rho_v (z shifted half cell)
	xs = xb - xi;//normalize distance by grid spacing
	zs = zb - (zi+0.5*dz);
	//d2z_dx2 = 0.;//set it to be 0 to avoid instability
	xb -= ( xs + zs*dz_dx )/(1. + dz_dx*dz_dx + zs*d2z_dx2);
      }
      xs /= dx;
      zs /= dz;
      eta = sqrt(xs*xs + zs*zs);    
      if(eta<len){//modify the medium within 10 stepsizes
	if(zs>0) eta = -eta;
	/* average over resistivity for vertical direction */
	tmp = rho_w + heaviside(eta)*(rho_f-rho_w);
	rho_v[iz][iy][ix] = tmp;
      }

      /*------- (ix,iz)=rho_h[ix+0.5, iz] ----*/
      xb = xi + 0.5*dx;
      for(i=0; i<niter; i++){
	getval(x, y, nh, dy_dx, d2y_dx2, ss, &xb, 1, &zb, &dz_dx, &d2z_dx2);
	//compute distance to (xi+0.5*dx, zi) to obtain rho_h (x shifted half cell)
	xs = xb - (xi+0.5*dx);
	zs = zb - zi;
	//d2z_dx2 = 0.;
	xb -= ( xs + zs*dz_dx )/(1. + dz_dx*dz_dx + zs*d2z_dx2);
      }
      xs /= dx;
      zs /= dz;
      eta = sqrt(xs*xs + zs*zs);    
      if(eta<len && fabs(zi-600)<=200){//modify the medium within 10 stepsizes
	if(zs>0) eta = -eta;
	/* average over conductivity for horizontal direction */
	tmp = 1./rho_w + heaviside(eta)*(1./rho_f - 1./rho_w);
	rho_h[iz][iy][ix] =1./tmp;
      }
	
    }//end for ix
  }//end for iz
  
  zmin = 1000;
  float rho1 = 1;
  float rho2 = 2;
  float rho3 = 10;
  float rho4 = 50;
  for(iy=0; iy<ny; iy++){
    for(ix=0; ix<nx; ix++){
      xi = ix*dx + xmin;
      for(iz=0; iz<nz; iz++){
	rho_v[iz][iy][ix] = rho_v[iz][0][ix];
	rho_h[iz][iy][ix] = rho_h[iz][0][ix];
	
	zi = x3nu[iz];
	if(zi>zmin){/* No homogenization below the seabed */
	  if(zi<1500){
	    rho_h[iz][iy][ix] = rho1;
	  }else	if(zi<2000){
	    rho_h[iz][iy][ix] = rho2;
	  }else{
	    rho_h[iz][iy][ix] = rho3;
	  }
	}	
	zi = 0.5*(x3nu[iz] + x3nu[MIN(iz+1,nz-1)]);
	if(zi>zmin){
	  if(zi<1500){
	    rho_v[iz][iy][ix] = rho1;
	  }else	if(zi<2000){
	    rho_v[iz][iy][ix] = rho2;
	  }else{
	    rho_v[iz][iy][ix] = rho3;
	  }
	}//end if

	if(fabs(xi)<2000 && zi>=2300 && zi<2500) rho_v[iz][iy][ix] = rho4;
	if(fabs(xi+0.5*dx)<2000 && zi>=2300 && zi<2500) rho_h[iz][iy][ix] = rho4;	  
      }
    }
  }//end for iy


  //check x-z section from 3D resistivity cube
  fp = fopen("check_rho11", "wb");
  iy = 0;
  for(ix=0; ix<nx; ix++)
    for(iz=0; iz<nz; iz++){
      tmp = rho_h[iz][iy][ix];
      fwrite(&tmp, sizeof(float), 1, fp);
    }
  fclose(fp);

  fp = fopen("check_rho33", "wb");
  iy = 0;
  for(ix=0; ix<nx; ix++)
    for(iz=0; iz<nz; iz++)
      fwrite(&rho_v[iz][iy][ix], sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("grid_rho11", "w");
  iy = 0;
  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){
      fprintf(fp, "%g %g %g\n", x3nu[iz], x1nu[ix], rho_h[iz][iy][ix]);
    }
    fprintf(fp, "\n");
  }
  for(iz=0; iz<nz; iz++){
    for(ix=0; ix<nx; ix++){
      fprintf(fp, "%g %g %g\n", x3nu[iz], x1nu[ix], rho_h[iz][iy][ix]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  //output the models
  fp = fopen("rho11", "wb");
  fwrite(&rho_h[0][0][0], nz*nx*ny*sizeof(float), 1, fp);
  fclose(fp);
  fp = fopen("rho22", "wb");
  fwrite(&rho_h[0][0][0], nz*nx*ny*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("rho33", "wb");
  fwrite(&rho_v[0][0][0], nz*nx*ny*sizeof(float), 1, fp);
  fclose(fp);

  float ox = -8000;
  float ex = 8000;
  nx = 81;
  dx = (ex - ox)/(nx-1);
  fp = fopen("acquisition.txt", "w");
  fprintf(fp, "%f \t %f \t %f \t %d \t %d %d\n", 0., 0., 300.0, 0, 0, 0);
  for(ix=0; ix<nx; ix++){
    xb = ix*dx + ox;
    getval(x, y, nh, dy_dx, d2y_dx2, ss, &xb, 1, &zb, &dz_dx, &d2z_dx2);

    fprintf(fp, "%f \t %f \t %f \t %d \t %d \t %d\n", xb, 0., zb-20, 0, 0, 1);
  }
  fclose(fp);


  fp = fopen("resistor_top.txt", "w");
  fprintf(fp, "%f \t %f\n", -2000., 2300.);
  fprintf(fp, "%f \t %f\n", 2000., 2300.);
  fclose(fp);

  fp = fopen("resistor_bottom.txt", "w");
  fprintf(fp, "%f \t %f\n", -2000., 2500.);
  fprintf(fp, "%f \t %f\n", 2000., 2500.);
  fclose(fp);

  free3float(rho_h);
  free3float(rho_v);

  free1float(x);
  free1float(y);
  free1float(x1nu);
  free1float(x2nu);
  free1float(x3nu);

}

/*----------------------------------------------------------------*/
void create_1d_model()
{
  int nh;
  int iz, ix, iy, i;
  float zi, zs, tmp, eta;
  float ***rho_h, ***rho_v;
  FILE *fp;

  float xmin = -10000;
  float xmax = 10000;
  float ymin = -10000;
  float ymax = 10000;
  float zmin = 0;
  float zmax = 5000;
  
  int len = 4;
  int nx = 101;
  int ny = 101;
  int nz = 101;

  float dx = (xmax-xmin)/(nx-1);
  float dy = (ymax-ymin)/(ny-1);
  float dz = 20;

  float *x1nu = alloc1float(nx);
  float *x2nu = alloc1float(ny);
  float *x3nu = alloc1float(nz);

  /*--------------------------------------*/
  /* step 1: generate NU grid coordinates */
  /*--------------------------------------*/
  for(ix=0; ix<nx; ix++) x1nu[ix] = ix*dx + xmin;
  for(iy=0; iy<ny; iy++) x2nu[iy] = iy*dy + ymin;
  for(iz=0; iz<nz; iz++) x3nu[iz] = iz*dz + zmin;


  iz = 30;
  zmin = iz*dz;
  nh = nz-iz;
  float *x3right = alloc1float(nh);
  eta = create_nugrid(nh-1, zmax-zmin, dz, x3right);
  printf("stretching factor=%g\n", eta);
  for(i=0; i<nh; i++) {
    x3nu[iz+i] = x3nu[iz]+x3right[i];
    //printf("z=%g\n", x3nu[iz+i]);
  }
  free1float(x3right);
  
  fp = fopen("x1nu", "wb");
  fwrite(x1nu, nx*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("x2nu", "wb");
  fwrite(x2nu, ny*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("x3nu", "wb");
  fwrite(x3nu, nz*sizeof(float), 1, fp);
  fclose(fp);

  /*--------------------------------------*/
  /* step 2: create seafloor horizon      */
  /*--------------------------------------*/
  /* nh = 5001; */
  /* float dhx = (xmax-xmin)/(nh-1);//grid spacing of horizon */
  /* float *hx = alloc1float(nh); */
  /* float *hz = alloc1float(nh); */
  /* for(ix=0; ix<nh; ix++){ */
  /*   hx[ix] = ix*dhx; */
  /*   hz[ix] = 330; */
  /* } */

  /*--------------------------------------------*/
  /* step3: asign values to rho_v and rho_h     */
  /*--------------------------------------------*/
  rho_h = alloc3float(nx, ny, nz);
  rho_v = alloc3float(nx, ny, nz);
  float rho1, rho2;
  iy = 0;
  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){

      zi = x3nu[iz];
      /*--------------------------------------------------*/
      if(zi<=325) tmp = 0.3;
      else if(zi<=1025) tmp = 1;
      else if(zi<=1525) tmp = 2;
      else tmp = 4;
      rho_h[iz][iy][ix] = tmp;
      rho_v[iz][iy][ix] = tmp;

      zi = 0.5*(x3nu[iz] + x3nu[MIN(iz+1,nz-1)]);
      /*--------------------------------------------------*/
      rho1 = 0.3;
      rho2 = 1;
      zs = zi-325;
      eta = zs/dz;
      if(fabs(eta)<len) {
	//tmp = heaviside(eta)*kaiser_window(eta, len);//Kaiser window to suppress ripples
	/* average over resistivity for vertical direction */
	rho_v[iz][iy][ix] = rho1 + heaviside(eta)*(rho2 -rho1);
	if(iy==0 && ix==0) printf("eta=%g\n", eta);// kaiser_window(eta, len));
      }

      /*--------------------------------------------------*/
      rho1 = 1;
      rho2 = 2;
      zs = zi-1025;
      //i = find_index(nz, x3nu, 1025);
      eta = zs/dz; //(x3nu[i+1]-x3nu[i]);
      if(fabs(eta)<len) {
	//tmp = heaviside(eta)*kaiser_window(eta, len);//Kaiser window to suppress ripples
      	/* average over resistivity for vertical direction */
      	rho_v[iz][iy][ix] = rho1 + heaviside(eta)*(rho2 -rho1);
      	if(iy==0 && ix==0) printf("eta=%g \n", eta);
      }

      /*--------------------------------------------------*/
      rho1 = 2;
      rho2 = 4;
      zs = zi-1525;
      //i = find_index(nz, x3nu, 1525);
      eta = zs/dz; //(x3nu[i+1]-x3nu[i]);
      if(fabs(eta)<len) {
	//tmp = heaviside(eta)*kaiser_window(eta, len);//Kaiser window to suppress ripples
      	/* average over resistivity for vertical direction */
      	rho_v[iz][iy][ix] = rho1 + heaviside(eta)*(rho2 -rho1);
      	if(iy==0 && ix==0) printf("eta=%g \n", eta);
      }

      if(iy==0 && ix==0) printf("z=%g rho=%g\n", zi, rho_v[iz][iy][ix]);
    }
  }

  for(iy=0; iy<ny; iy++){
    for(iz=0; iz<nz; iz++){
      for(ix=0; ix<nx; ix++){
	rho_v[iz][iy][ix] = rho_v[iz][0][ix];
	rho_h[iz][iy][ix] = rho_h[iz][0][ix];
      }
    }
  }//end for iy


  //output the models
  fp = fopen("rho11", "wb");
  fwrite(&rho_h[0][0][0], nz*nx*ny*sizeof(float), 1, fp);
  fclose(fp);
  fp = fopen("rho22", "wb");
  fwrite(&rho_h[0][0][0], nz*nx*ny*sizeof(float), 1, fp);
  fclose(fp);
  fp = fopen("rho33", "wb");
  fwrite(&rho_v[0][0][0], nz*nx*ny*sizeof(float), 1, fp);
  fclose(fp);

  free3float(rho_h);
  free3float(rho_v);

  free1float(x1nu);
  free1float(x2nu);
  free1float(x3nu);

}

void create_land_model()
{
  int ix, iy, iz, i;
  float xs1, xs2;
  float dx, dy, dz;
  float xmin, ymin, zmin;
  float xmax, ymax, zmax;
  int nx, ny, nz, nh;
  int nxleft, nyleft, nxright, nyright;
  float *x1nu, *x2nu, *x3nu;
  float *x1left, *x2left, *x3left, *x1right, *x2right, *x3right;
  FILE *fp;
  float r, xi, yi, zi, eta;


  int len = 3;
    
  xmin = -3000;
  xmax = 3000;
  ymin = -3000;
  ymax = 3000;
  zmin = 0;
  zmax = 1500;  

  xs1 = 0;
  xs2 = 0;
  
  nx = 101;
  ny = 101;
  nz = 61;

  dx = 60;
  dy = 60;
  dz = 25;

  nxleft = (xs1-xmin)/(xmax-xmin)*(nx-1)+1;
  nyleft = (xs2-ymin)/(ymax-ymin)*(ny-1)+1;
  nxright = nx - nxleft + 1;
  nyright = ny - nyleft + 1;
  
  x1nu = alloc1float(nx);
  x2nu = alloc1float(ny);
  x3nu = alloc1float(nz);

  x1left = alloc1float(nxleft);
  x2left = alloc1float(nyleft);
  x1right = alloc1float(nxright);
  x2right = alloc1float(nyright);

  r = create_nugrid(nxleft-1, xs1-xmin, dx, x1left);
  printf("r1left=%g\n", r);
  r = create_nugrid(nyleft-1, xs2-ymin, dy, x2left);  
  printf("r2left=%g\n", r);
  r = create_nugrid(nxright-1, xmax-xs1, dx, x1right);
  printf("r1right=%g\n", r);
  r = create_nugrid(nyright-1, ymax-xs2, dy, x2right);
  printf("r2right=%g\n", r);

  
  for(ix=0; ix<nxleft; ix++)    x1nu[nxleft-1-ix] = xs1 - x1left[ix];
  for(ix=0; ix<nxright; ix++)   x1nu[nxleft-1+ix] = xs1 + x1right[ix];
  for(iy=0; iy<nyleft; iy++)    x2nu[nyleft-1-iy] = xs2 - x2left[iy];
  for(iy=0; iy<nyright; iy++)   x2nu[nyleft-1+iy] = xs2 + x2right[iy];
  for(iz=0; iz<nz; iz++)        x3nu[iz] = iz*dz + zmin;


  x1nu[0] = xmin;
  x1nu[nx-1] = xmax;
  x2nu[0] = ymin;
  x2nu[ny-1] = ymax;
  x3nu[0] = zmin;
  x3nu[nz-1] = zmax;


  for(i=0; i<nz; i++) {
    printf("z=%g\n", x3nu[i]);
  }
    
  fp = fopen("x1nu", "wb");
  fwrite(x1nu, nx*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("x2nu", "wb");
  fwrite(x2nu, ny*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("x3nu", "wb");
  fwrite(x3nu, nz*sizeof(float), 1, fp);
  fclose(fp);

  free1float(x1left);
  free1float(x2left);
  free1float(x1right);
  free1float(x2right);

  /*--------------------------------------------*/
  /* step3: asign values to rho_v and rho_h     */
  /*--------------------------------------------*/
  float ***rho11 = alloc3float(nx, ny, nz);
  float ***rho22 = alloc3float(nx, ny, nz);
  float ***rho33 = alloc3float(nx, ny, nz);

  float xs, ys;
  float rho0 = 5;
  float rho1 = 1;
  float rho2 = 100;

  for(iz=0; iz<nz; iz++){
    for(iy=0; iy<ny; iy++){
      yi = iy*dy + ymin;
      for(ix=0; ix<nx; ix++){
	xi = ix*dx + xmin;
	
	rho11[iz][iy][ix] = rho0;
	rho22[iz][iy][ix] = rho0;
	rho33[iz][iy][ix] = rho0;


	/* zi = x3nu[iz]; */
	/* if(zi>=500 && zi<1000){ */
	/*   xs = xi+0.5*dx; */
	/*   ys = yi; */
	/*   if(xs>-1000 && xs<=0 && fabs(ys)<=1000) rho11[iz][iy][ix] = rho1; */
	/*   if(xs>0 && xs<=1000 && fabs(ys)<=1000) rho11[iz][iy][ix] = rho2; */
	  
	/*   xs = xi; */
	/*   ys = yi+0.5*dy; */
	/*   if(xs>-1000 && xs<=0 && fabs(ys)<=1000) rho22[iz][iy][ix] = rho1; */
	/*   if(xs>0 && xs<=1000 && fabs(ys)<=1000) rho22[iz][iy][ix] = rho2; */
	/* } */
	/* zi = 0.5*(x3nu[iz] +x3nu[MIN(iz+1,nz-1)]); */
	/* if(zi>=500 && zi<1000){ */
	/*   if(xi>-1000 && xi<=0 && fabs(yi)<=1000) rho33[iz][iy][ix] = rho1; */
	/*   if(xi>0 && xi<=1000 && fabs(yi)<=1000) rho33[iz][iy][ix] = rho2; */
	/* } */

      }
    }
  }//end for iy
  

  //output the models
  fp = fopen("rho11", "wb");
  fwrite(&rho11[0][0][0], nz*nx*ny*sizeof(float), 1, fp);
  fclose(fp);
  fp = fopen("rho22", "wb");
  fwrite(&rho22[0][0][0], nz*nx*ny*sizeof(float), 1, fp);
  fclose(fp);
  fp = fopen("rho33", "wb");
  fwrite(&rho33[0][0][0], nz*nx*ny*sizeof(float), 1, fp);
  fclose(fp);


  free3float(rho11);
  free3float(rho22);
  free3float(rho33);

  free1float(x1nu);
  free1float(x2nu);
  free1float(x3nu);
}

  
int main(int argc, char *argv[])
{
  int opt = 3;

  if(opt==1) create_1d_model();
  if(opt==2) create_tilted_seafloor();
  if(opt==3) create_land_model();
  
}
