/* Pengliang Yang, Rune Mittet Nov. 20, 2021
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cstd.h"

const double PI = 3.1415926535897932384626;


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

/* nz_top, nz_middle, nz_bottom 3 parts included */
void create_canonical_reservoir()
{
  int ix, iy, iz, i;
  float xs1, xs2;
  float dx, dy, dz;
  float xmin, ymin, zmin;
  float xmax, ymax, zmax;
  int nx, ny, nz;
  int nxleft, nyleft, nzleft, nxright, nyright, nzright;
  float *x1nu, *x2nu, *x3nu;
  float *x1left, *x2left, *x3left, *x1right, *x2right, *x3right;
  FILE *fp;
  float r, xi, yi, zi, eta, tmp;
  float zseabed, zb;
  
  xmin = -9000;
  xmax = 9000;
  ymin = -9000;
  ymax = 9000;
  zmin = 0;
  zmax = 3500;
  
  xs1 = 0;
  xs2 = 0;
  zseabed = 1000;
  
  nx = 91;
  ny = 91;
  nz = 91;
  dx = 200;
  dy = 200;
  dz = 25;

  
  nxleft = (xs1-xmin)/(xmax-xmin)*(nx-1)+1;
  nyleft = (xs2-ymin)/(ymax-ymin)*(ny-1)+1;
  nxright = nx - nxleft + 1;
  nyright = ny - nyleft + 1;
  nzleft = zseabed/dz +1;
  nzright = nz - nzleft + 1;

  x1nu = alloc1float(nx);
  x2nu = alloc1float(ny);
  x3nu = alloc1float(nz);

  x1left = alloc1float(nxleft);
  x2left = alloc1float(nyleft);
  x1right = alloc1float(nxright);
  x2right = alloc1float(nyright);
  x3left = alloc1float(nzleft);
  x3right = alloc1float(nzright);

  r = create_nugrid(nxleft-1, xs1-xmin, dx, x1left);
  printf("r1left=%g\n", r);
  r = create_nugrid(nyleft-1, xs2-ymin, dy, x2left);  
  printf("r2left=%g\n", r);
  r = create_nugrid(nxright-1, xmax-xs1, dx, x1right);
  printf("r1right=%g\n", r);
  r = create_nugrid(nyright-1, ymax-xs2, dy, x2right);
  printf("r2right=%g\n", r);
  r = create_nugrid(nzleft-1, zseabed-zmin, dz, x3left);
  printf("r3left=%g\n", r);
  r = create_nugrid(nzright-1, zmax-zseabed, dz, x3right);
  printf("r3right=%g\n", r);
  
  for(ix=0; ix<nxleft; ix++)    x1nu[nxleft-1-ix] = xs1 - x1left[ix];
  for(ix=0; ix<nxright; ix++)   x1nu[nxleft-1+ix] = xs1 + x1right[ix];
  for(iy=0; iy<nyleft; iy++)    x2nu[nyleft-1-iy] = xs2 - x2left[iy];
  for(iy=0; iy<nyright; iy++)   x2nu[nyleft-1+iy] = xs2 + x2right[iy];
  for(iz=0; iz<nzleft; iz++)    x3nu[nzleft-1-iz] = zseabed - x3left[iz];
  for(iz=0; iz<nzright; iz++)   x3nu[nzleft-1+iz] = zseabed + x3right[iz];

  x1nu[0] = xmin;
  x1nu[nx-1] = xmax;
  x2nu[0] = ymin;
  x2nu[ny-1] = ymax;
  x3nu[0] = zmin;
  x3nu[nz-1] = zmax;

  for(iz=0; iz<nz; iz++) printf("z=%g\n", x3nu[iz]);

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
  free1float(x3left);
  free1float(x1right);
  free1float(x2right);
  free1float(x3right);

  /*--------------------------------------*/
  /* step 2: create seafloor horizon      */
  /*--------------------------------------*/
  int **ibathy = alloc2int(nx, ny);
  for(iy=0; iy<ny; iy++){
    yi = iy*dy + ymin;
    for(ix=0; ix<nx; ix++){
      xi = ix*dx + xmin;

      tmp = zseabed - 100.*cos(2.*PI*xi/(xmax-xmin));
      ibathy[iy][ix] = find_index(nz, x3nu, tmp);//+len to keep effective medium below bathy
    }
  }
  fp = fopen("ibathy", "wb");
  fwrite(&ibathy[0][0], nx*ny*sizeof(int), 1, fp);
  fclose(fp);

  fp = fopen("topo.txt", "w");
  for(ix=0; ix<nx; ix++){
    xi = ix*dx + xmin;

    tmp = zseabed - 100.*cos(2.*PI*xi/(xmax-xmin));
    fprintf(fp, "%e \t %e\n", xi, tmp);
  }
  
  fclose(fp);
  
  /*--------------------------------------------*/
  /* step3: asign values to rho_v and rho_h     */
  /*--------------------------------------------*/
  float ***rho = alloc3float(nx, ny, nz);

  float rho1 = 0.3;
  float rho2 = 1.5;
  float rho3 = 2.5;
  float resistor1 = 10;
  float resistor2 = 100;

  for(iz=0; iz<nz; iz++){
    zi = x3nu[iz];
    for(iy=0; iy<ny; iy++){
      yi = iy*dy + ymin;
      for(ix=0; ix<nx; ix++){
	xi = ix*dx + xmin;

	zb = zseabed - 100.*cos(2.*PI*xi/(xmax-xmin));//seafloor bathymetry
	if(zi<zb) rho[iz][iy][ix] = rho1;
	else      rho[iz][iy][ix] = rho2 + (rho3-rho2)*(zi-zb)/(zmax-zb);
	//if(zi>zsalt) rho[iz][iy][ix] = rho3;
	
	//rhoistor 1
	if(zi>=1200 && zi<=1250 && fabs(xi)<500 && yi>-500 && yi<1000) rho[iz][iy][ix] = resistor1;
	//rhoistor 2
	tmp = 3000;//disk rhoervoir model, radius
	xi -= -500;
	if(zi>=2200 && zi<=2350 && xi*xi + yi*yi < tmp*tmp)  rho[iz][iy][ix] = resistor2;
      }

    }
  }
  fp = fopen("rho", "wb");
  fwrite(&rho[0][0][0], nx*ny*nz*sizeof(float), 1, fp);
  fclose(fp);

  for(iz=0; iz<nz; iz++){
    zi = x3nu[iz];
    for(iy=0; iy<ny; iy++){
      yi = iy*dy + ymin;
      for(ix=0; ix<nx; ix++){
	xi = ix*dx + xmin;

	zb = zseabed - 100.*cos(2.*PI*xi/(xmax-xmin));//seafloor bathymetry
	if(zi<zb)  rho[iz][iy][ix] = rho1;
	else       rho[iz][iy][ix] = rho2;// + (rho3-rho2)*(zi-zb)/(zmax-zb);
      }

    }
  }
  fp = fopen("rho_init", "wb");
  fwrite(&rho[0][0][0], nx*ny*nz*sizeof(float), 1, fp);
  fclose(fp);

  free3float(rho);
  free1float(x1nu);
  free1float(x2nu);
  free1float(x3nu);

}



int main(int argc, char *argv[])
{
  create_canonical_reservoir();

}
