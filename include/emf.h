#ifndef emf_h
#define emf_h

typedef struct {
  int mode; //mode=0, modeling; mode=1, FWI; mode=2, FWI gradient only
  int verb;/* verbose display */
  int reciprocity;
  
  float omega0;
  int n1, n2, n3, nb, ne, nbe, n123;
  int n1pad, n2pad, n3pad, n123pad;
  float d1, d2, d3, dt;

  int nchsrc, nchrec;//number of active src/rec channels for E and H
  char **chsrc, **chrec;

  int nfreq;//number of frequencies
  float *freqs, *omegas;//a list of frequencies
  float f0;//dominant frequency for the wavelet
  float freqmax;
  float Glim;
  float waterdepth;

  float ***rho11, ***rho22, ***rho33;//2 normal + 1 transverse resistivities
  float rhomin, rhomax; //mimum and maximum resistivity
  float vmin, vmax;//minimum and maximum velocities of the EM wave
  
  int airwave;//1=modeling with sea surface; 0=no sea surface
  int rd1, rd2, rd3;/* half length/radius of the interpolation operator */
  
  int n1fft, n2fft;
  float _Complex ***sH1kxky, ***sH2kxky;
  float ***sE12kxky;

  int ncorner;
  int gnopt;//1=Gauss-Newton inversion

  int nugrid;//nonuniform grid in z-axis
  float dx3_start, dx3_end;
  float *x3nu;//gridding over x3
  float *x3n, *x3s; //padded x3 for nonstaggered and staggered coordinates
  float **v3, **v3s;//FD coefficients for computing 1st derivative
  float **u3, **u3s;//weights for computing field
  

  int nt; //number of time steps in total
  float *stf;//source time function
  int **ibathy;//index of bathymetry along z direction

  float ***inveps11, ***inveps22, ***inveps33;//fictitous domain dielectric permittivity
  float *a1, *a2, *a3, *b1, *b2, *b3;

  float ***E1, ***E2, ***E3, ***H1, ***H2, ***H3;
  float ***curlE1, ***curlE2, ***curlE3, ***curlH1, ***curlH2, ***curlH3;
  float ***memD2H3, ***memD3H2, ***memD3H1, ***memD1H3, ***memD1H2, ***memD2H1;
  float ***memD2E3, ***memD3E2, ***memD3E1, ***memD1E3, ***memD1E2, ***memD2E1;

  /* Gauss-Newton inversion part */
  float ***gn_E1, ***gn_E2, ***gn_E3, ***gn_H1, ***gn_H2, ***gn_H3;
  float ***gn_curlE1, ***gn_curlE2, ***gn_curlE3, ***gn_curlH1, ***gn_curlH2, ***gn_curlH3;
  float ***gn_memD2H3, ***gn_memD3H2, ***gn_memD3H1, ***gn_memD1H3, ***gn_memD1H2, ***gn_memD2H1;
  float ***gn_memD2E3, ***gn_memD3E2, ***gn_memD3E1, ***gn_memD1E3, ***gn_memD1E2, ***gn_memD2E1;
  float ***v_rho11, ***v_rho22, ***v_rho33;

  
  /* computing box bound */
  int *i1min_fwd, *i1max_fwd;
  int *i2min_fwd, *i2max_fwd;
  int *i3min_fwd, *i3max_fwd;
  int *i1min_adj, *i1max_adj;
  int *i2min_adj, *i2max_adj;
  int *i3min_adj, *i3max_adj;

  float ***rho;//effective resistivity used for inversion
  float ***mask;//inversion mask to determine where should be updated
  float ***delta_emf;//uncertainty of electromagnetic fields(EMF)
  float ***obj;//significant misfit
  float _Complex ***dcal_fd, ***dobs_fd, ***dres_fd;
  float ***dres_td;


  /* fields in size of the model size * nfreq */
  float _Complex ****fwd_E1,****fwd_E2, ****fwd_E3;
  float _Complex ****fwd_H1,****fwd_H2, ****fwd_H3;
  float _Complex ****adj_E1,****adj_E2, ****adj_E3;
  float _Complex ****adj_H1,****adj_H2, ****adj_H3;
  float _Complex **expfactor;


  float _Complex **Ea, **Eb, **Ha, **Hb;

  int addnoise;
  int invmask;
  float offset_start;
  int *offset_ok;
  float *offset_weight;
  float noisefloorE, noisefloorH;
  float amp_perc;/* 0.03=3% deviation from the true value */
  float delta_phi;
} emf_t; /* type of electromagnetic field (emf)  */

#endif
