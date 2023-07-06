/*--------------------------------------------------------------------------
 * ( f(x0))  (1  x0-x (x0-x)^2 ... (x0-x)^n) ( f(x)      )
 * ( f(x1)) =(1  x1-x (x1-x)^2 ... (x1-x)^n) ( f^1(x)    )
 * ( ...  )  (...                          ) ( ...       )
 * ( f(xn))  (1  xn-x (xn-x)^2 ... (xn-x)^n) ( f^n(x)/n! )
 *           -------------------------------
 *            V^T (V=Vandermonde matrix)
 * Given the vector f=(f(x0), f(x1), ..., f(xn))^T and Vandermonde matrix 
 * V(x0, x1, ..., xn), the solution of Vandermonde matrix inversion V^T a =f
 * gives:
 *  (a0)  (f (x)    )  (w0 w1 ... wn) ( f(x0) )
 *  (a1)= (f'(x)    ) =(            ) ( f(x1) )
 *  (.)     ...        (            ) ( ...   )
 *  (an)  (f^n(x)/n!)  (            ) ( f(xn) )
 *                     --------------
 *                      V^{-1}
 * Reference: Golub 1996 Matrix computation. section 4.6 Vandermonde system
 *--------------------------------------------------------------------------
 *
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 *-------------------------------------------------------------------------*/
void vandermonde(int n, float *x, float *a, float *f)
{
  int i, k;

  /* calculate Newton representation of the interpolating polynomial */
  for(i=0; i<=n; ++i) a[i] = f[i];

  for(k=0; k<n; k++){
    for(i=n; i>k; i--){
      a[i] = (a[i]-a[i-1])/(x[i]-x[i-k-1]);
    }
  }

  for(k=n-1; k>=0; k--){
    for(i=k; i<n; i++){
      a[i] -= a[i+1]*x[k];
    }
  }
}


