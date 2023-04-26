/* find the index k in x[] such that x[k]<= val <x[k+1] 
 *
 * Copyright (c) 2020 Pengliang Yang. All rights reserved.
 * Email: ypl.2100@gmail.com
 */
int find_index(int n, float *x, float val)
{
  /*assume x[] has been sorted ascendingly */
  /* int i; */
  
  /* for(i=0; i<n; i++){ */
  /*   if(val<x[i]) break; */
  /* } */

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
