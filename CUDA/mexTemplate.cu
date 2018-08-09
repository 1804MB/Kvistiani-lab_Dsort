/*
 * Example of how to use the mxGPUArray API in a MEX file.  This example shows
 * how to write a MEX function that takes a gpuArray input and returns a
 * gpuArray output, e.g. B=mexFunction(A).
 *
 * Copyright 2012 The MathWorks, Inc.
 */
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <math.h>
#include <stdint.h>
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cstdlib>
#include <algorithm>
#include <iostream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
__global__ void average_snips(const int *Params, const int *st, const int *id, const float *dataraw, float *WU, float const pm){
  int nt0, tidx, tidy, bid, ind, NT, Nchan,N_spike;
  float xsum = 0.0f; 
  NT            = (int) Params[0];
  N_spike       = (int) Params[2];
  nt0           = (int) Params[3];
  Nchan         = (int) Params[4];
  
  tidx 		= threadIdx.x;
  tidy 		= threadIdx.y;
  bid 		= blockIdx.x;
  
  for(ind=0; ind<N_spike;ind++)
      if (id[ind]==bid){
		  tidy 		= threadIdx.y;
		  while (tidy<Nchan){	
			xsum = dataraw[st[ind]+tidx + NT * tidy]; 
  			WU[tidx+tidy*nt0 + nt0*Nchan * bid] = pm*WU[tidx+tidy*nt0 + nt0*Nchan * bid] +  (1-pm)* xsum;
			tidy+=blockDim.y;
		  }
	  }
}
//////////////////////////////////////////////////////////////////////////////////////////


/*
 * Host code
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    /* Declare input variables*/
  int *Params, *d_Params;
  int blocksPerGrid, N_spike, nt0;
  float *dpm;
  float pm;
  /* Initialize the MathWorks GPU API. */
  mxInitGPU();

  /* read Params and copy to GPU */
  Params        = (int*) mxGetData(prhs[0]);
  dpm           = (float*)mxGetData(prhs[5]);
  blocksPerGrid	= (int) Params[1];
  N_spike       = (int) Params[2];
  nt0           = (int) Params[3];
  pm            = (float)dpm[0];
  cudaMalloc(&d_Params,      sizeof(int)*mxGetNumberOfElements(prhs[0]));
  cudaMemcpy(d_Params,Params,sizeof(int)*mxGetNumberOfElements(prhs[0]),cudaMemcpyHostToDevice);
  
  /* collect input GPU variables*/
  mxGPUArray const *dataraw;
  const float      *d_dataraw;
  dataraw       = mxGPUCreateFromMxArray(prhs[1]);
  d_dataraw     = (float const *)(mxGPUGetDataReadOnly(dataraw));

  float *d_dWU;
  mxGPUArray *dWU;
  dWU       = mxGPUCopyFromMxArray(prhs[4]);
  d_dWU     = (float *)(mxGPUGetData(dWU));
  
  
  /* allocate new GPU variables*/
  int *d_st,*d_id;
  int *t,*id;
  t = (int*)mxGetData(prhs[2]);
  id = (int*)mxGetData(prhs[3]);

  cudaMalloc(&d_st,    N_spike * sizeof(int));
  cudaMalloc(&d_id,    N_spike * sizeof(int));

  
  cudaMemcpy(d_st,t,N_spike *   sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(d_id,id,N_spike *   sizeof(int),cudaMemcpyHostToDevice);

  
  dim3 block(nt0, 1024/nt0);
  average_snips<<<blocksPerGrid,block>>>(d_Params, d_st, d_id, d_dataraw, d_dWU,pm);
 
  plhs[0] 	= mxGPUCreateMxArrayOnGPU(dWU);
  int   *x,*idt;
  const mwSize dimst[] 	= {N_spike ,1}; 
  plhs[1] = mxCreateNumericArray(2, dimst, mxINT32_CLASS, mxREAL);
  plhs[2] = mxCreateNumericArray(2, dimst, mxINT32_CLASS, mxREAL);
  x = (int*) mxGetData(plhs[1]);
  cudaMemcpy(x,   d_st, N_spike * sizeof(int), cudaMemcpyDeviceToHost);
  idt = (int*) mxGetData(plhs[2]);
  cudaMemcpy(idt,   d_id, N_spike * sizeof(int), cudaMemcpyDeviceToHost);

  cudaFree(dpm);
  cudaFree(d_st);
  cudaFree(d_id);
  cudaFree(d_Params);
  mxGPUDestroyGPUArray(dataraw);
  mxGPUDestroyGPUArray(dWU); 
  
}