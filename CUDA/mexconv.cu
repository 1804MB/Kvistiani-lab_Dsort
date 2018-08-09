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

const int Nthreads = 1024;
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void	Conv1D(const int *Params, const float *data, const float *W, float *conv_sig){    
  __shared__ float  sW[81], sdata[(Nthreads+81)]; 
  float x,a,b;
  int tid, nt0, tid0, bid, i, nid, NT, Nfilt;

  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  Nfilt    	=   (int) Params[1];
  NT      	=   (int) Params[0];
  nt0       =   (int) Params[2];
  
  if(tid<nt0)
      sW[tid]= W[tid%nt0 + (bid + Nfilt * (tid/nt0))* nt0];
  __syncthreads();
  
  tid0 = 0;
  while (tid0<NT-Nthreads-nt0+1){
	  if (tid<nt0) sdata[tid%nt0 + (tid/nt0)*(Nthreads+nt0)] = 
			data[tid0 + tid%nt0+ NT*(bid + Nfilt*(tid/nt0))];
     
          sdata[tid + nt0+nid*(Nthreads+nt0)] = data[nt0+tid0 + tid+ NT*(bid +nid*Nfilt)];	  
	  __syncthreads();
      
	  
          x = 0.0f;
          a = 0.0f;
          b = 0.0f;
		  #pragma unroll 4
          for(i=0;i<nt0;i++){
              a    += sW[i + nid*nt0]*sW[i + nid*nt0];
              b    += sdata[i+tid + nid*(Nthreads+nt0)]*sdata[i+tid + nid*(Nthreads+nt0)];
              x    += sW[i + nid*nt0] * sdata[i+tid + nid*(Nthreads+nt0)];
           }
             x = x/(sqrt(a*b));
	  
      conv_sig[tid0  + tid + NT*bid]   = x;
      
      tid0+=Nthreads;
      __syncthreads();
  }
}
/*
 * Host code
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    /* Declare input variables*/
  int *Params, *d_Params;
  int blocksPerGrid, NT;
  int const threadsPerBlock = 1024;

  /* Initialize the MathWorks GPU API. */
  mxInitGPU();
 mxGPUArray *B;
  /* read Params and copy to GPU */
  Params        = (int*) mxGetData(prhs[0]);
  NT            = (int) Params[0];
  blocksPerGrid	= (int) Params[1];
  cudaMalloc(&d_Params,      sizeof(int)*mxGetNumberOfElements(prhs[0]));
  cudaMemcpy(d_Params,Params,sizeof(int)*mxGetNumberOfElements(prhs[0]),cudaMemcpyHostToDevice);
  /* collect input GPU variables*/
  mxGPUArray const  *W, *data; 
  const float      *d_W,*d_data;
  float *d_dout;
  const mwSize dimst1[] 	= {NT,blocksPerGrid}; 
  B = mxGPUCreateGPUArray(2,dimst1,mxSINGLE_CLASS, mxREAL,MX_GPU_DO_NOT_INITIALIZE);
  d_dout = (float *)(mxGPUGetData(B));

  W             = mxGPUCreateFromMxArray(prhs[2]);
  d_W        	= (float const *)(mxGPUGetDataReadOnly(W));
  data        	= mxGPUCreateFromMxArray(prhs[1]);
  d_data        = (float const *)(mxGPUGetDataReadOnly(data));
  
    
   Conv1D<<<blocksPerGrid,threadsPerBlock>>>(d_Params, d_data, d_W, d_dout);  
   plhs[0] = mxGPUCreateMxArrayOnGPU(B);

  
  cudaFree(d_Params); 
  mxGPUDestroyGPUArray(B);
  mxGPUDestroyGPUArray(data);
  mxGPUDestroyGPUArray(W); 
 
  
}