% mexGPUall. For these to complete succesfully, you need to configure the
% Matlab GPU library first (see README files for platform-specific
% information)
    mexcuda -largeArrayDims mexconv.cu    
%     mexcuda -largeArrayDims mexTemplate.cu

% If you get uninterpretable errors like "An unexpected error occurred during CUDA execution", and if you are using Pascal GPUs 
% (GTX 10xx series), it might be necessary to upgrade to Matlab 2017a / CUDA 8.0. 


