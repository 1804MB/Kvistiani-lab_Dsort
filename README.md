# Kvistiani-lab_Dsort
Dsort spike sorter
Installation

Dsort requires CUDA in order to run. This is because the normalized cross correlation routine is coded in CUDA (no cpu version is provided with Dsort).
  
What to do to run Dsort?

1. Install matlab R2018a 
2. Install cuda 9.0 https://developer.nvidia.com/cuda-90-download-archive. 
3. Install visual studio 2013 https://visualstudio.microsoft.com/vs/older-downloads/
4. Download Dsort
5. In Matlab go to Dsort folder then subfolder CUDA, and run the mexGPUall.m, 
this will compile the CUDA part of Dsort. For more information about CUDA and MATLAB 
check https://se.mathworks.com/help/distcomp/run-mex-functions-containing-cuda-code.html

Use Dsort

1. Open the config_file.m adapt the code to your settings, the settings  are described in the comments of the config.file
2. Run the Master file

Output
The results are saved in  a mat file called rez.mat. The saved variables are

ops            contains all the parameter provided by user from the config file
xc			   x coordinate of the probe
yc             y coordiante of the probe
xcoord         idem
ycoord         idem
connected      active channels (1) inactive channels(0)
Wrot           Whitenning matrix
temp           contain number of batch
st3            detected spikes using the NCC without post-processing
dWU            Template used for detection (Number of sampling point by Number of Channel by number of template matrix)
W              Temporal component (dWU = W*d*U, here W= W*d)
U              Spatial component
Weight         Weight used for normalization (Nrank matrix)
st             Detected spikes after post-processing
Chan           Best channels for each template (dWU)
Chan_score     Score of each template using the 1st component of matrix U
score          NCC results between template used for normalization
Merge_cluster  Summary of which template has been merged, how many spikes, best channels
Best_Channel_M Value of each detected spikes at its peak on all channels
template       Recalculated template
M_template     Merged template
M_std_template Standard deviation of each merged template
PC             Projection of each spikes onto its merged_template and best channels 
PETH           Auto-correlogram of each cluster

Post-processor
The user has the possibility to run a post-processor to vizualize and clean the cluster found by Dsort.
To do so, run the file post_master.

This displays for each cluster 4 figures:
1. Projection of each spikes onto its merged_template and best channels (determined via the matrix U)
2. The histogram of NCC score
3. The spike amplitude for the best channels
4. The average waveform with standard deviation + the PETH

The displayed prompt allow the user to choose different option:
- press 'Ok' to keep cluster
- enter '1' to delete cluster
- enter '2' to split the cluster, by default the cluster will be split in 2 using a guassian distribution on the best channel. 
The user can specify the number of cluster and the channel to use .
- enter '3' to clean cluster using the amplitude. The user needs to specify the amplitude threshold and the channel to use.
- enter '4' to clean the cluster using the normalized cross correlation score. The user needs to specify the NCC threshold to use.
 
 After dealing with all the cluster the results is saved in a different mat file called "rezD.mat". 
 The user can run the post-master as much as neccesary.
