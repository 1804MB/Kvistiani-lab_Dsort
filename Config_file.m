			
ops.datatype            = 'openEphys';  % binary ('dat', 'bin') or 'openEphys'
ops.root                = 'C:\Users\au551108\Documents\light\2017-08-31_09-36-30blue'; % 'openEphys' only: where raw files are
ops.fbinary             = strcat(ops.root,'\example.dat'); % will be created for 'openEphys'		
ops.chanMap             = strcat(ops.root,'\chanMapT.mat'); % make this file using createChannelMapFile.m		

ops.nt0                 = 60;           % number of sampling point for the waveform		
ops.nt0min              = 20;           % position of the peak of the template
ops.fs                  = 30000;        % sampling rate	
ops.NchanTOT            = 32;           % total number of channels 
ops.Nchan               = 32;           % number of active channels 
ops.chan_per_group      = 4;            %Number of channels per group
ops.Nb_group            = ceil(ops.Nchan/ops.chan_per_group);%number of groups
ops.Nfilt               = 4;           % number of clusters to use 
% options for channel whitening		
ops.whitening           = 'noSpikes'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)		
ops.whiteningRange      = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)	
ops.scaleproc = 200;
		
% other options for controlling the model and optimization		
ops.Nrank               = 3;     % matrix rank of spike template model (3)	
ops.variance            = 0.9;   %variance to capture by the svd on average
ops.maxFR               = 20000; % maximum number of spikes to extract per batch (20000)		
ops.fshigh              = 300;   % Frequency for high pass filtering	
ops.slow                = 5000;  % Frequency of low pass filter
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection		
ops.NT                  = 32*1024+ ops.ntbuff;% this is the batch size (try decreasing if out of memory) 		
% for GPU should be multiple of 32 + ntbuff				
ops.criteria_NCC         = 0.9;          % criteria to detect spikes from the normalized cross correlation [number between -1 and 1]	
ops.max_itera            = 0.05;         % Max number of iteration 
%Criteria that control the post-processor
ops.Threshold = 0.85;                        % Control the criteria to merge clusters
ops.Chan_criteria = 0.2;                     %Criteria used to select the best channels

% options for initializing spikes from data	
%Import external template
ops.template = 'P'  ;                                      %'Y' (imported), 'G' (artificial), 'P' for tetrode peak to peak, 'N' (Kilosort initialization)
ops.template_file = strcat(ops.root,'\S_template.mat')  ;   %Name of the file containaing the template matrix, only used of previous option is 'Y'
ops.freeze = 'Y';                                          %update the template (N) or not (Y)
%********************Initilialization criteria
%if ops.template = 'P'  ;  
ops.T_crit          = 4;        %value for spike detection thresholding using option ops.Template = P;
ops.Ns_min          = 100;      %minimum nmber of spikes to detect
%ops.template = 'N'  ;  , from original Kilosort
ops.Th              = [4 6];    % threshold for detecting spikes on template-filtered data ([6 12 12])		
ops.lam             = [5 20 20];   % large means amplitudes are forced around the mean ([10 30 30])	
ops.spkTh           = -4;      % spike threshold in standard deviations (4)		
ops.loc_range       = [3  4];  % ranges to detect peaks; plus/minus in time and channel ([3 1])		
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])		
ops.maskMaxChannels = 4;       % how many channels to mask up/down ([5])		
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)		
ops.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)	
dd                  = load('PCspikes2.mat'); % you might want to recompute this from your own data		
ops.wPCA            = dd.Wi(:,1:7);   % PCs 		
s = zeros(ops.nt0,ops.Nrank);
for i=1:ops.Nrank
    s(:,i)=spline(1:length(ops.wPCA),ops.wPCA(:,i),1:ops.nt0);
end
ops.wPCA =s;
ops.ForceMaxRAMforDat   = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.

%Bit per microvolt conversion
ops.bitmVolt = 0.1950;


