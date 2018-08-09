%*************************************************************************%
%*************************************************************************%
%PreprocessData: Filter, whiten the raw data 
%and provide if needed the pc projectionsused to initialize the template
%Simplified and commented by Madeny belkhiri 29/01/18.
%*************************************************************************%
%*************************************************************************%
function [rez, DATA] = preprocessData(ops)


%Load the channel map information
if ~isempty(ops.chanMap)
    if ischar(ops.chanMap)
        load(ops.chanMap);
        try
            chanMapConn = chanMap(connected>1e-6);
            xc = xcoords(connected>1e-6);
            yc = ycoords(connected>1e-6);
        catch
            chanMapConn = 1+chanNums(connected>1e-6);
            xc = zeros(numel(chanMapConn), 1);
            yc = [1:1:numel(chanMapConn)]';
        end
        ops.Nchan    = getOr(ops, 'Nchan', sum(connected>1e-6));
        ops.NchanTOT = getOr(ops, 'NchanTOT', numel(connected));
        if exist('fs', 'var')
            ops.fs       = getOr(ops, 'fs', fs);
        end
    else
        chanMap = ops.chanMap;
        chanMapConn = ops.chanMap;
        xc = zeros(numel(chanMapConn), 1);
        yc = [1:1:numel(chanMapConn)]';
        connected = true(numel(chanMap), 1);      
        
        ops.Nchan    = numel(connected);
        ops.NchanTOT = numel(connected);
    end
%Create a channel map if not provided by user
else
    chanMap  = 1:ops.Nchan;
    connected = true(numel(chanMap), 1);
    
    chanMapConn = 1:ops.Nchan;    
    xc = zeros(numel(chanMapConn), 1);
    yc = [1:1:numel(chanMapConn)]';
end

%Check which channels are actived
if exist('kcoords', 'var')
    kcoords = kcoords(connected);
else
    kcoords = ones(ops.Nchan, 1);
end

%Number of active Channels
NchanTOT = ops.NchanTOT;
%Size of each batch
NT       = ops.NT ;
%Save in the rez file the config data
rez.ops         = ops;

%Coordinates on the probe
rez.xc = xc;
rez.yc = yc;
if exist('xcoords')
   rez.xcoords = xcoords;
   rez.ycoords = ycoords;
else
   rez.xcoords = xc;
   rez.ycoords = yc;
end
rez.connected   = connected;
rez.ops.chanMap = chanMap;
rez.ops.kcoords = kcoords; 

%Directory of the recording
d = dir(ops.fbinary);
%number of sampling point to read
ops.sampsToRead = floor(d.bytes/NchanTOT/2);

%allocate the memory based on the size of the recording
if ispc
    dmem         = memory;
    memfree      = dmem.MemAvailableAllArrays/8;
    memallocated = min(ops.ForceMaxRAMforDat, dmem.MemAvailableAllArrays) - memfree;
    memallocated = max(0, memallocated);
else
    memallocated = ops.ForceMaxRAMforDat;
end
nint16s      = memallocated/2;

%The definition of NTbuff is not clear, why 4*ops.ntbuff [MB?]
NTbuff      = NT + 4*ops.ntbuff;
Nbatch      = ceil(d.bytes/2/NchanTOT /(NT-ops.ntbuff)); %Number of batch
Nbatch_buff = floor(4/5 * nint16s/rez.ops.Nchan /(NT-ops.ntbuff)); % factor of 4/5 for storing PCs of spikes
Nbatch_buff = min(Nbatch_buff, Nbatch);

%Filtering coefficients butterworth of 3rd order
[b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.slow/ops.fs]*2, 'bandpass');

%Start routine to read the raw data file for filtering and covariance matrix calculation
fid = fopen(ops.fbinary, 'r');
ibatch = 0;
Nchan = rez.ops.Nchan;
%Covariance matrix
CC = gpuArray.zeros( Nchan,  Nchan, 'single');
nPairs = gpuArray.zeros( Nchan,  Nchan, 'single');

%Initialize the matrix DATA that will containt the filtered data
if ~exist('DATA', 'var')
    DATA = zeros(NT, rez.ops.Nchan, Nbatch_buff, 'int16');
end

isproc = zeros(Nbatch, 1);
%Filter the data batch per batch, 
%so while batch are not empty do
while 1

    ibatch = ibatch + 1;
	%Deleted ibatch + ops.nSkipCov because we want to calculate the whitening matrix for each batch
	%Don't see why skiping some batch would be usefull 
    
	
	%Offset the data point depending of the batch
    offset = max(0, 2*NchanTOT*( (NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff) );
    if ibatch==1
        ioffset = 0;
    else
        ioffset = ops.ntbuff;
    end
    fseek(fid, offset, 'bof');
	%read the data Number of channel * size of batch buffer
    buff = fread(fid, [NchanTOT NTbuff], '*int16');
    
    %If nothing to read stop
    if isempty(buff)
        break;
    end
	%Size of the current batch, will be different for the last batch
    nsampcurr = size(buff,2);
	
	%For the last batch we remap the matrix, sicne the dimension is smaller than NTbuff
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
    end
	
	%Load the batch data on the gpu
    dataRAW = gpuArray(buff);

    dataRAW = dataRAW';
    dataRAW = single(dataRAW);
	%Inly keep the data of the active channels
    dataRAW = dataRAW(:, chanMapConn);
    
	%Forward filtering using butterworth coefficient
    datr = filter(b1, a1, dataRAW);
    datr = flipud(datr);
	%Backward filtering using butterworth coefficient
    datr = filter(b1, a1, datr);
    datr = flipud(datr);
    
%Detect spikes in order to uncorrelate the signal to build a better whitening matrix
%Calculation of the covariance matrix
            smin      = my_min(datr, ops.loc_range, [1 2]);
            sd = std(datr, [], 1);
			%Simple peak detection threshold based on the standart deviation
			%ops.spkTh provided by user
            peaks     = single(datr<smin+1e-3 & bsxfun(@lt, datr, ops.spkTh * sd));
            blankout  = 1+my_min(-peaks, ops.long_range, [1 2]);
            smin      = datr .* blankout;
            CC        = CC + (smin' * smin)/NT; %covariance matrix, addition because we operate per batch
            nPairs    = nPairs + (blankout'*blankout)/NT;
    
        DATA(:,:,ibatch) = gather_try(int16( datr(ioffset + (1:NT),:))); %Save the filtered matrix datr (GPU) in DATA 

end
%CC = CC / ceil((Nbatch-1)/ops.nSkipCov); %Deletion of CC / ceil((Nbatch-1)/ops.nSkipCov) because we don't skip any batch 
CC = CC / ceil(Nbatch-1);% Divide by number of batch
nPairs = nPairs/ibatch;
CC = CC ./nPairs;%Why divide by nPairs? [MB?]

fclose(fid);
%close file



%Using SVD we get the whitening Matrix Wrot

%If only whiten some channels
if ops.whiteningRange<Inf
    ops.whiteningRange = min(ops.whiteningRange, Nchan);
    Wrot = whiteningLocal(gather_try(CC), yc, xc, ops.whiteningRange);
	
else %If whiten all channels
    [E, D] 	= svd(CC);
    D = diag(D);
    eps 	= 1e-6; %added for problem of singular matrix
    Wrot 	= E * diag(1./(D + eps).^.5) * E';
end


Wrot = Wrot *ops.scaleproc;
%Apply the whitening matrix to the filtered data and get the PC projection of the spikes
%used to initialize the template used in the next function FitTemplate.m
for ibatch = 1:Nbatch
	datr = single(gpuArray(DATA(:,:,ibatch)));
	%Whitening of the data
    datr    = datr * Wrot;
    %Save the filtered on whitened data 
    DATA(:,:,ibatch) = gather_try(datr);   
end

Wrot        = gather_try(Wrot);
%save the whitening matrix in the rez file
rez.Wrot    = Wrot;

rez.temp.Nbatch = Nbatch;
rez.temp.Nbatch_buff = Nbatch_buff;

