function [dWU] = Template_building(DATA,ops)
[~, Nchan ,Nbatch] = size(DATA);   
Nfilt = ops.Nfilt;
% Deleted the option to choose from data or not, now only fromdata.
% thus, the initiale template will always be calculated from the data.
if(ops.template =='Y')
    dWU0   = load(ops.template_file);
    dWU    = dWU0.T;
elseif(ops.template=='G') % or it generates templates based on config file
    [WUinit] = init_template(ops,Nchan);
    dWU    = WUinit(:,:,1:Nfilt);
else
    
i0  = 0;
wPCA = ops.wPCA;


uproj = zeros(1e6,  size(wPCA,2) * Nchan, 'single');
%Apply the whitening matrix to the filtered data and get the PC projection of the spikes
%used to initialize the template used in the next function FitTemplate.m
for ibatch = 1:Nbatch
    dataRAW = single(DATA(:,:,ibatch))/ops.scaleproc;
	%Division by the magical number, why did we use it before? [MB?] 
    %If you wish to initial the template from the data, otherwise skipped

        % find isolated spikes
        [row, col, ~] = isolated_peaks(dataRAW, ops.loc_range, ops.long_range, ops.spkTh,ops);
        
        % find their PC projections
        uS = get_PCproj(dataRAW, row, col, wPCA, ops.maskMaxChannels);
        uS = permute(uS, [2 1 3]);
        uS = reshape(uS,numel(row), Nchan * size(wPCA,2));
       
        if i0+numel(row)>size(uproj,1)
            uproj(1e6 + size(uproj,1), 1) = 0;            
        end
        %Save the PC projection on the CPU 
        uproj(i0 + (1:numel(row)), :) = gather_try(uS);
        i0 = i0 + numel(row);
    %Save the filtered on whitened data     
end


uproj(i0+1:end, :) = []; 
    WUinit = optimizePeaks(ops,uproj,Nchan);%does a scaled kmeans to determine the  Nfilt clusters template
    dWU    = WUinit(:,:,1:Nfilt);
end

[dWU,top] = SVD_topchan(dWU,ops.Chan_criteria);
 dWU = alignWU(dWU,top, ops.nt0min); %make sure that the deep is set at nt0min vital for substracting the template from signal