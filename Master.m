%************************************************************************************%
%************************************************************************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Master file for D.sort%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************************************************************************************%
%************************************************************************************%
%Done by Madeny Belkhiri from Duda Kvitsiani lab, Aarhus university,DANDRITE
%Version 3.0, date:18-03-2018
%Based on Kilosort spike sorting algorithm, GNU GENERAL PUBLIC LICENSE
%************************************************************************************%
%************************************************************************************%
clear all;
close all;
rng('default');
tic; % start timer
gpuDevice(1); % initialize GPU (will  erase any existing GPU arrays)
%run the config file filled by the user
Config_file;
% 
%If openephys format, need to convert into binary saved in example.dat
if strcmp(ops.datatype , 'openEphys')
    ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
    fprintf('Time %3.2f. minutes after converting from openephys to binary... \n', toc/60);
    %Get the conversion factor for bit in microvolt
    file =  strcat(ops.root,'\100_CH1.continuous');
    [~, ~, info] = load_open_ephys_data_faster(file);
    ops.bitmVolt = info.header.bitVolts;
end
%************************************************************************************%
%*******************************Preprocess data *************************************%
%***************************Filtering and whitening**********************************%
[rez, DATA] = preprocessData(ops);
fprintf('Time %3.2f minutes after preproccessing... \n', toc/60);
%************************************************************************************%
% fit templates and extract timestamps iteratively
%************************************************************************************%
k = 0;   %iterate for channel group
irun = 0;%iterate to accumulate detected spikes
kt = 0;  %iterate number of cluster
%initialize templates
%  dWU = zeros(ops.nt0,ops.Nchan,ops.Nfilt*ops.Nb_group);
 
%get spikes for each group

for ig = 1:8
     fprintf('Calculation for group: %d... \n',ig); %, toc
     channel = 1+k:1:ops.chan_per_group+k;
     DATAg = DATA(:,channel,:);
     %get initial templates
     [T] = Template_building3(DATAg,ops);
     
     if(isempty(nonzeros(T)))
        Ncl = 0;
        fprintf('No Templates found for group %d, you might want to lower ops.T_crit or ops.spkTh\n',ig); %, toc
     else
        %fit templates and extract timestamps, iteratively
        [st3]= NCC_overlap(DATAg, T,ops);  
        [~,~,Ncl] = size(T);
        dWU(:,channel,1+kt:Ncl+kt) = T;
        if(~isempty(st3))
            if(irun>0)
                incr = kt;
            else
                incr = 0;
            end
            st3(:,3) = st3(:,3) + incr;
            [L,~] = size(st3);
            rez.st3(irun+(1:L),:) = st3;
            irun = irun + L;            
        end
     end
     k = k + ops.chan_per_group;
     kt = kt + Ncl;
 end 
%************************************************************************************%
%Delete template that did not detect enough spikes
ops.Nfilt = kt;
%check number of templates actually used
list =[];
cluster_left = 1:ops.Nfilt;
for i=1:ops.Nfilt
    spike_id = find(rez.st3(:,3)==i);
    if(isempty(spike_id))
       list = [list;i];
    end    
end
%delete unused templates
dWU(:,:,list) =[];
cluster_left(list)= [];
ops.Nfilt = size(dWU,3);
for i=1:ops.Nfilt
    spike_id = find(rez.st3(:,3)==cluster_left(i));
    rez.st3(spike_id,3) = i;   
end
%************************************************************************************%
% Save template
[dWU,W, U, Weight] = SVD_template(dWU,ops.Nrank,ops.Chan_criteria );
rez.dWU = dWU;
rez.W = W;
rez.U = U;
rez.Weight = Weight;
rez.ops = ops;
%Vizualize all template used for detection dWU= W*d*U(SVD), in fact  W=W*d
figure;
subplot(1,2,1)
plot(W(:,:,1))
xlim([0 ops.nt0-1]);
title('W(:,:,1)')            
subplot(1,2,2)
imagesc(U(:,:,1))
title('U(:,:,1)')

%order spikes by timestamps
[~, isort]      = sort(rez.st3(:,1), 'ascend');
rez.st3         = rez.st3(isort,:);
%*****************************Post-process the data******************************%
%copy the detected spikes information in different matrix in order to
%post-process them and conserve the initial  values.
rez.st = rez.st3;
%************************************************************************************%
%Get the best channels for each clusters
[rez] = Chan_cluster(rez,ops.Chan_criteria);
% %************************************************************************************%
% %Merging Templates
% %************************************************************************************%
[rez] = merging(rez,ops.Threshold);
fprintf('Time %3.2f minutes after merging... \n', toc/60);
%************************************************************************************%
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3')
%************************************************************************************%
% %************************************************************************************%
[rez] = Waveform(rez,DATA);
fprintf('Time %3.2f minutes after mean waveform calculation... \n', toc/60);
%************************************************************************************%
%Calculate projection of each spikes onto its template
%************************************************************************************%
%Calculate waveforms
% [~,~,U,~] = SVD_template(rez.M_template,rez.ops.Nrank,rez.ops.Chan_criteria );
% rez.U = U;
% [rez] = Chan_cluster_M(rez,rez.ops.Chan_criteria);
% M_clust = length(unique(rez.st(:,end)));
% for j=1:M_clust
%        i_chan = rez.Chan{j};
%        rez.Merge_cluster{j,5} = i_chan;   
% end
[rez] = projection(rez,DATA);
%************************************************************************************%
%************************************************************************************%
%Calculate auto-correlogram 
%************************************************************************************%
[rez] = autocorrelogram(rez);
%************************************************************************************%
%Save the data
%************************************************************************************%
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3')
fprintf('Time %3.2f minutes for Full calculation... \n', toc/60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THE END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
