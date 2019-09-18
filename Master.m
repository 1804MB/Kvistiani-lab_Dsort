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
session = dir;
session = session(3:end);

for i = 1:length(session) % [8 11 18 21 24 26 27]*
sessionpath = [cd,'\',session(i).name];
ops = Config_file(sessionpath);
% 

%If openephys format, need to convert into binary saved in example.dat
if strcmp(ops.datatype , 'opeEphys') % dat opeEphys
    ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
    fprintf('Time %3.2f. minutes after converting from openephys to binary... \n', toc/60);
    %Get the conversion factor for bit in microvolt
    file =  strcat(ops.root,'\100_CH17.continuous');
%     [~, ~, info] = load_open_ephys_data_faster(file);
%     ops.bitmVolt = info.header.bitVolts;
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


for ig = 1:ops.Nb_group 
     fprintf('Calculation for group: %d... \n',ig); %, toc
      channel = 1+k:1:ops.chan_per_group+k;
%      channel = 13:16;
     DATAg = DATA(:,channel,:);
     %get initial templates
     [T] = Template_building(DATAg,ops);
%      if(ops.snr == 'O')
%         [T] = estimate_snr(rez,T);
%      end
     if(isempty(nonzeros(T)))
        Ncl = 0;
        fprintf('No Templates found for group %d, you might want to lower ops.T_crit or ops.spkTh\n',ig); %, toc
     else
        %fit templates and extract timestamps, iteratively
        [st,T] = NCC_overlap(DATAg, T,ops);  
%         [st,To] = check_template(DATAg,st,T,ops); 
        [~,~,Ncl] = size(T);
        dWU(:,channel,1+kt:Ncl+kt) = T;
        if(~isempty(st))
            if(irun>0)
                incr = kt;
            else
                incr = 0;
            end
            st(:,3) = st(:,3) + incr;
            [L,~] = size(st);
            rez.st(irun+(1:L),:) = st;
            irun = irun + L;            
        end
     end
     k = k + ops.chan_per_group;
     kt = kt + Ncl;
 end 
%************************************************************************************%
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
[~, isort]      = sort(rez.st(:,1), 'ascend');
rez.st         = rez.st(isort,:);
%*****************************Post-process the data******************************%
%************************************************************************************%
%Get the best channels for each clusters
[rez] = Chan_cluster(rez,ops.Chan_criteria);
% %************************************************************************************%
% %Merging Templates
% %************************************************************************************%
[rez] = merging(rez,ops.Threshold);
fprintf('Time %3.2f minutes after merging... \n', toc/60);
%************************************************************************************%
save(fullfile(ops.root,  'rezI.mat'), 'rez', '-v7.3')
%************************************************************************************%
% %************************************************************************************%
[rez] = Waveform(rez,DATA);
fprintf('Time %3.2f minutes after mean waveform calculation... \n', toc/60);
%************************************************************************************%
%Calculate projection of each spikes onto its template
%************************************************************************************%
[rez] = projection(rez,DATA);
%************************************************************************************%
%************************************************************************************%
%Calculate auto-correlogram 
%************************************************************************************%
[rez] = autocorrelogram(rez);
%************************************************************************************%
[rez]=probability(rez);
%Save the data
%************************************************************************************%
save(fullfile(ops.root,  'rezI.mat'), 'rez', '-v7.3')
fprintf('Time %3.2f minutes after projection_autoccorelation calculation... \n', toc/60);
% Spliting, cleaning and final merging clusters
%************************************************************************************%
[rez] = cleaning_clusters(rez,DATA);
save(fullfile(ops.root,  'rezF.mat'), 'rez', '-v7.3')
fprintf('Time %3.2f minutes for Full calculation... \n', toc/60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THE END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
end
