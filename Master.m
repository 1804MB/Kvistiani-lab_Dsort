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
rootdir = cd;
session = dir;
session = session(3:end);
sessions2cluster =  [1];
for i =  sessions2cluster % [8 11 18 21 24 26 27]*
    sessionpath = [cd,'\',session(i).name];
    ops = Config_file(sessionpath);
    
    
%     If openephys format, need to convert into binary saved in example.dat
    if strcmp(ops.datatype , 'openEphys') % dat opeEphys
        ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
        fprintf('Time %3.2f. minutes after converting from openephys to binary... \n', toc/60);
        %Get the conversion factor for bit in microvolt
        file =  strcat(ops.root,'\100_CH17.continuous');
        %     [~, ~, info] = load_open_ephys_data_faster(file);
        %     ops.bitmVolt = info.header.bitVolts;
    end
    
    %************************************************************************************%
    %********************if template is precomputed from another session*****************%
    %************************************************************************************%
    
    if ops.template == 'Y'
        ops.Nb_group = 1;
        ops.chan_per_group = 32;
        ops.template_file = 'O:\ST_Duda\Maria\MM003\2019-02-21_16-03-54\rezFMerged.mat';
        ops.criteria_NCC = 0.85;
        ops.Ncycle = 4;
    end
    
    %*************************************************************************************%
    
    
    
    %************************************************************************************%
    %*******************************Preprocess data *************************************%
    %***************************Filtering and whitening**********************************%
    [rez, DATA] = preprocessData(ops);
    fprintf('Time %3.2f minutes after preproccessing... \n', toc/60);
    %************************************************************************************%
    % load presorted templates if they exist
    %************************************************************************************%
    if exist([ops.root,'\rezJunk.mat'])
         
        
    
%     if exist([ops.root,'\rezFinal.mat'])
%         load ([ops.root,'\rezFinal.mat'])
%     elseif exist([ops.root,'\rezFMerged.mat'])
%         load ([ops.root,'\rezFMerged.mat'])
%         ops.Ncycle = 0;
%     elseif exist([ops.root,'\rezF2.mat'])
%         load ([ops.root,'\rezF2.mat'])
%         ops.Ncycle = 5;   
%     elseif exist([ops.root,'\rezF.mat'])
%         load ([ops.root,'\rezF.mat'])
%         ops.Ncycle = 7;
%     elseif exist([ops.root,'\rezI.mat'])
%         load ([ops.root,'\rezI.mat'])
%         ops.Ncycle = 10;
    else
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
        rez.st         =  rez.st(isort,:);
        %*****************************Post-process the data******************************%
        %************************************************************************************%
        %Get the best channels for each clusters
        [rez] = Chan_cluster(rez,ops.Chan_criteria);
        % %************************************************************************************%
        % %Merging Templates
        % %************************************************************************************%
        [rez] = merging(rez,0.99);
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
    end
    
    %************************************************************************************%
    %************************* clean, merge and delete clusters**************************%
    rez.ops = ops;
    rez = CleaningMergingDeletingClusters(rez,DATA);
    %************************************************************************************%
    
    %If rez file already exists continue only cleaning and merging steps
    %************************************************************************************%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THE END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Split, merge and delete clusters%%%%%%%%%%%%%%%%%%%%%
function rez = CleaningMergingDeletingClusters(rez,DATA)
% Spliting and deleting clusters that have low amplitude
%************************************************************************************%
[rez] = cleaning_clusters(rez,DATA);
save(fullfile(rez.ops.root,  'rezF.mat'), 'rez', '-v7.3')
fprintf('Time %3.2f minutes after cleaning... \n', toc/60);

%   merging based on WV correlation. repeat untill there is nothing to merge
%************************************************************************************%
NClust = length(unique(rez.st(:,end)));
i = 1;
k = 1;
Threshold = 0.9 ; %rez.ops.Threshold;
while i > 0
    [rez] = merging_T(rez,Threshold,k);
    [rez] = recompute(rez,DATA);
    NClust = [NClust length(unique(rez.st(:,end)))];
    save(fullfile(rez.ops.root,  'rezFMerged.mat'), 'rez', '-v7.3')
    i = NClust(end-1) - NClust(end);
    k = k + 0;
end
Threshold = 0.9 ; %rez.ops.Threshold;
[rez] = merging_S(rez,Threshold);
[rez] = recompute(rez,DATA);
% fprintf('Time %3.2f minutes after merging... \n', toc/60);

%   Delete ovewrlapping clusters
%************************************************************************************%

% [rez]=Remove_T(rez,Threshold);
[rez]=Remove_S(rez,Threshold);

[rez] = recompute(rez,DATA);
save(fullfile(rez.ops.root,  'rezFinalS.mat'), 'rez', '-v7.3')
fprintf('Time %3.2f minutes for Full calculation... \n', toc/60);
close all;
end
