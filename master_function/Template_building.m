function [dWU] = Template_building(DATA,ops)
[NT, Nchan ,Nbatch] = size(DATA);
Nfilt = ops.Nfilt;
% Deleted the option to choose from data or not, now only fromdata.
% thus, the initiale template will always be calculated from the data.
if(ops.template =='Y')
    %     dWU0   = load(ops.template_file);
    load(ops.template_file);
%   dWU    = S_template;
    dWU    = rez.M_template;
elseif(ops.template=='G') % or it generates templates based on config file
    [WUinit] = init_template(ops,Nchan);
    dWU    = WUinit(:,:,1:Nfilt);
elseif(ops.template=='P') % %initilization using peak to peak
    dWU = [];
    X  = [];
    WV = [];
    irun = 0;
    P = randperm(Nbatch);
    Nbatch_use = ceil(ops.Nbatch_ratio*Nbatch);
    P = P(1:Nbatch_use);
    for ibatch = 1:Nbatch_use
        data   = single(DATA(:,:,P(ibatch)))/ops.scaleproc;
        dat    = -single(DATA(:,:,P(ibatch)))/ops.scaleproc;
        for k=1:Nchan
            [~, inx] = findpeaks(dat(ops.nt0:end-ops.nt0,k), 'MinPeakDistance',ops.nt0);
            inx = inx+ops.nt0;
            stdd = mad(dat(:,k));
            time = inx(dat(inx,k)>ops.T_crit*stdd);
            el = 1:NT;
            el(time) =  [];
            dat(el,k) = 0;
        end
        [tr,~] = max(dat,[],2);
        [~, inx] = findpeaks(tr, 'MinPeakDistance',ops.nt0);
        if(~isempty(inx))
            ind = int32(inx);
            ar =  data(ind,:);
            ind = int32(inx-ops.nt0min);
            inds = repmat(ind', ops.nt0, 1) + repmat(int32(1:ops.nt0)', 1, numel(ind));
            ard =  data(inds,:);
            ard = reshape(ard,[size(inds) Nchan]);
            ard = permute(ard,[2,1,3]);
            ars = (reshape(ard,length(inx),Nchan*ops.nt0));
            X(irun + (1:numel(inx)),:) = ar;
            WV(irun + (1:numel(inx)),:) = ars;
            irun  = irun +numel(inx);
        end
    end
    if(length(X)>ops.Ns_min)
        NK  = ops.Nfilt;
        idx = kmeans(X,NK,'Replicates',20,'MaxIter',250);
        %calculate the mean of each template
        cvg = 0;
        i   = 1;
        while(cvg==0)
            id = find(idx==i);
            if(length(id)>ops.Ns_min)
                wave = WV(idx==i,:);
                Mean_wave = mean(wave);
                PCP = wave*mean(wave)';
                if(PCP<4000)
                    PCP = interp1(1:length(PCP),PCP,1:ceil(length(PCP)/4000):length(PCP));
                end
                [~, pv] = HartigansDipSignifTest(sort(PCP), length(PCP));
                if(pv<ops.p_val_dip_test )
                    dWU = reshape(Mean_wave,[ops.nt0 Nchan]);
                    Nrank = 3;
                    [~,W, ~, ~,~] = SVD_template(dWU, Nrank,ops.Chan_criteria) ;
                    WVP = reshape(WV(id,:), [length(id) ops.nt0 Nchan]);
                    WVP = permute(WVP, [2 1 3]);
                    coefs = reshape(squeeze(W(:,1,1:Nrank))' *reshape(WVP, ops.nt0, []) , Nrank, numel(id), Nchan);
                    PC = permute(coefs, [3 1 2]);
                    %use the kurtosis to guess on which channel to make the
                    %cut
                    sig_for = zeros(Nchan,Nrank);
                    for k=1:Nchan
                        for ij = 1:Nrank
                            sig_for(k,ij) =  kurtosis(PC(k,ij,:));
                        end
                    end
                    [~,ichan]=min(min(sig_for,[],2));
                    options = statset('MaxIter',500);
                    warning('off','stats:gmdistribution:FailedToConvergeReps');
                    try   % try catch loop to avoid problems with GMM.  line 95 ends
                        GMModel = fitgmdist(squeeze(PC(ichan,:,:))',2,'Options',options,'Replicates',100,'RegularizationValue',0.01);
                        clusterX = cluster(GMModel,squeeze(PC(ichan,:,:))');
                        id2 = find(clusterX==2);
                        idx(id(id2)) = length(unique(idx))+1;
                    catch         
                        i = i +1;
                        NK = length(unique(idx));
                        if(i>NK)
                            cvg = 1;
                        end
                    end
                    
                else
                    i = i +1;
                    NK = length(unique(idx));
                    if(i>NK)
                        cvg = 1;
                    end
                end
            else
                if(isempty(id))
                    cvg =1;
                end
                WV((id),:) = [];
                idx(id) = [];
                list_left = unique(idx);
                for jk = 1:length(list_left)
                    idL = find(idx==list_left(jk));
                    idx(idL) = jk;
                end
            end
        end
        NK = length(unique(idx));
        dWU =[];
        for i=1:NK
            dWU(i,:) = mean(WV(idx==i,:));
        end
        dWU = reshape(dWU,[NK ops.nt0 Nchan]);
        dWU = permute(dWU,[2 3 1]);
        
    end
else%Get template using Kilosort initialization
    i0  = 0;
    wPCA = ops.wPCA;
    uproj = zeros(1e6,  size(wPCA,2) * Nchan, 'single');
    %Apply the whitening matrix to the filtered data and get the PC projection of the spikes
    %used to initialize the template used in the next function FitTemplate.m
    for ibatch = 1:Nbatch
        dataRAW = single(DATA(:,:,ibatch))/ops.scaleproc;
        %find isolated spikes
        [row, col, ~] = isolated_peaks(dataRAW, ops.loc_range, ops.long_range, ops.spkTh,ops);
        % find their PC projections
        uS = get_PCproj(dataRAW, row, col, wPCA, ops.maskMaxChannels,ops.nt0min);
        uS = permute(uS, [2 1 3]);
        uS = reshape(uS,numel(row), Nchan * size(wPCA,2));
        
        if i0+numel(row)>size(uproj,1)
            uproj(1e6 + size(uproj,1), 1) = 0;
        end
        %Save the PC projection on the CPU
        uproj(i0 + (1:numel(row)), :) = gather_try(uS);
        i0 = i0 + numel(row);
    end
    uproj(i0+1:end, :) = [];
    %does a scaled kmeans to determine the  Nfilt clusters template
    WUinit = optimizePeaks(ops,uproj,Nchan);
    dWU    = WUinit(:,:,1:Nfilt);
end
if(~isempty(dWU))
    [dWU,top] = SVD_topchan(dWU,ops.Chan_criteria);
    dWU = alignWU(dWU,top, ops.nt0min); %make sure that the deep is set at nt0min vital for substracting the template from signal
end