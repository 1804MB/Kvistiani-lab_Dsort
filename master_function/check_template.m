function  [st,T]=check_template(DATA,st,T,ops)
M_clust = length(unique(st(:,3)));
nt0  = ops.nt0;
[NT,Nchan,~] = size(DATA);

Nbatch = length(unique(st(:,end)));
batch = unique(st(:,end));
Nrank = 3;%default value for projection hard coded
%calculate the mean of each template
irun = 0;
if(length(st(:,1))>ops.Ns_min)
    for j=1:Nbatch
        id = find(st(:,end)==batch(j));
        time = double(st(id,1));
        if(batch(j)==1)
            offset = 0;
        else
            offset = ops.ntbuff;
        end
        ind = time - ops.nt0min + offset - (NT - ops.ntbuff)* (double(batch(j))-1);   
        ind = int32(ind);
        inds = repmat(ind', nt0, 1) + repmat(int32(1:nt0)', 1, numel(ind));
        dat = single(DATA(:,:,batch(j)))/ops.scaleproc;
        ar  = dat(inds(:),:);        
        ar  = reshape(ar,[size(inds) Nchan]);
        ar = permute(ar,[2,1,3]);
        ar = (reshape(ar,length(time),Nchan*ops.nt0));
        WV(irun + (1:numel(time)),:) = ar;
        irun  = irun +numel(time);     
    end

    idx = int32(st(:,3));
        %calculate the mean of each template
        for i=1:M_clust
            id = find(idx==i);
            if(length(id)>ops.Ns_min)
                wave = WV(idx==i,:);
                Mean_wave = mean(wave);
                PCP = wave*mean(wave)';
                [~, pv] = HartigansDipSignifTest(sort(PCP), length(PCP));
                [N,edge]=histcounts(PCP);
                figure;
                bar(edge(1:end-1),N);
                title(['Final distribution, cluster=',num2str(i), 'pv:', num2str(pv)]) 
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
                    GMModel = fitgmdist(squeeze(PC(ichan,:,:))',2,'Options',options,'Replicates',100);
                    clusterX = cluster(GMModel,squeeze(PC(ichan,:,:))');
                    id2 = find(clusterX==2);
                    idx(id(id2)) = length(unique(idx))+1;
                end
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
            T  = dWU;
            st(:,3) = idx;
end
 