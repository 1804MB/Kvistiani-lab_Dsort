function  [rez,chan] = clean_cluster(rez)
M_clust = length(unique(rez.st(:,end)));
chan = NaN(1,2*M_clust);
Nrank = 3;
del = [];
for i =1:M_clust
    Spike_id = find(rez.st(:,end)==i);
    channel  = rez.Merge_cluster{i,5};
    Nchan = length(channel);
    PC = rez.PC{i};
    pv = zeros(Nchan,Nrank);
    if(length(Spike_id)>2000)
        for k=1:Nchan
            for j=1:Nrank
                amplitude = squeeze(PC(k,j,:));
                amplitude = interp1(1:length(amplitude),amplitude,1:ceil(length(amplitude)/2000):length(amplitude));
                [~, pv(k,j)] = HartigansDipSignifTest(sort(amplitude), length(amplitude));
            end
        end
        [~,ichan]=min(min(pv,[],2));
        if(pv(ichan,1)<rez.ops.p_val_dip_test ||pv(ichan,2)<rez.ops.p_val_dip_test ||pv(ichan,3)<rez.ops.p_val_dip_test )
            options = statset('MaxIter',500);
            warning('off','stats:gmdistribution:FailedToConvergeReps');
            GMModel = fitgmdist(squeeze(PC(ichan,:,:))',2,'Options',options,'Replicates',100,'RegularizationValue',0.01);
            clusterX = cluster(GMModel,squeeze(PC(ichan,:,:))');
            id1 = find(clusterX==1);
            id2 = find(clusterX==2);
            plot_pc_init_two(PC,Nchan,i,id1,id2)
            rez.st(Spike_id(id2),end) = max(unique(rez.st(:,end)))+1; % length was instead of max
            chan(i) = channel(ichan); chan(max(unique(rez.st(:,end)))+1) = channel(ichan);
        else
            chan(i)  = channel(ichan); %rez.Merge_cluster{i,5}(1);
        end
    else
        chan(i)  = rez.Merge_cluster{i,5}(1);
    end
end
chan(find(isnan(chan))) = [];


