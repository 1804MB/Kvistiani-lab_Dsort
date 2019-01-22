function plot_amp_ncc(rez,i)
     Spike_id = find(rez.st(:,end)==i);
     SpikePETH_M =  rez.PETH(:,i);
     Nchan = length(cell2mat(rez.Merge_cluster(i,5)));
    %Get histogram of spikes amplitudes
     amplitudes = rez.st(Spike_id,5);
     [Ncount,edgesS] = histcounts(amplitudes);figure;

     subplot(1,2,1)
     bar(edgesS(1:end-1),Ncount(1:end))
     ylabel('Spike count')
     xlabel('Amplitude in \muV')   
     ncc_score =rez.st(Spike_id,4);
     [Ncc,edgesSc] = histcounts(ncc_score,min(ncc_score):0.001:max(ncc_score));
     subplot(1,2,2)
     bar(edgesSc(1:end-1),Ncc(1:end))
     lim = edgesS(end-1);
     ylabel('Spike count')
     xlabel('NCC score')
     xlim([min(ncc_score) 1])      
     title(['Cluster:', num2str(i)]);
