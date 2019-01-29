function plot_ncc(rez,i)
     Spike_id = find(rez.st(:,end)==i);
     ncc_score =rez.st(Spike_id,4);
     [Ncc,edgesSc] = histcounts(ncc_score,min(ncc_score):0.001:max(ncc_score));
     figure
     bar(edgesSc(1:end-1),Ncc(1:end))
     ylabel('Spike count')
     xlabel('NCC score')
     xlim([min(ncc_score) 1])      
     title(['Cluster:', num2str(i)]);
