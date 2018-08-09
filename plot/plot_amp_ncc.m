function plot_amp_ncc(rez)
M_clust = length(unique(rez.st(:,end)));
for i=1:M_clust
     id = find(rez.st(:,end)==i);
     amp = rez.st(id,5);
     ncc_score =rez.st(id,4);
     [Ncount,edgesS] = histcounts(amp,min(amp):10:max(amp));
    [Ncc,edgesSc] = histcounts(ncc_score,min(ncc_score):0.005:max(ncc_score));
     figure;
     subplot(1,2,1)
      bar(edgesS(1:end-1),Ncount(1:end))
     lim = edgesS(end-1);
     ylabel('Spike count')
     xlabel('amplitude in mV')
   
     subplot(1,2,2)
     bar(edgesSc(1:end-1),Ncc(1:end))
     lim = edgesSc(end-1);
     ylabel('Spike count')
     xlabel('score')
     xlim([0.8 lim])      
     title(['Cluster:', num2str(i)]);
 end