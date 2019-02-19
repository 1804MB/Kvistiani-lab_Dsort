function plot_amp(rez,i)
     Spike_id = find(rez.st(:,end)==i);     
     f=figure;
     channel  = rez.Merge_cluster{i,5};
     Nchan = length(channel);
     for k=1:Nchan
       amplitudes = rez.best_channel_M(Spike_id,channel(k));
       [Ncount,edgesS] = histcounts(amplitudes);
       subplot(Nchan,1,k)   
       bar(edgesS(1:end-1),Ncount(1:end))
       ylabel('Spike count')
       xlabel(['Amplitude (\muV) of Channel:',num2str(channel(k))])   
     end
     movegui(f,'northeast');
     hold off
