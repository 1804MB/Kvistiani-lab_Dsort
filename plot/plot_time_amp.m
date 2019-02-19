function plot_time_amp(rez,i)
     Spike_id = find(rez.st(:,end)==i);     
     f=figure;
     channel  = rez.Merge_cluster{i,5};
     Nchan = length(channel);
     for k=1:Nchan
       subplot(Nchan,1,k)   
       plot(rez.st(Spike_id,2),rez.best_channel_M(Spike_id,channel(k)),'.')
       xlabel('Time in ms')
       ylabel('Amplitude (\muV) ')   
       title(['Channel:',num2str(channel(k))])
     end
     movegui(f,'southeast');
     hold off
