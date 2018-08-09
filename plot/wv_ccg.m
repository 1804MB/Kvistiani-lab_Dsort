function wv_ccg(rez,i)
     Spike_id = find(rez.st(:,end)==i);
     SpikePETH_M =  rez.PETH(:,i);
     Nchan = length(cell2mat(rez.Merge_cluster(i,5)));
    %Get histogram of spikes amplitudes
     amplitudes = (rez.st(Spike_id,5));
     [Ncount,edgesS] = histcounts(amplitudes,min(amplitudes):10:max(amplitudes));
     h=figure;
      movegui(h,'south');
     subplot(1+ceil(Nchan/4),4,[1,2])
     bar(-50:1:48,SpikePETH_M(1:99,1))
     title(['PETH of cluster:',num2str(i),' with:', num2str(length(Spike_id)), ' spikes and refrac:', num2str(SpikePETH_M(51,1))])
     ylabel('Spike count')
     xlabel('Time (ms)')
     xlim([-50 50])
     hold on;
     subplot(1+ceil(Nchan/4),4,[3,4])
     bar(edgesS(1:end-1),Ncount(1:end))
     lim = edgesS(end-1);
     ylabel('Spike count')
     xlabel('Amplitude in mV')     
     Channel = rez.Merge_cluster{i,5};
     for k = 1:Nchan
        amean = [squeeze(rez.M_template(:,Channel(k), i))'; squeeze(rez.M_template(:,Channel(k), i))'];
        astd = [squeeze(rez.M_std_template(:,Channel(k), i))'; squeeze(rez.M_std_template(:,Channel(k), i))'];       
       subplot(1+ceil(Nchan/4),4,k+4)
       stdshade_sorting(amean, astd, 0.3, 'b')

       title([ 'Channel: ', num2str(Channel(k))])
       xlim([0 rez.ops.nt0]);
     end 

