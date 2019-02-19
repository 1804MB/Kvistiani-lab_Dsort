function run_peth(cluster1,cluster2,rez)

edges = -50:1:50;
    Spike_id = find(rez.st(:,3)==cluster1);
    M1 = rez.st(Spike_id,2);
    Spike_id2 = find(rez.st(:,3)==cluster2);
    M2 = rez.st(Spike_id2,2); 
    TrialN = length(M2);    
    C = zeros(TrialN,length(edges)-1);
    SpikePETH_M = zeros(length(edges)-1,1);
    %calculate  PETH
    for iTrial = 1:TrialN
         allspikes = M1 - M2(iTrial);
         C(iTrial,:) = histcounts(allspikes,edges);
         SpikePETH_M(:,1) = SpikePETH_M(:,1)  + C(iTrial,:)';
    end
figure;
     bar(-50:1:48,SpikePETH_M(1:99,1))
     title(['PETH of cluster:',num2str(cluster1),' and:', num2str(cluster2)])
     ylabel('Spike count')
     xlabel('Time (ms)')
     xlim([-50 50])
     hold on;