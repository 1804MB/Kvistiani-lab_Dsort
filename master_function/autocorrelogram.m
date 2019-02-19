function [rez]=autocorrelogram(rez)
edges =-50:1:50;
M_clust = length(unique(rez.st(:,end)));

rez.PETH = zeros(100,M_clust);
    for i=1:M_clust
    Spike_id = find(rez.st(:,end)==i);
    M1 = rez.st(Spike_id,2);
    M2 = M1; 
    TrialN = length(M2);    
    C = zeros(TrialN,length(edges)-1);
    SpikePETH_M = zeros(length(edges)-1,1);
    %calculate  PETH
    for iTrial = 1:TrialN
         allspikes = M1 - M2(iTrial);
         C(iTrial,:) = histcounts(allspikes,edges);
         SpikePETH_M(:,1) = SpikePETH_M(:,1)  + C(iTrial,:)';
    end
     SpikePETH_M(51,1) =  SpikePETH_M(51,1) - length(M1);
     rez.PETH(:,i) = squeeze(SpikePETH_M(:,1)); 
    end
