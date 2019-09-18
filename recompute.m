function [rez] = recompute(rez,DATA)
%recalculate waveforms, projections and autocorrelogram
%Calculate waveforms
cluster = unique(rez.st(:,end));
M_clust = length(unique(rez.st(:,end)));
for j=1:M_clust
       spike_id = find(rez.st(:,end)==cluster(j));  
       rez.st(spike_id,end) =  j;
end
%************************************************************************************%
fprintf('Recalculate waveforms\n');
[rez] = Waveform_M(rez,DATA);
[~,~,U,~] = SVD_template(rez.M_template,rez.ops.Nrank,rez.ops.Chan_criteria );
rez.U = U;
[rez] = Chan_cluster_M(rez,rez.ops.Chan_criteria);

cluster = unique(rez.st(:,end));
M_clust = length(unique(rez.st(:,end)));
fprintf('Number of clusters %d\n',M_clust); 
rez.Merge_cluster = cell(M_clust,5);
    for j=1:M_clust
       spike_id = find(rez.st(:,end)==cluster(j));
       rez.Merge_cluster{j,1} = unique(rez.st(spike_id,3));
       rez.Merge_cluster{j,2} = length(spike_id);
       rez.Merge_cluster{j,3} = j;
       rez.Merge_cluster{j,4} = 1/(length(spike_id)/(rez.ops.NT*rez.temp.Nbatch/(rez.ops.fs))); 
       el = rez.Merge_cluster{j,3};
       i_chan = rez.Chan{el(1)};
       rez.Merge_cluster{j,5} = i_chan;   
       rez.st(spike_id,end) =  j;
    end
fprintf('Recalculate projections\n');
%Calculate probability of each spikes to belong to its cluster
[rez] = projection(rez,DATA);
%************************************************************************************%
%Vizualize waveforms, auto-correlogram and amplitude histogram
%************************************************************************************%
[rez]=autocorrelogram(rez);
[rez]=probability(rez);
end