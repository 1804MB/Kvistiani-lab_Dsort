function [rez] = merge_cluster(rez,cluster1,cluster2)

id = find(rez.st(:,end)==cluster2);
rez.st(id,end) = cluster1;   
rez.Merge_cluster(cluster2,:) = [];

M_clust = length(unique(rez.st(:,end)));
for i=1:M_clust
   spike_id = find(rez.st(:,end)==i);
   rez.Merge_cluster{i,1} = unique(rez.st(spike_id,3));
   rez.Merge_cluster{i,2} = length(spike_id);
   rez.Merge_cluster{i,3} = i;
   rez.Merge_cluster{i,4} = 1/(length(spike_id)/(rez.ops.NT*rez.temp.Nbatch/(rez.ops.fs))); 
   el = rez.Merge_cluster{i,1};
   i_chan = rez.Chan{el(1)};
   rez.Merge_cluster{i,5} = i_chan;   
end
   
