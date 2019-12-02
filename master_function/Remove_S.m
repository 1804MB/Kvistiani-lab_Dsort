function [rez]=Remove_S(rez,Threshold)

%Variable initialization
SF = rez.ops.fs;                        %sampling frequency to convert from sampling rate to second
Nclusters = length(unique(rez.st(:,end)));              %Number of initial clusters from kilosort

NG_clus = Nclusters;             %Number of clusters
nt0 = rez.ops.nt0;
Chan = rez.Chan;
SpikeCountThreshold = 1000;
%extract the 4th largest channels for each template


%********************************************************************%
%*******************Calculate cross-correlograms**********************%
%********************************************************************%
%choose overlapping waveforms and delete one that has less spikes

h = waitbar(0, 'Delete redundant templates');
del = [];
for i= 1:NG_clus
    Progress = (i/NG_clus);
    waitbar(Progress)
    if  length(find(rez.st(:,end)==i)) < 1000
        rm = find(rez.st(:,end)==i);
        rez.st(rm,:) = [];
        rez.best_channel_M(rm,:)=[];
        del = [del;rez.Merge_cluster{i,1}];
    end
end
    
   

close(h)
% MERGING CLUSTERS (THIRD PART): SUPPORT VECTOR MACHINES
% Compute the number of spikes needed to support the vector that can maximize
%the margin between two clusters that have the same channels for largest projection


M_clust =length(unique(rez.st(:,end))) %number of merged clusters
old_Clusters = unique(rez.st(:,end));
for i=1:M_clust
    spike_id = find(rez.st(:,end) ==old_Clusters(i));
    rez.st(spike_id,end) = i;
end
%Assign a cluster number to the merged clusters
Merge_cluster = cell(M_clust,5);

for i=1:M_clust
    spike_id = find(rez.st(:,end) ==i);
    Merge_cluster{i,1} = unique(rez.st(spike_id,end));
    Merge_cluster{i,2} = length(spike_id);
    Merge_cluster{i,3} = i;
    Merge_cluster{i,4} = 1/(length(spike_id)/(rez.ops.NT*rez.temp.Nbatch/(SF)));
    el = Merge_cluster{i,1};
    i_chan = Chan{el(1)};
    Merge_cluster{i,5} = i_chan;
end

rez.Merge_cluster = Merge_cluster;