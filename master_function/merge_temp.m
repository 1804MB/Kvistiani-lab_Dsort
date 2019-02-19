function [rez]=merge_temp(rez,Threshold)

%Variable initialization
SF = rez.ops.fs;                        %sampling frequency to convert from sampling rate to second
Nclusters = length(unique(rez.st(:,3)));              %Number of initial clusters from kilosort
NG_clus = Nclusters;             %Number of clusters
nt0 = rez.ops.nt0;
Chan = rez.Chan;
%extract the 4th largest channels for each template

% MERGING CLUSTERS (THIRD PART): SUPPORT VECTOR MACHINES
% Compute the number of spikes needed to support the vector that can maximize 
%the margin between two clusters that have the same channels for largest projection
rez.score = zeros(NG_clus,NG_clus);
h = waitbar(0, 'NCC for merging...');
for i=1:NG_clus
    id =find(rez.st(:,end)==i);
     Chan1 = Chan{i};
     if(length(Chan1)==length(Chan2) )     
                if(Chan1(1)==Chan2(1))
       
        Progress = (i/NG_clus);
        waitbar(Progress)
        as = rez.dWU(:,Chan1,i);
        as = reshape(as,[],length(Chan1)*nt0);
        for j=i:NG_clus
            id2 =find(rez.st(:,end)==j);
 
                if(isempty(id2)||overlap(i,j)==0)
                else
                    bs = rez.dWU(:,Chan1,j);
                    bs = reshape(bs,[],length(Chan1)*nt0);
                    score = xcorr(as,bs,0,'coeff');
                    rez.score(i,j) = score;
                    if(score>Threshold)
                        rez.st(id2,end) = i;
                        
                    end
                end
        end
    end
end
close(h)

M_clust =length(unique(rez.st(:,end))) %number of merged clusters
old_Clusters = unique(rez.st(:,end));
for i=1:M_clust
    spike_id = find(rez.st(:,end) ==old_Clusters(i));
    rez.st(spike_id,end) = i;
end
%Assign a cluster number to the merged clusters
Merge_cluster = cell(M_clust,5);
Nbatch = length(unique(rez.st(:,end-1)));


for i=1:M_clust
   spike_id = find(rez.st(:,end) ==i);
   Merge_cluster{i,1} = unique(rez.st(spike_id,3));
   Merge_cluster{i,2} = length(spike_id);
   Merge_cluster{i,3} = i;
   Merge_cluster{i,4} = 1/(length(spike_id)/(rez.ops.NT*rez.temp.Nbatch/(SF))); 
   el = Merge_cluster{i,1};
   i_chan = Chan{el(1)};
   Merge_cluster{i,5} = i_chan;   
end

rez.Merge_cluster = Merge_cluster;



