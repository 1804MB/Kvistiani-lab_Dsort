function [rez]=Remove_T(rez,Threshold)

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
for i= 1:NG_clus
    Progress = (i/NG_clus);
    waitbar(Progress)
    id =find(rez.st(:,end)==i);
    %   Chan1 = Chan{i};
    Chan1 = rez.Merge_cluster{i,end};
    if(~isempty(id))
        cl1 = rez.st(id,2);
        as = rez.M_template(:,:,i); % use all channels
        as = reshape(as,[],rez.ops.NchanTOT*nt0);
        %         as = rez.M_template(:,Chan1,i); % use only top channels
        %         as = reshape(as,[],length(Chan1)*nt0);
        for j=1:NG_clus
            id2 =find(rez.st(:,end)==j);
            if(~isempty(id2))
                if(i ~= j)
                    %                   Chan2 = Chan{j};
                    Chan2 = rez.Merge_cluster{j,end};
                    cl2 = rez.st(id2,2);
                    Raster = zeros(length(cl2),100); % initialize CCG
                    bs = rez.M_template(:,:,j); % use all channels for merging
                    bs = reshape(bs,[],rez.ops.NchanTOT*nt0);
                    %                       bs = rez.M_template(:,Chan2,j); % use only top channels for merging
                    %                       bs = reshape(bs,[],length(Chan2)*nt0);
                    score = xcorr(as,bs,0,'coeff');
                    rez.score(i,j) = score;
                    for icell = 1:length(cl2)
                        allspikes = cl1 - cl2(icell);
                        spikes = ceil(allspikes(allspikes > -50 & allspikes < 50) + 50);
                        Raster(icell,spikes) = 1;
                    end
                    MinSpikeCount = min(length(find(rez.st(:,end)==i)),length(find(rez.st(:,end)==j)));
                    if score > Threshold && length(find(rez.st(:,end)== i)) >= length(find(rez.st(:,end)== j))
                        rez.st(id2,:) = [];
                    elseif sum(Raster(:,50)) > MinSpikeCount/5 && length(find(rez.st(:,end)== i)) >= length(find(rez.st(:,end)== j)) || length(find(rez.st(:,end)== j)) < SpikeCountThreshold
                        rez.st(id2,:) = [];
                    end
                    %                     if length(Chan1) == length(Chan2)
                    %                         if Chan1 == Chan2 & sum(Raster(:,50)) > mean(sum(Raster)) & length(find(rez.st(:,end)== i)) >= length(find(rez.st(:,end)== j))
                    %                             rez.st(id2,:) = [];
                    %                         end
                    %                     end
                end
            end
        end
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