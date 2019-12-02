function [rez]=merging_T(rez,Threshold,k)

%Variable initialization
SF = rez.ops.fs;                        %sampling frequency to convert from sampling rate to second
Nclusters = length(unique(rez.st(:,end)));              %Number of initial clusters from kilosort

NG_clus = Nclusters;             %Number of clusters
nt0 = rez.ops.nt0;
Chan = rez.Chan;
%extract the 4th largest channels for each template


%********************************************************************%
%*******************Calculate cross-correlograms**********************%
%********************************************************************%
N = 200;          %Number of iteration (number of shuffling)
pour = 1;           %Value for the confidence interval calculated via percentile function in %
tcor = zeros(N,1);  %Vector contianing the number of spikes violating the refractory period at each shuffled iteration
%Shuffle interval between [-s_int s_int] in ms
s_int = 50;
overlap =  zeros(NG_clus,NG_clus);

h = waitbar(0, 'Shuffled data for refractory period and NCC');
for i=1:NG_clus
    Progress = (i/NG_clus);
    waitbar(Progress)
    id =find(rez.st(:,end)==i);
    if(~isempty(id))
        cl1 = rez.st(id,2);
        as = rez.M_template(:,:,i);
        as = reshape(as,[],rez.ops.NchanTOT*nt0);
        counter = 0;
        for j=i:NG_clus
            id2 =find(rez.st(:,end)==j);
            if(~isempty(id2))
                if(i ~= j)
                    cl2 = rez.st(id2,2);
                    Raster = zeros(length(cl2),100);
                    bs = rez.M_template(:,:,j); % use all channels for merging
                    bs = reshape(bs,[],rez.ops.NchanTOT*nt0);
                    score = xcorr(as,bs,0,'coeff');
                    rez.score(i,j) = score;
                    if(score > Threshold)
                        for icell = 1:length(cl2)
                            allspikes = cl1 - cl2(icell);
                            spikes = ceil(allspikes(allspikes > -50 & allspikes < 50) + 50);
                            Raster(icell,spikes) = 1;
                        end
                        if  sum(Raster(:,50)) < k  && counter < 1 %mean(sum(Raster))/100
                            rez.st(id2,end) = i;
                            counter = counter + 1;
                        end
                    end
                end
            end
        end
    end
end
close(h)
% MERGING CLUSTERS (THIRD PART): SUPPORT VECTOR MACHINES
% Compute the number of spikes needed to support the vector that can maximize
%the margin between two clusters that have the same channels for largest projection


M_clust = length(unique(rez.st(:,end))) %number of merged clusters
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