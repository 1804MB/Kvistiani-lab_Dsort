function [rez]=merging(rez,Threshold)

%Variable initialization
SF = rez.ops.fs;                        %sampling frequency to convert from sampling rate to second
Nclusters = length(unique(rez.st(:,3)));              %Number of initial clusters from kilosort
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

h = waitbar(0, 'Shuffled data for refractory period');
for i=1:NG_clus
    Progress = (i/NG_clus);
    waitbar(Progress)
    id =find(rez.st(:,3)==i);
    if(~isempty(id))
    cl1 = rez.st(id,2);  
    Chan1=Chan{i};
    for j=i:NG_clus
         Chan2=Chan{j};
        id2 =find(rez.st(:,3)==j);
        if(~isempty(id2))
        if(i==j)
            overlap(i,j) = 1;
            overlap(j,i) = 1;
        else
            if(length(Chan1)==length(Chan2) )     
                if(Chan1(1)==Chan2(1))
                    cl2 = rez.st(id2,2);  
                    t0 = [cl1;cl2];
                    count_s = length(t0); % Total number of spikes           
                    for k = 1:N
                        shuffle = t0 + randi([-s_int s_int],length(t0),1);
                        tcor(k) =length(unique(shuffle));
                    end
                    Init_value = count_s - length(unique(t0)); %Number of spikes violating the refractory period from the real data
                    CI_low = count_s - prctile(tcor,pour) ;   %Number of spikes violating the refractory period using the the confidence interval  
                    if(Init_value<CI_low)
                        overlap(i,j) = 1;
                        overlap(j,i) = 1;
                    else
                        overlap(i,j) = 0;
                        overlap(j,i) = 0;
                    end 
                end
            end
        end
        end
    end
    end
end
close(h)
% MERGING CLUSTERS (THIRD PART): SUPPORT VECTOR MACHINES
% Compute the number of spikes needed to support the vector that can maximize the margin between two clusters that have the same channels for largest projection
rez.st(:,end+1) = rez.st(:,3);
rez.score = zeros(NG_clus,NG_clus);
h = waitbar(0, 'NCC for merging...');
for i=1:NG_clus
    id =find(rez.st(:,end)==i);
     Chan1 = Chan{i};
    if(isempty(id))
    else
       
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



