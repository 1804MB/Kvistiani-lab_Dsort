function  [rez] = Waveform_M(rez,DATA)
M_clust = length(unique(rez.st(:,end)));
nt0 = rez.ops.nt0;
Nchan = rez.ops.Nchan;
%calculate the average template of the merged clusters
M_template = zeros(rez.ops.nt0,Nchan,M_clust);
M_std_template = zeros(rez.ops.nt0,Nchan,M_clust);
suma = zeros(nt0,Nchan,M_clust);    %To calculate the average template
std_sum = zeros(nt0,Nchan,M_clust);    %To calculate the average template
nb_spikes = zeros(M_clust,1);       %Number of spikes per cluster

batch = length(unique(rez.st(:,end-1)));
irun = 0;
h = waitbar(0, 'Computing mean waveforms...');
%calculate the mean of each template
    for j=1:batch
        Progress = (j/batch);
        waitbar(Progress)
        id = rez.st(:,end-1)==j;
        time = double(rez.st(id,1));
        id_M = int32(rez.st(id,end));
         id_T = int32(rez.st(id,3));
        if(j==1)
            offset = 0;
        else
            offset = rez.ops.ntbuff;
        end
        ind = time - rez.ops.nt0min + offset - (rez.ops.NT - rez.ops.ntbuff)* (double(j)-1);   
        ind = int32(ind);
        inds = repmat(ind', nt0, 1) + repmat(int32(1:nt0)', 1, numel(ind));
        dat = single(DATA(:,:,j))/rez.ops.scaleproc;
        ar = dat(inds(:),:);
        ar = reshape(ar,[size(inds) Nchan]);
        rez.best_channel_M(irun + (1:numel(time)), :) =  squeeze(ar(rez.ops.nt0min,:,:)-ar(1,:,:))*rez.ops.bitmVolt*rez.ops.scaleproc;
        for il=1:M_clust
            idd = find(id_M ==il);                     
            suma(:,:,il) =  suma(:,:,il) + squeeze(sum(ar(:,idd,:),2));
            nb_spikes(il) = nb_spikes(il) + length(idd);
        end
         irun = irun + numel(time);
    end
    close(h)
    %calculate mean waveform for each cluster
    for i=1:M_clust
        M_template(:,:,i) = suma(:,:,i)/nb_spikes(i);
    end
    
 h = waitbar(0, 'Computing std deviation...');   
  for j=1:batch
              Progress = (j/batch);
        waitbar(Progress)
        id = rez.st(:,end-1)==j;
        time = double(rez.st(id,1));
        id_M = int32(rez.st(id,end));
        if(j==1)
            offset = 0;
        else
            offset = rez.ops.ntbuff;
        end
        ind = time - rez.ops.nt0min + offset - (rez.ops.NT - rez.ops.ntbuff)* (double(j)-1);   
        ind = int32(ind);
        dat = single(DATA(:,:,j))/rez.ops.scaleproc;
        for ite =1:M_clust
                    id_cluster = id_M==ite;  
                    time_sub = ind(id_cluster);
                    if(~isempty(time_sub))
                        inds = repmat(time_sub', nt0, 1) + repmat(int32(0:nt0-1)', 1, numel(time_sub));
                        subs = dat(inds(:),:);
                        subs = reshape(subs,[size(inds) Nchan]);
                        subs = permute(subs,[1 3 2]);
                        subs =  (subs-squeeze(M_template(:,:,ite))).*(subs-squeeze(M_template(:,:,ite)));    
                        std_sum(:,:,ite) = std_sum(:,:,ite) + sum(subs,3);
                    else
                    end
        end
  end  
  close(h)
  %calculate std waveform for each cluster
  for i=1:M_clust
      M_std_template(:,:,i) =  sqrt(std_sum(:,:,i)/(nb_spikes(i)-1));
  end
          

rez.M_template = single(M_template*rez.ops.bitmVolt*rez.ops.scaleproc);
rez.M_std_template = single(M_std_template*rez.ops.bitmVolt*rez.ops.scaleproc);
