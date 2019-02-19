function  [rez] = Waveform_T(rez,DATA)
M_T     = length(unique(rez.st(:,3)));
temp = unique(rez.st(:,3));
nt0 = rez.ops.nt0;
Nchan = rez.ops.Nchan;
%calculate the average template of the merged clusters
template = zeros(rez.ops.nt0,Nchan,M_T);
sumaT = zeros(nt0,Nchan,M_T);    %To calculate the average template
nb_spikesT = zeros(M_T,1);       %Number of spikes per cluster
batch = length(unique(rez.st(:,end-1)));
irun = 0;
h = waitbar(0, 'Computing Template...');
%calculate the mean of each template
    for j=1:batch
        Progress = (j/batch);
        waitbar(Progress)
        id = rez.st(:,end-1)==j;
        time = double(rez.st(id,1));
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
        for ik=1:M_T
            iddT = find(id_T ==temp(ik));                     
            sumaT(:,:,ik) =  sumaT(:,:,ik) + squeeze(sum(ar(:,iddT,:),2));
            nb_spikesT(ik) = nb_spikesT(ik) + length(iddT);
        end
         irun = irun + numel(time);
    end
    close(h)

    
        %calculate mean waveform for each cluster
    for i=1:M_T
        template(:,:,i) = sumaT(:,:,i)/nb_spikesT(i);
    end
    rez.template =  template*rez.ops.bitmVolt*rez.ops.scaleproc;

