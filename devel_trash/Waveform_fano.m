function  [rez]=Waveform_fano(rez,N_clust,DATA)
NchanTOT = rez.ops.NchanTOT;

%calculate the average template of the merged clusters
template = zeros(NchanTOT,rez.ops.nt0,N_clust);
std_template = zeros(NchanTOT,rez.ops.nt0,N_clust);
Fano = zeros(NchanTOT,N_clust);

for i=1:N_clust
    id = find(rez.st3(:,3)==i);
    batch = rez.st3(id,6);
    WAVEGPU = NaN(rez.ops.nt0,NchanTOT,length(id));
    for j=1:length(batch)
        dat = DATA(:,:,batch(j));
        time = double(rez.st(id(j),1));
        if(batch(j)==1)
            offset = 0;
        else
            offset = rez.ops.ntbuff;
        end
        ind = time - rez.ops.nt0min + offset - (rez.ops.NT - rez.ops.ntbuff)* (double(batch(j))-1);   
        ind = int32(ind);
        dataRAW = gpuArray(dat);
        dataRAW = single(dataRAW);
        dataRAW = single(dataRAW) / rez.ops.scaleproc;
        datSp       = dataRAW(ind:ind+int32(rez.ops.nt0)-1, :);
        WAVEGPU(:,:,j) = gather(datSp);
    end
    WAVEGPU = permute(WAVEGPU,[2,1,3]);
    template(:,:,i) = mean(WAVEGPU,3);
    std_template(:,:,i) = std(WAVEGPU,0,3);
    [a,~] = max(-squeeze(template(:,:,i)));
    [~,c] = max(a);
    Fano(:,i) = std_template(:,c,i).*std_template(:,c,i)./template(:,c,i);
end



rez.template = template;
rez.std_template = std_template;
rez.Fano = -Fano;
% Templates
%  for i = 1:N_clust
%     figure;   
%      for k = 1:NchanTOT
%         amean = [squeeze(template(k, :, i)); squeeze(template(k, :, i))];
%         astd = [squeeze(std_template(k, :, i)); squeeze(std_template(k, :, i))];
%         
%        subplot(ceil(NchanTOT/3),3,k)
%        stdshade_sorting(amean, astd, 0.3, 'b')
%        title([ 'Channel: ', num2str(k), 'Fano: ', num2str(Fano(k,i))])
%        xlim([0 60]);
%      end 
%      [ax,h3]=suplabel(['Cluster: ', num2str(i)],'t'); 
%      set(h3,'FontSize',22) 
%           
% end