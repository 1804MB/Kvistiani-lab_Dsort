function [rez] = remove_cluster(rez,SpkThr,ChC,NccThr)
del = [];

for i = 1:length(unique(rez.st(:,end)))
    try
        ncc_score = rez.st(rez.st(:,end)==i,4);
        [Ncc,~] = histcounts(ncc_score,min(ncc_score):0.01:max(ncc_score));
        NumChannel = length(rez.Merge_cluster{i,end});
    catch; Ncc = zeros(1,10);
    end
    [~, inx] = max(Ncc);
        
    ChannelC = length(rez.Merge_cluster{i});
    
    if  isempty( rez.Merge_cluster{i,end})
        rm = rez.st(:,end)==i;
        rez.st(rm,:) = [];
        rez.best_channel_M(rm,:)=[];
        del = [del;rez.Merge_cluster{i,1}];
        fprintf('Cluster %d deleted\n',i);
        
%     elseif min(rez.M_template(:,rez.Merge_cluster{i,end}(1),i)) < - 400 && ChannelC < ChC
%             continue
       
    elseif  min(rez.M_template(:,rez.Merge_cluster{i,end}(1),i) > SpkThr) || ChannelC > ChC || inx < NccThr
        rm = rez.st(:,end)==i;
        rez.st(rm,:) = [];
        rez.best_channel_M(rm,:)=[];
        del = [del;rez.Merge_cluster{i,1}];
        fprintf('Cluster %d deleted\n',i);
    end
end