<<<<<<< HEAD
function [rez] = remove_cluster(rez,SpkThr,ChC,NccThr,RefrThr)
del = [];

PETH = rez.PETH;

for i = 1:length(unique(rez.st(:,end)))
    ViolRefr = (100*mean(PETH(51,i)))/mean(PETH(:,i));    
=======
function [rez] = remove_cluster(rez,SpkThr,ChC,NccThr)
del = [];

for i = 1:length(unique(rez.st(:,end)))
>>>>>>> b40eda2ba9d8488d51e4fdab208592d3fc36a9f0
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
<<<<<<< HEAD
   
    elseif  min(rez.M_template(:,rez.Merge_cluster{i,end}(1),i) > SpkThr) || ChannelC > ChC || inx < NccThr || ViolRefr > RefrThr
%         if ~isempty(rez.st(:,end)==i)
=======
       
    elseif  min(rez.M_template(:,rez.Merge_cluster{i,end}(1),i) > SpkThr) || ChannelC > ChC || inx < NccThr
>>>>>>> b40eda2ba9d8488d51e4fdab208592d3fc36a9f0
        rm = rez.st(:,end)==i;
        rez.st(rm,:) = [];
        rez.best_channel_M(rm,:)=[];
        del = [del;rez.Merge_cluster{i,1}];
        fprintf('Cluster %d deleted\n',i);
<<<<<<< HEAD
%         end
=======
>>>>>>> b40eda2ba9d8488d51e4fdab208592d3fc36a9f0
    end
end