function [rez] = cleaning_amp(rez,channel,clusterN,Th)
    id_i = find(rez.st(:,end)==clusterN);
    del_id = find(rez.best_channel_M(id_i,channel)>Th);
    id = id_i(del_id);
    rez.st(id,:) = [];
    rez.best_channel_M(id,:)=[];



 