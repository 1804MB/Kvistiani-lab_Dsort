function [rez] = cleaning_ncc(rez,clusterN,Th)
    id_i = find(rez.st(:,end)==clusterN);
    del_id = find(rez.st(id_i,4)<Th);
    id = id_i(del_id);
    rez.st(id,:) = [];
    rez.best_channel_M(id,:)=[];



 