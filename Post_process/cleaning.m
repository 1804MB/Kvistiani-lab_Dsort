function [rez] = cleaning(rez,clusterN,channel,Th)
id_i = find(rez.st(:,end)==clusterN);
i_chan=  rez.Merge_cluster{clusterN,5};
cha_id= find(i_chan==channel);
PDF = rez.PDF{clusterN,cha_id};
del_id = find(PDF<Th);
id = id_i(del_id);
rez.st(id,:) = [];


 