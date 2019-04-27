function [rez] = cleaning(rez,i,channel,Th)

del = [];

    
    id_i = find(rez.st(:,end)==i);
    i_chan=  rez.Merge_cluster{i,5};
    cha_id= find(i_chan==channel);
    
        PDF = rez.PDF{i,cha_id};
        del_id = find(PDF<Th);
        id = id_i(del_id);
        rez.st(id,:) = [];
end

