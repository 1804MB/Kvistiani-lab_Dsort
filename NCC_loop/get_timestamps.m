function [time] = get_timestamps(beta,nt0)        

blocks = 0;
indexs = 1;
while(blocks==0)
    over =0;
    toto = beta(indexs:indexs+nt0-1);
    [ac,index] =max(toto);
    index = index -1;
    if(ac==0)
        indexs = indexs + nt0;
    else
        ind_o = index;
        while(over==0)
            i0 = indexs+ind_o;
            toto = beta(i0:i0+nt0-1);
            [~,ind_n] =max(toto);
            if(ind_n==1)   
                ind_o = ind_o + ind_n -1;
                over =1;
            else
                ind_o = ind_o + ind_n - 1;
            end
         end                        
         i0 = indexs+ind_o;
         beta(indexs:i0-1) =0;
         beta(i0+1:i0+nt0-1) =0;
         indexs = i0 + nt0; 
     end
     if(indexs+nt0>length(beta))
        blocks = 1;
     end
end
time = find(beta>0);