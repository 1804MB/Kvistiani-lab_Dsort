function [rez] = Chan_template(rez,criteria)

N_clust = length(unique(rez.st(:,3)));
U = rez.U(:,:,1);
Chan = cell(N_clust,1);
Chan_score = cell(N_clust,1);
for i =1:N_clust
     U(:,i) = U(:,i)/max(U(:,i));
     id = find(U(:,i)>=criteria);
     [~,chan_order] = sort(U(id,i),'descend');
     Chan{i} =  id(chan_order);
     Chan_score{i} = U(id(chan_order),i);
end
rez.Chan =  Chan;
rez.Chan_score = Chan_score;
   