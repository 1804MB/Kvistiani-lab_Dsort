function [Chan] = Chan_cluster2(U,N_clust,criteria)
U = U(:,:,1);
Chan = cell(N_clust,1);
for i =1:N_clust
     U(:,i) = U(:,i)/max(U(:,i));
     id = find(U(:,i)>=criteria);
     [~,chan_order] = sort(U(id,i),'descend');
     Chan{i} =  id(chan_order);
end

   