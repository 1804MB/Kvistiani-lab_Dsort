function [rez] = new_cluster(rez,clusterN,channel,Nclust,M_clust)
id_i = find(rez.st(:,end)==clusterN);
PC = rez.PC{clusterN};
i_chan=  rez.Merge_cluster{clusterN,5};
cha_id= find(i_chan==channel);
pcs = squeeze(PC(cha_id,:,:))';
M_template = max(unique(rez.st(:,3)))+1;
% Nclust = 2;
GMModel = fitgmdist(pcs,Nclust ,'Replicates',100);
[clusterX,nlogl,P,logpdf,M] = cluster(GMModel ,pcs);
% clusterX = kmeans(pcs(:,:),2)
RC = {'b.','r.','g.','y.','m.','c.','k.','bo'};
figure;
for i=1:Nclust  
    id = find(clusterX==i);
    plot3(pcs(clusterX==i,1),pcs(clusterX==i,2),pcs(clusterX==i,3),RC{i},'MarkerSize',2)
    hold on
    if i > 1
    id_s = find(clusterX==i);
    id_s = id_i(id_s);
    M_clust = M_clust + 1;
    rez.st(id_s,end) = M_clust;   
    rez.st(id_s,3) = M_template;  
    end
end

prompt = 'press enter ';
x = input(prompt);