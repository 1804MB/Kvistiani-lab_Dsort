function  [rez]=check_template(DATA,st,T,ops)
M_clust = length(unique(st(:,3)));
nt0  = ops.nt0;
Nchan = 
batch = length(unique(rez.st(:,end-1)));
Nrank = 3;%default value for projection hard coded
[~,W, ~, ~,~] = SVD_template(T, Nrank,rez.ops.Chan_criteria) ;  
%calculate the mean of each template
irunm = zeros(1,M_clust);
irun =0;
rez.PC=cell(M_clust,1);
rez.distPC=cell(M_clust,1);
    for j=1:batch
        id = rez.st(:,end)==j;
        time = double(rez.st(id,1));
        id_M = int32(rez.st(id,end));
        if(j==1)
            offset = 0;
        else
            offset = rez.ops.ntbuff;
        end
        ind = time - rez.ops.nt0min + offset - (rez.ops.NT - rez.ops.ntbuff)* (double(j)-1);   
        ind = int32(ind);
        inds = repmat(ind', nt0, 1) + repmat(int32(1:nt0)', 1, numel(ind));
        dat = single(DATA(:,:,j))/rez.ops.scaleproc;
        ar = dat(inds(:),:);        
        ar = reshape(ar,[size(inds) Nchan]);
        for il=1:M_clust
            channel = rez.Merge_cluster{il,5};
            idd = find(id_M ==il);   
            coefs = reshape(squeeze(W(:,il,1:Nrank))' * reshape(ar(:,idd,channel), nt0, []), Nrank, numel(idd), length(channel));
            coefs = permute(coefs, [3 1 2]);
            ard = permute(ar(:,idd,channel),[2,1,3]);
            ard = (reshape(ard,length(idd),length(channel)*rez.ops.nt0));   
            wave =  rez.M_template(:,channel,il)/(rez.ops.bitmVolt*rez.ops.scaleproc);
            wave =  reshape(wave,[rez.ops.nt0*length(channel) 1 ]);
            rez.distPC{il} =  cat(1,rez.distPC{il},ard * wave);
            
            rez.PC{il} = cat(3,rez.PC{il},coefs);
            irunm(il) = irunm(il) + numel(idd);
        end
        irun  = irun +numel(id_M);
    end

    
%     %save pcfeature of each merged group by identical channels order
% 
for i=1:M_clust
    PC = rez.distPC{i}; 
    [N,edge]=histcounts(PC);
    figure;
    bar(edge(1:end-1),N);
    title(['Final distribution, cluster=',num2str(i)]) 
end

 