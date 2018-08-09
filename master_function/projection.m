function  [rez]=projection(rez,DATA)
M_clust = length(unique(rez.st(:,end)));
nt0  = rez.ops.nt0;
Nchan = rez.ops.Nchan;
batch = length(unique(rez.st(:,end-1)));
Nrank = 3;%default value for projection hard coded
 [~,W, ~, ~,~] = SVD_template(single(rez.M_template/(rez.ops.bitmVolt*rez.ops.scaleproc)), Nrank,rez.ops.Chan_criteria) ;  
%calculate the mean of each template
irunm = zeros(1,M_clust);
irun =0;
rez.PC=cell(M_clust,1);
    for j=1:batch
        id = rez.st(:,end-1)==j;
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
            idd = find(id_M ==il);   
            coefs = reshape(squeeze(W(:,il,1:Nrank))' * reshape(ar(:,idd,:), nt0, []), Nrank, numel(idd), Nchan);
            coefs = permute(coefs, [3 1 2]);
            rez.PC{il} = cat(3,rez.PC{il},coefs);
            irunm(il) = irunm(il) + numel(idd);
        end
        irun  = irun +numel(id_M);
    end

    
    %save pcfeature of each merged group by identical channels order

for i=1:M_clust
    PC = rez.PC{i}; 
    chan = rez.Merge_cluster{i,5};
    rez.PC{i} = PC(chan,:,:);    
end

 