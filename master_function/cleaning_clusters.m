function rez = cleaning_clusters(rez,DATA)
    M_clust = length(unique(rez.st(:,end)));
    dims = [1 40]; 
    del = [];
    opts.Resize = 'on';
    opts.WindowStyle= 'normal';
    ops.Interpreter = 'none';
    
    Ncycle = rez.ops.Ncycle;
    ChC = rez.ops.ChC;         % set to 10 
    SpkThr = rez.ops.SpkThr;   % -150;
    Ncc_Thr = rez.ops.Ncc_Thr; % [0 0 2 3];
    pdf_Thr = rez.ops.pdf_Thr;     % [0.05 0.02 0.01 0];

    for i =1:Ncycle
        [rez,chan] = split_cluster(rez);
        [rez] = recompute(rez,DATA);
%       [rez] = cleaning(rez,chan,pdf_Thr(i));
        [rez] = remove_cluster(rez,SpkThr(i),ChC,Ncc_Thr(i));  % -150 3
        [rez] = recompute(rez,DATA);
        close all
    end
    [rez] = merging_M(rez,0.99);
    [rez] = recompute(rez,DATA);
    close all
end