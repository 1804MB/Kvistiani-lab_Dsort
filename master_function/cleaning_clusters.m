function rez = cleaning_clusters(rez,DATA)
<<<<<<< HEAD
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
RefrThr = rez.ops.RefrThr;

%
%     for i =1:Ncycle
% li = 11;
i = 1;

while i < Ncycle
    NclustB =  length(unique(rez.st(:,end)));
    [rez,chan] = split_cluster(rez);
    [rez] = recompute(rez,DATA);
    NclustA = length(unique(rez.st(:,end)));
    %       [rez] = cleaning(rez,chan,pdf_Thr(i));
    [rez] = remove_cluster(rez,SpkThr(i),ChC(i),Ncc_Thr(i),RefrThr(i));  % -150 3
    [rez] = recompute(rez,DATA);
%   i = (NclustA -NclustB)*100/(NclustB);
    i = i + 1;
    disp('cleaning cycle'); disp(i);
    close all
end
%     [rez] = merging_M(rez,0.99);
%     [rez] = recompute(rez,DATA);
%     close all
=======
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
>>>>>>> b40eda2ba9d8488d51e4fdab208592d3fc36a9f0
end