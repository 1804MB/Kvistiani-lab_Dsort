clear all
session= dir;
session = session(3:end);
RootDir = cd ;

for isess = 11:length(session)
%     try
        sessionpath = [RootDir,'\',session(isess).name];
        cd(sessionpath)
        if exist('rezDMMerged.mat') || exist('rezDMerged.mat') 
            continue
        end
        load rezI
        [rez, DATA] = preprocessData(rez.ops);
        load rezI
        M_clust = length(unique(rez.st(:,end)));
        dims = [1 40];
        opts.Resize = 'on';
        opts.WindowStyle= 'normal';
        ops.Interpreter = 'none';
        Ncycle = 8;
        del = [];
        ChC = 10;
        SpkThr = [-100 -150 -150 -150 -150 -150 -150 -150];
        Ncc_Thr = [0 0 2 3 3 3 3 3];
        Thr_pdf = [0.05 0.02 0.01 0];
        for i =1:Ncycle
            [rez,chan] = split_cluster(rez);
            [rez] = recompute(rez,DATA);
            %         for iCluster = 1:length(unique(rez.st(:,end)))
            %             [rez] = cleaning(rez,iCluster,chan(iCluster),Thr_pdf(i));
            %         end
            [rez] = remove_cluster(rez,SpkThr(i),ChC,Ncc_Thr(i));  % -150 3
            [rez] = recompute(rez,DATA);
            close all
        end
        [rez] = merging_M(rez,0.99);
        [rez] = recompute(rez,DATA);
        save(fullfile(sessionpath,  'rezF2.mat'), 'rez', '-v7.3')
        cd(RootDir)
        close all
        disp(isess)
%     catch
%         continue
%     end
end

