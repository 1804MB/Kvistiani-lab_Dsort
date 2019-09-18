M_clust = length(unique(rez.st(:,end)));
stop=0;
prompt = {'Automatic (0) or Manual (1)','Number of cycle:'};
title = 'Input';
dims = [1 40];
definput = {'','1'};
opts.Resize = 'on';
opts.WindowStyle= 'normal';
ops.Interpreter = 'none';
answer = inputdlg(prompt,title,dims,definput,opts);
mode = str2double(answer{1});
Ncycle = str2double(answer{2});
del = [];
ChC = 10;
SpkThr = -150;
Ncc_Thr = [0 0 2 3]; 
Thr_pdf = [0.01 0.005 0 0];
if(mode==0)
    for i =1:Ncycle
        [rez,chan] = split_cluster(rez);
        [rez] = recompute(rez,DATA);
%         for iCluster = 1:length(unique(rez.st(:,end)))
%             [rez] = cleaning(rez,iCluster,chan(iCluster),Thr_pdf(i));
%         end
        [rez] = remove_cluster(rez,SpkThr,ChC,Ncc_Thr(i));  % -150 3
        [rez] = recompute(rez,DATA);
    end
             [rez] = merging_M(rez,0.9);
             [rez] = recompute(rez,DATA);
             save(fullfile(rez.ops.root,  'rezFMerged.mat'), 'rez', '-v7.3')
else
    for i= 1:M_clust
        if(stop==0)
            plot_pc(rez,i)
            wv_ccg(rez,i)
            plot_ncc(rez,i)
            plot_amp(rez,i)
            plot_time_amp(rez,i)
            prompt = {'Action: empty (skip), 1 (delete),2 (cut), 3 (clean Amplitude), 4 (clean NCC), 5 (clean Prob)','Channel','Number of clusters to cut','Threshold for Cleaning (amplitude or NCC)','Stop'};
            title = 'Input';
            i_chan = rez.Merge_cluster{i,5};
            dims = [1 80];
            definput = {'',num2str(i_chan(1)),'2','','0'};
            
            opts.Resize = 'on';
            opts.WindowStyle= 'normal';
            ops.Interpreter = 'none';
            answer = inputdlg(prompt,title,dims,definput,opts);
            option = str2double(answer{1});
            stop = str2double(answer{5});
            if(isempty(option))
                fprintf('Cluster %d kept\n',i);
            elseif(option==1)
                rm = find(rez.st(:,end)==i);
                rez.st(rm,:) = [];
                rez.best_channel_M(rm,:)=[];
                del = [del;rez.Merge_cluster{i,1}];
                fprintf('Cluster %d deleted\n',i);
            elseif(option==2)
                channel = str2double(answer{2});
                Nclust  = str2double(answer{3});
                M_clust = max(unique(rez.st(:,end)));
                [rez] = new_cluster(rez,i,channel,Nclust,M_clust);
                fprintf('Cluster %d cut\n',i);
            elseif(option==3)
                channel = str2double(answer{2});
                Th = str2double(answer{4});
                [rez] = cleaning_amp(rez,channel,i,Th);
                fprintf('Cluster %d using Amplitude\n',i);
            elseif(option==4)
                Th = str2double(answer{4});
                [rez] = cleaning_ncc(rez,i,Th);
                fprintf('Cluster %d cleaned using NCC\n',i);
            elseif(option==5)
                channel = str2num(answer{2});
                Th = str2double(answer{4});
                [rez] = cleaning(rez,i,channel,Th);
                fprintf('Cluster %d cleaned\n',i);
                fprintf('Cluster %d kept\n',i);
            else
                fprintf('Cluster %d kept\n',i);
            end
            close all;
        else
        end
        
    end
    [rez] = recompute(rez,DATA);
    save(fullfile(rez.ops.root,  'rezF3.mat'), 'rez', '-v7.3')
end


%