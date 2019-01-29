M_clust = length(unique(rez.st(:,end)));
stop=0;
del = [];
for i=1:M_clust
    if(stop==0)
        plot_pc(rez,i)
        wv_ccg(rez,i)
        plot_ncc(rez,i)
        plot_amp(rez,i)
        prompt = {'Action: empty (skip), 1 (delete),2 (cut), 3 (clean Amplitude), 4 (clean NCC)','Channel','Number of clusters to cut','Threshold for Cleaning (amplitude or NCC)','Stop'};
        title = 'Input';
        i_chan = rez.Chan{i};
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
            [rez] = new_cluster(rez,i,channel,Nclust,M_clust);
            fprintf('Cluster %d split\n',i); 
        elseif(option==3)
            channel = str2double(answer{2});
            Th = str2double(answer{4});  
            [rez] = cleaning_amp(rez,channel,i,Th);
            fprintf('Cluster %d using Amplitude\n',i); 
        elseif(option==4)
            Th = str2double(answer{4});       
           [rez] = cleaning_ncc(rez,i,Th);
            fprintf('Cluster %d cleaned using NCC\n',i); 
        else
            fprintf('Cluster %d kept\n',i);
        end
    close all;
    else
    end
    
end
%recalculate waveforms, projections and autocorrelogram
%Calculate waveforms
cluster = unique(rez.st(:,end));
M_clust = length(unique(rez.st(:,end)));
for j=1:M_clust
       spike_id = find(rez.st(:,end)==cluster(j));  
       rez.st(spike_id,end) =  j;
end
%************************************************************************************%
fprintf('Recalculate waveforms\n');
[rez] = Waveform_M(rez,DATA);
[~,~,U,~] = SVD_template(rez.M_template,rez.ops.Nrank,rez.ops.Chan_criteria );
rez.U = U;
[rez] = Chan_cluster_M(rez,rez.ops.Chan_criteria);

cluster = unique(rez.st(:,end));
M_clust = length(unique(rez.st(:,end)));
fprintf('Number of clusters %d\n',M_clust); 
rez.Merge_cluster = cell(M_clust,5);
    for j=1:M_clust
       spike_id = find(rez.st(:,end)==cluster(j));
       rez.Merge_cluster{j,1} = unique(rez.st(spike_id,3));
       rez.Merge_cluster{j,2} = length(spike_id);
       rez.Merge_cluster{j,3} = j;
       rez.Merge_cluster{j,4} = 1/(length(spike_id)/(rez.ops.NT*rez.temp.Nbatch/(rez.ops.fs))); 
       el = rez.Merge_cluster{j,3};
       i_chan = rez.Chan{el(1)};
       rez.Merge_cluster{j,5} = i_chan;   
       rez.st(spike_id,end) =  j;
    end
fprintf('Recalculate projections\n');
%Calculate probability of each spikes to belong to its cluster
[rez] = projection(rez,DATA);
%************************************************************************************%
%Vizualize waveforms, auto-correlogram and amplitude histogram
%************************************************************************************%
[rez]=autocorrelogram(rez);
 save(fullfile(rez.ops.root,  'rezDs.mat'), 'rez', '-v7.3')
% 