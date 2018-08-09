M_clust = length(unique(rez.st(:,end)));
stop=0;
del = [];
for i=1:M_clust
    if(stop==0)
    plot_pc(rez,i)
    wv_ccg(rez,i)
    plot_template_mean(rez,i)
    plot_proba(rez,i)
    prompt = {'Action: empty (skip), 1 (delete),2 (cut), 3 (split), 4 (clean)','Channel','Number of clusters to cut','Template to separate','Threshold for Cleaning','Stop'};
    title = 'Input';
    dims = [1 80];
    el = rez.Merge_cluster{i,1};
    i_chan = rez.Chan{el(1)};
    definput = {'',num2str(i_chan(1)),'2','','0.01','0'};
    opts.Resize = 'on';
    opts.WindowStyle= 'normal';
    ops.Interpreter = 'none';

    answer = inputdlg(prompt,title,dims,definput,opts);
    option = str2num(answer{1});
    stop = str2num(answer{6});
    if(isempty(option))
        fprintf('Cluster %d kept\n',i);
    elseif(option==1)
        rm = find(rez.st(:,end)==i);
        rez.st(rm,:) = [];
        del = [del;rez.Merge_cluster{i,1}];
        fprintf('Cluster %d deleted\n',i);
    elseif(option==2)
         channel = str2num(answer{2});
         Nclust  = str2num(answer{3});
        [rez] = new_cluster(rez,i,channel,Nclust,M_clust);
        fprintf('Cluster %d split\n',i); 
    elseif(option==3)
        template = str2num(answer{4});
        [rez] = separate_cluster(rez,template,M_clust);
         fprintf('Cluster %d cut\n',i); 
    elseif(option==4)
        channel = str2num(answer{2});
        Th = str2double(answer{5});       
       [rez] = cleaning(rez,i,channel,Th);
        fprintf('Cluster %d cleaned\n',i); 
    else
         fprintf('Cluster %d kept\n',i);
    end
   close all;
    else
    end
    
end

fprintf('Recalculate waveforms and projections\n');
%recalculate waveforms, projections and autocorrelogram
%Calculate waveforms
%************************************************************************************%



N_clust =length(unique(rez.st(:,3)));
cluster_left = unique(rez.st(:,3));

for i=1:N_clust
        spike_id = find(rez.st(:,3)==cluster_left(i));
        rez.st(spike_id,3) = i;   
end
[rez]=Waveform_T(rez,DATA);
[~,~,U,~] = SVD_template(rez.template,2,ops.Chan_criteria );
rez.U = U;
[rez] = Chan_cluster(rez,ops.Chan_criteria);

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
       el = rez.Merge_cluster{j,1};
       i_chan = rez.Chan{el(1)};
       rez.Merge_cluster{j,5} = i_chan;   
       rez.st(spike_id,end) =  j;
    end
%Calculate probability of each spikes to belong to its cluster
[rez] = Waveform_M(rez,DATA);
[rez]=projection(rez,DATA);
[rez]=probability(rez);
%************************************************************************************%
%Vizualize waveforms, auto-correlogram and amplitude histogram
%************************************************************************************%
[rez]=autocorrelogram(rez);
 save(fullfile(ops.root,  'rezD.mat'), 'rez', '-v7.3')
% 