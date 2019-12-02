%%
close all
iChan =[25 28];
del = [];
for i = 1:length(unique(rez.st(:,end)))
    
    if rez.Merge_cluster{i,end}(1) >=iChan(1) && rez.Merge_cluster{i,end}(1) <= iChan(2)
        wv_ccg(rez,i)
        %         plot_ncc(rez,i)
    end
end

%   [ NaN 14 17  42 62 3  46 18 NaN 44 NaN]
%   [ 6 47  48  53 50 43 54 51 6   16 60]
%%
close all
Cluster1 =    [10];% [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];%[1  NaN  NaN  NaN NaN  NaN  NaN  35 37 NaN 37 NaN 39 NaN 41 42 43 NaN 45 NaN NaN NaN 48 49 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ]; %[1  2  5  5  5  9  19 21 24 28 47 49 49 49];
Cluster2 =    [22];%[1   2   3   4   6   7   8   9   10   11 13   25 29  31  32  33  34  35  36  38  39  40  42  43  44  45  46  47  48  49  50  51]; %[31 2    32   3   33   4    34   5  35 6   7  38  9  10  11 12 13 14  15 46  47  17  15 19 20  51  20  51  22  52  53  24  54  25  55  26  27  57  28  58  59  30  31  29  60 ];  %[22 52 16 17 18 11 46 36 53 54 50 51 32 26];
del = [];



for i = 1:length(Cluster2)
    inx1 = find(rez.st(:,end)==Cluster2(i));
    if ~isnan(Cluster1(i))
        rez.st(inx1,end) = Cluster1(i);
    end
    rm = find(rez.st(:,end)==Cluster2(i));
    rez.st(rm,:) = [];
    del = [del;rez.Merge_cluster{Cluster2(i),1}];
end

    %  [rez] = merging_M(rez,0.9);
    [rez] = recompute(rez,DATA);
    save(fullfile(rez.ops.root,  'rezFinal.mat'), 'rez', '-v7.3')
    

% for i = 1:length(Cluster2)
%     id = Cluster1(i);
%     inx1 = find(rez.st(:,end)==Cluster2(i));
%     if ~isnan(Cluster1(i))
%         rez.st(inx1,end) = Cluster1(i);
%     else
%         id2 = Cluster2(i);
%         
%         [~, inx1, inx2] = union(rez.st(id,2),rez.st(id2,2));
%         id3 = sort([id(inx1)' id2(inx2)']);
%         rez.st(id3,end) = Cluster1(i);
%         id = id3';
%         clear inx1 inx2
%         if ~isempty(rez.st(:,end) == Cluster2(i))
%             rm = rez.st(:,end) == Cluster2(i);
%             rez.st(rm,:) = [];
%             id = find(rez.st(:,end)==Cluster1(i));
%         end
%     end
% end
    

