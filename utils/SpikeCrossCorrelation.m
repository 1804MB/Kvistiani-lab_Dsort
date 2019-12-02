% clear all
% load rez
% load TT2_1.mat
<<<<<<< HEAD
%   [d, timestamps, info] = load_open_ephys_data(sprintf('100_CH1.continuous'));
clear title
k = 1;
Toler = 1;
for icell_1 = [10] %     4  8 19  16 15] 
     for icell_2 = [22] % 26 22 13  30 21]     
=======
clear title
for icell_1 = [7 8 9 23] %     4  8 19  16 15] 
    for icell_2 = [7 8 9 23] % 26 22 13  30 21]     
>>>>>>> b40eda2ba9d8488d51e4fdab208592d3fc36a9f0
        Cluster1 = icell_1;
        Cluster2 = icell_2;
        T1 = rez.st(rez.st(:,end)==Cluster1,2);
        T2 = rez.st(rez.st(:,end)==Cluster2,2);
<<<<<<< HEAD
        
%         load('TT5_01.mat'); 
% %         T2 = double(tSpikes*10000000) - timestamps(1)*1000;
%         T2 = tSpikes*1000 - timestamps(1)*1000;
%         T2 = round(T2);
% 
%         load('TT5_02.mat'); 
% %         T2 = double(tSpikes*10000000) - timestamps(1)*1000;
%         T1 = tSpikes*1000 - timestamps(1)*1000;
%         T1 = round(T1);

%         T4 = round((tSpikes*1000 - timestamps(1)*1000)/10);
%         T4 = T4*10;
%          T2 = T4;
=======
>>>>>>> b40eda2ba9d8488d51e4fdab208592d3fc36a9f0
%         T1 = rez.st(inx2,2) ;
%         T2 = rez.st(inx1,2) ;

%         T1 = rez.st(rez.st(:,end)==Cluster1,2);
%         T2 = (double(tSpikes))/10 -timestamps(1)*1000;

%         TS1 = load('TT2_2');
%         TS2 = load('TT3_1');
%         T1 = TS1.tSpikes/10;
%         T2 = TS2.tSpikes/10;
%         T1 = double(T1);
%         T2 = double(T2);
<<<<<<< HEAD
OverlapT1 = [];
OverlapT2 = cell(1,length(T2));

=======
>>>>>>> b40eda2ba9d8488d51e4fdab208592d3fc36a9f0
        Raster = zeros(length(T2),100);
        for i = 1:length(T2)
            allspikes = T1 -T2(i);
            spikes = ceil(allspikes(allspikes > -50 & allspikes < 50) + 50);
            Raster(i,spikes) = 1;
<<<<<<< HEAD
            % find overlapping spikes

             OverlapT1 = [OverlapT1 T1(abs(T1-T2(i)) <= Toler)'];
%            OverlapT2{i} = T2(find(abs(T1-T2(i) <= 1)));
        end
        
=======
        end
>>>>>>> b40eda2ba9d8488d51e4fdab208592d3fc36a9f0
        figure
        bar(sum(Raster))
        title([icell_1 icell_2])
    end
<<<<<<< HEAD
     [DiffT1, inxT1] = setdiff(T1,OverlapT1);
     [DiffT2, inxT2] = setdiff(T2,OverlapT1);
     k = k + 1;
=======
>>>>>>> b40eda2ba9d8488d51e4fdab208592d3fc36a9f0
end

