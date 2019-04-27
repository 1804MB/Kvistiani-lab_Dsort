% clear all
% load rez
% load TT2_1.mat
clear title
for icell_1 = [7 8 9 23] %     4  8 19  16 15] 
    for icell_2 = [7 8 9 23] % 26 22 13  30 21]     
        Cluster1 = icell_1;
        Cluster2 = icell_2;
        T1 = rez.st(rez.st(:,end)==Cluster1,2);
        T2 = rez.st(rez.st(:,end)==Cluster2,2);
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
        Raster = zeros(length(T2),100);
        for i = 1:length(T2)
            allspikes = T1 -T2(i);
            spikes = ceil(allspikes(allspikes > -50 & allspikes < 50) + 50);
            Raster(i,spikes) = 1;
        end
        figure
        bar(sum(Raster))
        title([icell_1 icell_2])
    end
end

