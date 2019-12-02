%specify clusters to load
clear inx
Cluster1 = 3;
Cluster2 = 19;

% inx1 = find(rez.st(:,end)==Cluster1); inx1 = inx1(inxT1);
% inx2 = find(rez.st(:,end)==Cluster1); inx2 = inx2(inxTO1);
% inx{1} = inx1;
% inx{2} = inx2;

inx1 =  find(rez.st(:,4) >= 0);
inx2 =  find(rez.st(:,4) >= 0);
inx{1} = intersect(find(rez.st(:,end) == Cluster1),inx1);
inx{2} = intersect(find(rez.st(:,end) == Cluster2),inx2);



clear inx1 inx2
% inx{2} = find(rez.st(:,end)== Cluster2 & rez.st(:,4) > 0.9 & rez.st(:,5) < -100);
% inx{1} = find(rez.st(:,end) == Cluster1 & rez.st(:,4) > 0.9 & rez.st(:,5) < -100);
Sorting = 'DSort'; %'MClust' 'DSort'
% specify channels for the tetrode

[d, timestamps, info] = load_open_ephys_data(sprintf('100_CH%d.continuous',1));
data = zeros(32,round(length(d)/10));
[b1, a1] = butter(3, [rez.ops.fshigh/rez.ops.fs,rez.ops.slow/rez.ops.fs]*2, 'bandpass');

for i = [5:8] %1:32
    [d, timestamps, info] = load_open_ephys_data(sprintf('100_CH%d.continuous',i));
    d = filtfilt(b1, a1, d);
    data(i,:) = d(1:size(data,2),1);
    clear d;
end
% dataW = data'*rez.Wrot;


dataW = data';
% for i = 1:32
%     dataW(i,:) = dataW(i,:) - mean(dataW);
% end
% dataW = dataW';
% dataW = data';
SelectCh = [5:8];

% data1 = dataW(:,29);
% data2 = dataW(:,30);
% data3 = dataW(:,31);
% data4 = dataW(:,32);



% plot mean and std of the wv
YLim = [-800 200];
XLim = [ 0 61];

% switch Sorting
%     case 'DSort'
%         for iInx = 1:length(inx)
%             inxS  = inx{iInx};
%             D = zeros(length(inxS),4,length(XLim(1):XLim(2)-1));
%             sel_inx = 1:length(inxS);
%             for i = 2:length(sel_inx)-1
%                 for k = 1:length(SelectCh)
%                     D(i,k,:)= dataW(rez.st(inxS(sel_inx(i)),1)-20:rez.st(inxS(sel_inx(i)),1)+40,SelectCh(k));
%                 end
%             end
%             
%             figure
%             for k = 1:4
%                 subplot(2,2,k)
%                 %                 stdshade_sorting(mean(squeeze(D(:,k,:))),std(squeeze(D(:,k,:))),'Color','b')
%                 errorbar(mean(squeeze(D(:,k,:))), std(squeeze(D(:,k,:)))/2 ,'Color','b')
%                 %                   plot([1:61],mean(squeeze(D(:,k,:))),'Color','b')
%                 
%                 axis([XLim YLim])
%             end
%         end
        
        
%         for i = 1:4
%             for k = 1:2
%                 figure;
%                 plot(dataW(rez.st(inx{1},1),17),dataW(rez.st(inx{1}),20),'*')
                
                
                
                
                
                sel_inx = 1:100; %26:27;
                inx1 = inx{1};
                inx2 = inx{2};
                
%                 inx1 = DiffT1{1}*30;
%                 inx2 = DiffT2{1}*30;
%                 inx1 = OverlapT1*30;
                SpkWnd = 120;
                Offset = SpkWnd + 30;
                Factor = 500;
%                 figure
                clear title
                for i = 2:length(sel_inx)
%                      h = figure;
%                     movegui(h,'south');
%                     subplot(1,2,1)
                    f = 0;
                    for ich = 1:4
                        figure(2)
                        subplot(4,1,ich) 
                       plot(dataW(rez.st(inx1(sel_inx(i),1))-SpkWnd:rez.st(inx1(sel_inx(i),1))+SpkWnd,SelectCh(ich)) - dataW(rez.st(inx1(sel_inx(i),1))-Offset,SelectCh(ich)) + Factor*f, 'b')
                       hold on 
                       plot(dataW(rez.st(inx2(sel_inx(i),1))-SpkWnd:rez.st(inx2(sel_inx(i),1))+SpkWnd,SelectCh(ich)) - dataW(rez.st(inx2(sel_inx(i),1))-Offset,SelectCh(ich)) + Factor*f,'r')

                       %plot(eval([sprintf('data%d(rez.st(inx1(sel_inx(i)),1)-SpkWnd:rez.st(inx1(sel_inx(i)),1)+SpkWnd)',ich) sprintf('-data%d(rez.st(inx1(sel_inx(i)),1)-SpkWnd)',ich)]),'Linewidth',2)
%                         plot(dataW(inx1(i) -SpkWnd:inx1(i) + SpkWnd,SelectCh(ich)) - dataW(inx1(i)- Offset,SelectCh(ich)) + Factor*f, 'b')

                         hold on
                         axis([ 0 240 -250 100])

%                         title(['Prob =' num2str(rez.st(inx1(i),4)),'Ampl = ' num2str(rez.st(inx1(i),5))])
%                         f = f + 1;
                    end
                     f = 0;
%                        subplot(1,2,2)
%                  h = figure;
%                  movegui(h,'north');
%                      for ich = 1:4
%                       
%                         
%                          subplot(4,1,ich)
%                          plot(dataW(rez.st(inx2(sel_inx(i),1))-SpkWnd:rez.st(inx2(sel_inx(i),1))+SpkWnd,SelectCh(ich)) - dataW(rez.st(inx2(sel_inx(i),1))-Offset,SelectCh(ich)) + Factor*f,'r')
% %                         %plot(eval([sprintf('data%d(rez.st(inx2(sel_inx(i)),1)-SpkWnd:rez.st(inx2(sel_inx(i)),1)+SpkWnd)',ich) sprintf('-data%d(rez.st(inx2(sel_inx(i)),1)-SpkWnd)',ich)]),'Linewidth',2)                        
% %                         plot(dataW(inx2(i) -SpkWnd:inx2(i) + SpkWnd,SelectCh(ich)) - dataW(inx2(i)- Offset,SelectCh(ich)) + Factor*f, 'b')
%                         hold on
% %                          title(['Prob =' num2str(rez.st(inx2(i),4)),'Ampl = ' num2str(rez.st(inx2(i),5))])
% %                          f = f + 1;
%                            axis([ 0 600 -800 400])
%                      end
                end
                
%                 case 'MClust'
                     load('TT1_02.mat')
                    % tSpikes = TSAligned;
%                   inx{1}= double(tSpikes*3) - timestamps(1)*30000;
                    inx{1}=  tSpikes*30000 - timestamps(1)*30000;
                    
                    % load('TT2_2.mat')
                    load('TT1_03.mat');
                   inx{2} =  tSpikes*30000 - timestamps(1)*30000;
%                     inx{2} = double(tSpikes*3) - timestamps(1)*30000;
%                     % inx{2} = tSpikes; %*3;% double(tSpikes)*3 - timestamps(1)*30000;
%                     
%                     sel_inx = 1:length(inxS);
%                     for iInx = 1:length(inx)
%                         inxS  = inx{iInx};
%                         sel_inx = 1:length(inxS);
%                         D1 = zeros(length(sel_inx),61);
%                         D2 = zeros(length(sel_inx),61);
%                         D3 = zeros(length(sel_inx),61);
%                         D4 = zeros(length(sel_inx),61);
%                         
%                         for i = 2:length(sel_inx)-1
%                             D1(i,:) = data1(rez.st(inxS(sel_inx(i)),1)-20:rez.st(inxS(sel_inx(i)),1)+40);
%                             D2(i,:) = data2(rez.st(inxS(sel_inx(i)),1)-20:rez.st(inxS(sel_inx(i)),1)+40);
%                             D3(i,:) = data3(rez.st(inxS(sel_inx(i)),1)-20:rez.st(inxS(sel_inx(i)),1)+40);
%                             D4(i,:) = data4(rez.st(inxS(sel_inx(i)),1)-20:rez.st(inxS(sel_inx(i)),1)+40);
%                         end
%                         figure
%                         subplot(2,2,1)
%                         errorbar(mean(D1), std(D1)/2)
%                         axis([XLim YLim])
%                         subplot(2,2,2)
%                         errorbar(mean(D2), std(D2)/2,'Color','r')
%                         axis([XLim YLim])
%                         subplot(2,2,3)
%                         errorbar(mean(D3), std(D3)/2,'Color','g')
%                         axis([XLim YLim])
%                         subplot(2,2,4)
%                         errorbar(mean(D4), std(D4)/2,'Color','k')
%                         axis([XLim YLim])
%                     end
%                     
%                     %
                    sel_inx = 1:100;
                    inx1 = round(inx{1});
                    inx2 = round(inx{2});
                    SpkWnd = 120;
                    Factor = 500;
%                     figure
                    clear title
                    for i = 2:length(sel_inx)
%                         h = figure;
%                         movegui(h,'north');
%                         subplot(1,2,1)
                        f = 0;
             
                        for ich = 1:4  
                             figure(1)
                            subplot(4,1,ich)                         
                            plot(dataW(inx1(sel_inx(i),1)-SpkWnd:inx1(sel_inx(i),1)+SpkWnd,SelectCh(ich)) - dataW(inx1(sel_inx(i),1)-SpkWnd,SelectCh(ich)) + Factor*f, 'b')
                            hold on
                            plot(dataW(inx2(sel_inx(i),1)-SpkWnd:inx2(sel_inx(i),1)+SpkWnd,SelectCh(ich)) - dataW(inx2(sel_inx(i),1)-SpkWnd,SelectCh(ich)) + Factor*f, 'r')
                            hold on
                            axis([ 0 240 -250 100])
                           

%                             f = f + 1;
                        end
%                         f = 0;
%                         subplot(1,2,2)
  
%                         for ich = 1:4
%                             figure(2)
%                             subplot(4,1,ich) 
%                             plot(dataW(inx2(sel_inx(i),1)-SpkWnd:inx2(sel_inx(i),1)+SpkWnd,SelectCh(ich)) - dataW(inx2(sel_inx(i),1)-SpkWnd,SelectCh(ich)) + Factor*f, 'r')
%                             hold on
% %                              f = f + 1;
%                              axis([ 0 240 -500 500])
%                         end
                    end
%             end
