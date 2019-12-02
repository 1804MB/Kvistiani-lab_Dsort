%specify clusters to load
% load rezFinal
% ScalD = 10;
  SelectCh =[1:32];
% [d, timestamps, info] = load_open_ephys_data(sprintf('100_CH%d.continuous',1));
% [b1, a1] = butter(3, [rez.ops.fshigh/rez.ops.fs,rez.ops.slow/rez.ops.fs]*2, 'bandpass');
% inx1 = round(1:length(d)/30);
% inx2 = round([length(d)/2:(length(d)/2 + length(d)/30)]);
% inx3 = round([(length(d) - length(d)/30):length(d)]);
% d = d([inx1 inx2 inx3]);
% d = filtfilt(b1, a1, d);
% data = zeros(32,length(d));
% data(1,:) = d(1:size(data,2),1);
% 
% for i = 2 : 32 %length(SelectCh)
% [d, timestamps, info] = load_open_ephys_data(sprintf('100_CH%d.continuous',i));
% d = d([inx1 inx2 inx3]);
% d = filtfilt(b1, a1, d);
% data(i,:) = d(1:size(data,2),1);
% end
% 
% dataW = data'*rez.Wrot;
nt0 = rez.ops.nt0;
inx = find(rez.st(:,end) == 1);
WV = zeros(60,32,length(inx)); 
SpkWnd = [19 40];
Bit2uVolt = rez.ops.bitmVolt*1000;
MaxSpk = 100;

for ich = 1:32
    for i = 1:MaxSpk
        if intersect([inx1 inx2 inx3],inx(i))
        WV(:,ich,i) = dataW(rez.st(inx(i),1)-SpkWnd(1):rez.st(inx(i),1)+SpkWnd(2),SelectCh(ich)) - dataW(rez.st(inx(i),1)- SpkWnd(1),SelectCh(ich));
        WV(:,ich,i) = WV(:,ich,i)*Bit2uVolt;
        else
        WV(:,ich,i) = [];  
        end
    end
end
WV(:,:,MaxSpk + 1:end) = [];

for i = 1:100
as = WV(:,1:4,i); 
as = reshape(as,[],4*60); %rez.ops.NchanTOT*nt0);
                    
bs = rez.M_template(:,1:4,1); % use all channels for merging
bs = reshape(bs,[],4*60);%rez.ops.NchanTOT*nt0);
score(i) = xcorr(as,bs,0,'coeff');
end

figure(1)

for ich = 1:4
    for i = 1:100
        subplot(2,2,ich)
        plot(squeeze(WV(ich,i,:)))
        hold on
    end
%     figure(2)
%     subplot(2,2,ich)
%     plot(mean(squeeze(WV(ich,:,:)))); hold on;
end


