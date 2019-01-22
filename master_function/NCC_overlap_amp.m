%*************************************************************************%
%*************************************************************************%
%fitTemplates: Detect spike and calculate the templates
%Simplified and commented by Madeny belkhiri 30/01/18.
%*************************************************************************%
%*************************************************************************%

function [st3,all] = NCC_overlap_amp(DATA, dWU,ops)

[NT, Nchan ,Nbatch] = size(DATA);   

nt0     = ops.nt0;       % Number of points in sampling rate for a waveform
nt0min  = ops.nt0min;    % Strange constant added at the timesstamps why [MB?]
Nfilt 	= ops.Nfilt;                 % Number of Clusters provided by user via config file
Chan_criteria = ops.Chan_criteria;

last = 0;
st3 = [];

Nrank = get_rank(dWU,ops.variance);
ops.Nrank = Nrank;      % Rank to use for the SVD decomposition given by user via config file
rez.ops.Nrank = Nrank;
fprintf('Rank used for calculation:%d\n', Nrank)
dWU0 = dWU;
%Start the loop to detect spikes and optimize the template,
criteria_NCC = ops.criteria_NCC;
suma0 = zeros(nt0,Nchan,Nfilt);%To calculate the average template
nb_spikes0 =  zeros(Nfilt,1) ;%number of spikes
nbsp_tot = 0;%count cumulated number of spikes
stop =0;   
iteration =1;
all_channel = [];
irun = 0;
while(stop==0 || criteria_NCC>=1)
   if iteration>1
       last =1;
   end
   batch_score = zeros(Nbatch,NT);
   batch_id    = zeros(Nbatch,NT);
   nbsp = 0;
   %Decompose the template dWU (called K in the paper) in two matrix W*U using SVD.
   % Using svd one get dWU = W*D*U, in fact output W is W*D, D being the
   % diagonal matrix containing the eigenvalues
   %Chan_criteria is used to disregard irrelevant channels
  
   [dWU,W, U, Weight, Topchan] = SVD_template(dWU, Nrank,Chan_criteria );

    Weight = mean(Weight,2);
    W = reshape(W,nt0,Nrank*Nfilt);
   %Some plot of the template, number of spikes and mask
  if iteration==1
    figure;
    subplot(1,2,1)
    plot(W(:,:,1))
    xlim([0 nt0-1]);
    title('W(:,:,1)')            
    subplot(1,2,2)
    imagesc(U(:,:,1))
    title('U(:,:,1)')     
    drawnow;
  end
   Params = int32([NT Nfilt*Nrank nt0]); %Used for the convolution
        for i=1:Nbatch 
               time = [];
               id = [];
            %Take a batch of the whitened and filtered data
             dat = single(DATA(:,:,i))/ops.scaleproc;  
            %project in the spatial component base
             data = dat*U(:,:);

             % Call CUDA to execute convolution
             [tconv] = mexconv(Params,data,gpuArray(W));  
             tconv = gather(tconv);
             tconv (isnan(tconv ))=0;
             tclus= zeros(NT,Nfilt,'single');

             %Calculate the normalized cross-correlation using a weights of
             %eigenvalue of each components
             for ki=1:Nfilt
                for irank =1:Nrank
                     tclus(:,ki) =  tclus(:,ki) + tconv(:,ki+(irank-1)*Nfilt).*Weight(irank);
                end
             end
             tclus(end-ops.nt0:end,:)=0; % the first nt0 points ar useless since we do a NCC
             %0 all the point below the score threshold
             %find the higest Normalized-Cross-Correlation(NCC) value for each
             %sampling point
             [best,it]= max(tclus,[],2);
             batch_score(i,:) = best;
             batch_id(i,:)    = it;
             %Apply the criteria
             best(best<criteria_NCC)=0;
             %Search for the spikes timestamps if score goes above threshold
             %Find the highest NCC value every nt0 points (should correspond to 1 ms)
             [time] =  get_timestamps(best,nt0);
             %Get the corresponding Cluster id
             id = int32(it(time));
             %number of spikes
             nbsp = nbsp + length(time);
             time  = int32(time);
             Top_chan = Topchan(id);

            %Save the detected spikes and projections
            if(last==1)
                if ~isempty(time)
                        inds  = repmat(time', nt0, 1) + repmat(int32(1:nt0)', 1, numel(time));
                        try  datSp  = dat(inds(:), :);
                        catch
                        datSp = dat(inds(:), :);
                        end
                        datSp = reshape(datSp, [size(inds) Nchan]);
                        Best_channel = diag(squeeze(datSp(nt0min,:,Top_chan)-datSp(1,:,Top_chan)))*ops.bitmVolt*ops.scaleproc;
                        
                        all_channel(irun + (1:numel(time)),:) = squeeze(datSp(nt0min,:,:)-datSp(1,:,:))*ops.bitmVolt*ops.scaleproc;
                        irun  = irun + numel(time);
                        
                        if i==1
                            ioffset = 0;
                        else
                            ioffset = ops.ntbuff;
                        end
                        timesp   = time - ioffset;
                        %Saved variable STT in rez file(timestamps in
                        %sampling rate, timestamps in ms, cluter id,
                        %amplitude of spike, Cost function, batch number
                        STT = cat(2,double(timesp+nt0min) +(NT-ops.ntbuff)*(i-1),  round(( double(timesp+nt0min) +(NT-ops.ntbuff)*(i-1))/(ops.fs*1e-3)), double(id),double(best(time)),double(Best_channel),double(Top_chan),double(i*ones(length(time),1)));
                        st3 = cat(1, st3, STT);
                end
            end

%****************************************Start***********************************************%
%***********************Substatction of detected spikes using templates**********************%
%********************************************************************************************%
            if(iteration>1) 
                if(~isempty(time))
                    A = dWU(:,:,id);
                    Mask = logical(U(:,id));
                    inds = repmat(time', nt0, 1) + repmat(int32(0:nt0-1)', 1, numel(time));
                    subs = dat(inds(:),:);
                    if(length(time)>1)
                        subs = reshape(subs,[size(inds) Nchan]);
                        subs = permute(subs,[1 3 2]);
                        Align = squeeze(subs(nt0min,:,:)).*Mask-squeeze(A(nt0min,:,:)).*Mask;
                        Align = reshape(Align,1,Nchan,length(time));
                        subs = subs - A + Align ;
                        substract = permute(subs,[1 3 2]);
                        substract = reshape(substract,[length(time)*nt0 Nchan]);
                        dat(inds(:),:) = substract;
                        DATA(:,:,i) = int16(dat*ops.scaleproc);
                    else
                        Align = squeeze(subs(nt0min,:)).*Mask'-squeeze(A(nt0min,:)).*Mask';
                        subs = subs - A + Align ;
                        dat(inds(:),:) = subs;
                        DATA(:,:,i) = int16(dat*ops.scaleproc);
                    end  
                end
            else%if it is 1st iteration then update sum the detected waveform of each clusters
                 inds = repmat(time', nt0, 1) + repmat(int32(1:nt0)', 1, numel(time));
                 subs = dat(inds(:),:);
                 subs = reshape(subs,[size(inds) Nchan]);
                 for ite=1:Nfilt
                    id_cluster = find(id==ite);  
                    suma0(:,:,ite) =  suma0(:,:,ite) + squeeze(sum(subs(:,id_cluster,:),2));
                    nb_spikes0(ite,1) = nb_spikes0(ite,1) + length(id_cluster);
                 end                       
            end
%******************************************END************************************************%
%***********************Substatction of detected spikes using templates**********************%
          
        end%end of loop on batch number
        
        
        %reevaluate the average template for each clusters or preserved them
        if(ops.freeze =='Y')  %Option use to overwrite the template to preserve the template provided by user              
           dWU = dWU0;
          else
           if(iteration==1)
                for il=1:Nfilt
                    if(nb_spikes0(il,:)==0)
                    else
                     dWU(:,:,il) = suma0(:,:,il)/sum(nb_spikes0(il,:));
                    end
                end
                [dWU,top] = SVD_topchan(dWU,Chan_criteria);
                dWU = alignWU(dWU,top, nt0min); %make sure that the deep is set at nt0min vital for substracting the template from signal
            end
        end
        
       %user output
       text =['Iteration:',num2str(iteration),'   Number of spikes:',num2str(nbsp),'   Threshold:',num2str(criteria_NCC)];
       disp(text)
       
       %Calculation will be stop if no more spikes are detected or number
       %of spikes are below a precision provided by user (ops.max_itera)
       if(iteration>1)
            nbsp_tot = nbsp_tot +nbsp;
             if(nbsp/nbsp_tot<ops.max_itera ||nbsp==0)
                    stop =1;
             end           
       end
       text =['Cumulated Spikes:',num2str(nbsp_tot)];
       disp(text)        
       iteration = iteration + 1;
end
all = all_channel;





