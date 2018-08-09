%*************************************************************************%
%*************************************************************************%
%fitTemplates: Detect spike and calculate the templates
%Simplified and commented by Madeny belkhiri 30/01/18.
%*************************************************************************%
%*************************************************************************%

function rez = NCC_overlap(rez, DATA, uproj)

rng('default');
rng(1);

nt0             = rez.ops.nt0;       % Number of points in sampling rate for a waveform
nt0min         = rez.ops.nt0min;    % Strange constant added at the timesstamps why [MB?]
ops = rez.ops;                       % Load all the parameters from the config file and preprocess
Nbatch  = rez.temp.Nbatch;           % Number of Batch
Nfilt 	= ops.Nfilt;                 % Number of Clusters provided by user via config file
NT  	= ops.NT;                    % Number of element in a batch
Nrank   = ops.Nrank;                 % Rank to use for the SVD decomposition given by user via config file
Nchan 	= ops.Nchan;                 % Number Channels provided by user via config file
last = 0;
iperm = randperm(Nbatch);            % Randomize the batch id, basically shuffle the batch order 
st3 = [];
rez.st3 = [];

	

% Deleted the option to choose from data or not, now only fromdata.
% thus, the initiale template will always be calculated from the data.
if(ops.template =='Y') % look at the tempaltes i give
    dWU0   = load(ops.template_file);
    dWU    = dWU0.T;
    
elseif(ops.template=='G') % or it generates templates based on config file
    [WUinit] = init_template(ops);
    dWU    = WUinit(:,:,1:Nfilt);
else
    WUinit = optimizePeaks(ops,uproj);%does a scaled kmeans to determine the  Nfilt clusters template
    dWU    = WUinit(:,:,1:Nfilt);
end

%just the shuffled batch
miniorder = repmat(iperm, 1, 1);
fprintf('Time %3.0fs. Optimizing templates ...\n', toc)
load PCspikes
%Start the loop to detect spikes and optimize the template,
criteria_NCC = ops.criteria_NCC;
irun =0;


maxiter =10;  %20 nr of itreration you want to run


suma = zeros(nt0,Nchan,Nfilt);%To calculate the average template
nb_spikes =  zeros(Nfilt,maxiter) ;
suma0 = zeros(nt0,Nchan,Nfilt);%To calculate the average template
nb_spikes0 =  zeros(Nfilt,1) ;
nbsp_tot = 0;

for iteration =1:maxiter
   if iteration>1
       last =1;
   end

   nbsp = 0;
   % parameter update
   %Decompose the template dWU (called K in the paper) in two matrix W*U using
   %SVD.
   % Using svd one get dWU = W*D*U, in fact output W is W*D, D being the
   % diagonal matrix containing the eigenvalues
  
   [W, U, Weight] = SVD_template(dWU, Nrank);
    Weight = mean(Weight,2);
    W = reshape(W,nt0,Nrank*Nfilt);
   %Some plot of the template, number of spikes and mask
   figure;
   subplot(1,2,1)
   plot(W(:,:,1))
   xlim([0 nt0-1]);
   title('W(:,:,1)')            
   subplot(1,2,2)
   imagesc(U(:,:,1))
   title('U(:,:,1)')     
   drawnow;
   Params = int32([NT Nfilt*Nrank 1 nt0]); %Used for the convolution
        for i=1:Nbatch 
               time = [];
               id = [];
               score_s = [];
            % some of the parameters change with iteration number
            % select batch and load from RAM or disk
            ibatch = miniorder(i);
            %Take a random batch of the whitened and filtered data
             dat = single(DATA(:,:,ibatch))/rez.ops.scaleproc;  
             data  = zeros(NT, Nfilt, Nrank, 'single');
             for irank = 1:Nrank
                 data(:,:,irank) = dat * U(:,:,irank);
             end
             
             data = reshape(data, NT, Nfilt*Nrank);
             % Call CUDA to execute convolution
             [tconv] = mexconv2(Params,data,gpuArray(W));           
             tconv = gather(tconv);
             tconv (isnan(tconv ))=0;
             
             %Calculate the normalized cross-correlation using a weights of
             %eigenvalue of each components
             tclus= zeros(NT,Nfilt,'single');

             %Calculate the normalized cross-correlation using a weights of
             %eigenvalue of each components
             for ki=1:Nfilt
                for irank =1:Nrank
                     tclus(:,ki) =  tclus(:,ki) + tconv(:,ki+(irank-1)*Nfilt).*Weight(irank);
                end
             end
             %0 all the point below the score threshold
             clear tconv;
             %find the higest Normalized-Cross-Correlation(NCC) value for each
             %sampling point
             [best,it]= max(tclus,[],2);
             %Apply the criteria
             best(best<criteria_NCC)=0;
             %Search for the spikes timestamps if score goes above threshold
             %Find the highest NCC value every nt0 points (should correspond to 1 ms)
             [time] =  get_timestamps(best,nt0);
             itera = iteration*ones(length(time),1);
             %Get the corresponding Cluster id
             id = it(time);
             %number of spikes
             nbsp = nbsp + length(time);
             time  = int32(time);
             id = int32(id);

            
            %Save the detected spikes and projections
            if(last==1)
                if ~isempty(time)
                       
                        if ibatch==1
                            ioffset = 0;
                        else
                            ioffset = ops.ntbuff;
                        end
                        timesp   = time - ioffset;
                        %Saved variable STT in rez file(timestamps in
                        %sampling rate, timestamps in ms, cluter id,
                        %amplitude of spike, Cost function, batch number
                        STT = cat(2,double(timesp+nt0min) +(NT-ops.ntbuff)*(ibatch-1),  round(( double(timesp+nt0min) +(NT-ops.ntbuff)*(ibatch-1))/(rez.ops.fs*1e-3)), double(id),double(best(time)),double(itera),double(time),double(ibatch*ones(length(time),1)));
                        st3 = cat(1, st3, STT);
                end
                irun = irun + numel(time);
            end
             %Sum all the waveform for a given cluster and count number of
            %spikes over iterations
            if(iteration>1)
                 for ite =1:Nfilt
                    id_cluster = find(id==ite);  
                    time_sub = time(id_cluster);
                    id_U =  find(U(:,ite,1)>0);
                    if(~isempty(time_sub))
                        inds = repmat(time_sub', nt0, 1) + repmat(int32(0:nt0-1)', 1, numel(time_sub));
                        subs = dat(inds(:),:);
                        subs = reshape(subs,[size(inds) Nchan]);
                        suma(:,:,ite) =  suma(:,:,ite) + squeeze(sum(subs,2));
                        nb_spikes(ite,iteration) = nb_spikes(ite,iteration) + length(time_sub);
                        subs = permute(subs,[1 3 2]);

                        subs(:,id_U,:) = subs(:,id_U,:)-squeeze(dWU(:,id_U,ite))+(subs(nt0min,id_U,:)-dWU(nt0min,id_U,ite));
                        substract = permute(subs,[1 3 2]);
                        substract = reshape(substract,[length(id_cluster)*nt0 Nchan]);
                        dat(inds(:),:) = substract;
                       DATA(:,:,ibatch) = int16(dat*rez.ops.scaleproc);
                    else
                    end
                 end
                
            else
                inds = repmat(time', nt0, 1) + repmat(int32(1:nt0)', 1, numel(time));
                subs = dat(inds(:),:);
                subs = reshape(subs,[size(inds) Nchan]);
                for ite=1:Nfilt
                    id_cluster = find(id==ite);  
                    suma0(:,:,ite) =  suma0(:,:,ite) + squeeze(sum(subs(:,id_cluster,:),2));
                    nb_spikes0(ite,1) = nb_spikes0(ite,1) + length(id_cluster);
                end                       
            end
            
        end
                                               
%          dWU = zeros(nt0,Nchan,Nfilt);
%         %reevaluate the average template for each clusters
        if(ops.freeze =='Y')  %Option use to overwrite the template to preserve the template provided by user              
           dWU = dWU0;
           dWU = gpuarray(dWU);
        else
%            if(iteration>1)
% 
%                 for il=1:Nfilt
%                     if(sum(nb_spikes(il,:))==0)
%                     else
%                      dWU(:,:,il) = suma(:,:,il)/sum(nb_spikes(il,:));
%                     end
%                 end
%             else
%                 for il=1:Nfilt
%                     if(nb_spikes0(il,:)==0)
%                     else
%                      dWU(:,:,il) = suma0(:,:,il)/nb_spikes0(il,1);
%                     end
%                 end
%             end
        end

       text =['Iteration:',num2str(iteration),'   Number of spikes:',num2str(nbsp),'   Threshold:',num2str(criteria_NCC)];
       disp(text)
       if(iteration>1)
       nbsp_tot = nbsp_tot +nbsp;
       end
       text =['Cumulated Spikes:',num2str(nbsp_tot)];
       disp(text)

end
% Save template
[W, U, Weight] = SVD_template(dWU, Nrank);
W = reshape(W,nt0,Nrank*Nfilt);
rez.dWU = dWU;
rez.W = W;
rez.U = U;
figure;
subplot(1,2,1)
plot(W(:,:,1))
xlim([0 nt0-1]);
title('after W(:,:,1)')            
subplot(1,2,2)
imagesc(U(:,:,1))
title('after U(:,:,1)')
rez.Weight = Weight;
rez.nspikes     = nbsp;
[~, isort]      = sort(st3(:,1), 'ascend');
%Save timestamps
st3             = st3(isort,:);
rez.st3         = st3;
% sort best channels
[~, iNch]       = sort(abs(U(:,:,1)), 1, 'descend');
rez.iNeighPC    = iNch;

%Rearrange PC


% 
