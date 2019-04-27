function [T,rez] = estimate_snr(rez,T)
  [~,~, U, ~] = SVD_template(T,3,rez.ops.Chan_criteria);
  [~,~,Nclust] = size(T); 
  [Chan] = Chan_cluster2(U,Nclust,rez.ops.Chan_criteria);
  sigma =  rez.sigdata;
  mu  =  rez.mudata;
  Nbatch = length(sigma);
  nt0 = rez.ops.nt0;
  th = rez.ops.snr_th;
  snr_coeff_batch = zeros(Nclust,Nbatch);
  snr_mean = zeros(Nclust,1);
  N  =1000;
  list = [];
  for i=1:Nclust
        channel =  Chan{i};
        Template = reshape(T(:,channel,i),[length(channel)*nt0 1]);
    for ibatch = 1:Nbatch
      Noise =  sigma(channel,ibatch).*randn(length(channel),nt0,N) + mu(channel,ibatch);
      Noise =  reshape(Noise,[length(channel)*nt0 N]);
      snr_s = zeros(1,N);
      for j=1:N
        snr_s(j) = snr(Template,Noise(:,j));
      end
      snr_coeff_batch(i,ibatch) = mean(snr_s);
    end
    snr_mean(i) = mean(snr_coeff_batch(i,:));
    
    
if length(channel) > 4
    channel = channel(1:4);
end

    if(snr_mean<th(length(channel)))
       list  = [list,i];
    end
  end
   T(:,:,list) = [];
      rez.snr_coeff_batch  = snr_coeff_batch;
      rez.snr_mean  = snr_mean;
end