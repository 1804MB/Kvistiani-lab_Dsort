function [snr_coeff] = estimate_snr(rez,N)
  dWU = rez.M_template/rez.ops.bitmVolt;
  [~,~,Nclust] = size(dWU); 
  sigma =  rez.sigdata;
  mu  =  rez.mudata;
  Nbatch = length(sigma);
  nt0 = rez.ops.nt0;
   snr_coeff = zeros(Nclust,Nbatch);
  for i=1: Nclust
        channel =  rez.Merge_cluster{i,5};
        Template = reshape(dWU(:,channel,i),[length(channel)*nt0 1]);
    for ibatch = 1:Nbatch
      Noise =  sigma(channel,ibatch).*randn(length(channel),nt0,N) + mu(channel,ibatch);
      Noise =  reshape(Noise,[length(channel)*nt0 N]);
      snr_s = zeros(1,N);
      for j=1:N
        snr_s(j) = snr(Template,Noise(:,j));
      end
      snr_coeff(i,ibatch) = mean(snr_s);
    end
  end


end