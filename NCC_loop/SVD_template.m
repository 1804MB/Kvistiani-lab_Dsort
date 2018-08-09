function  [dWU,W, U, Weight,Topchan] = SVD_template(dWU, Nrank,criteria)

[nt0, Nchan ,Nfilt] = size(dWU);

W = zeros(nt0,  Nfilt, Nrank ,'single');
U = zeros(Nchan, Nfilt,Nrank, 'single');
Weight = zeros(Nrank,Nfilt,'single');
Topchan = zeros(Nfilt,1);
dWU(isnan(dWU)) = 0;
for i=1:Nfilt
    if(isempty(nonzeros(dWU(:,:,i))))
    else
     [Ws,Ds,Us] = svd(dWU(:,:,i));
     %To get the correct sign for the peak deep
     [~, imax] = max(abs(Ws(:,1)));
     Us(:,1) = -Us(:,1) * sign(Ws(imax,1));
     Ws(:,1) = -Ws(:,1) * sign(Ws(imax,1));
     eigen = diag(Ds);
     Weight(:,i) = (eigen(1:Nrank).^2)/sum(sum(Ds(:,:).^2)) ; 
     W(:,i,:) = Ws(:,1:Nrank)*Ds(1:Nrank,1:Nrank);
     
     U_s = Us;   
     U_s(:,1) = U_s(:,1)/max(Us(:,1));
     id = find(U_s(:,1)>=criteria);
     id_del = find(U_s(:,1)<criteria);
     Us(id_del,1:Nrank) = 0;
     
         
     U(:,i,:) = Us(:,1:Nrank);
     [~,chanl] = max(U(:,i,1));
     Topchan(i) = chanl; 
     dWU(:,:,i) = squeeze(W(:,i,:))*squeeze(U(:,i,:))';
    end
end
