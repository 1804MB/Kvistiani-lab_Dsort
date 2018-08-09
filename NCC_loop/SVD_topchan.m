function  [dWU,top] = SVD_topchan(dWU,criteria)

[~, ~ ,Nfilt] = size(dWU);

top = zeros(Nfilt,1);
dWU(isnan(dWU)) = 0;
for i=1:Nfilt
    if(isempty(nonzeros(dWU(:,:,i))))
    else
     [Ws,~,Us] = svd(dWU(:,:,i));
     %To get the correct sign for the peak deep
     [~, imax] = max(abs(Ws(:,1)));
     Us(:,1) = -Us(:,1) * sign(Ws(imax,1));
     [maxU,maxi] = max(Us(:,1));
     top(i) = maxi;
     Us(:,1) = Us(:,1)/maxU;
     id_del = find(Us(:,1)<criteria);
     dWU(:,id_del,i)= 0;
    end
end
