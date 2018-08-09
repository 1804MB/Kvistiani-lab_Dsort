function dWU = alignWU(dWU, top,nt0min)

[nt0 ,~, Nfilt] = size(dWU);
dWU(isnan(dWU)) = 0;
for i = 1:Nfilt
    if(isempty(nonzeros(dWU(:,:,i))))
    else
    [~, imax] = min(squeeze(dWU(:,top(i),i)));
    dmax = -(imax - nt0min);
    if(dmax>0)
        dWU((dmax + 1):nt0,top(i), i) = dWU(1:nt0-dmax,top(i), i);
    elseif(dmax<0)
        dWU(1:nt0+dmax,top(i), i) = dWU((1-dmax):nt0,top(i), i);
    else
    end
    end
end




