function  [Nrank] = get_rank(dWU, variance)

[~, ~ ,Nfilt] = size(dWU);
Rank = [];
j = 1;
for i=1:Nfilt
    if(isempty(nonzeros(dWU(:,:,i))))
    else
     j = j+1;
     [~,Ds,~] = svd(dWU(:,:,i));
     %To get the correct sign for the peak deep
     normsqS = sum(Ds.^2);                        %// total variance
     Rank(j) = find(cumsum(Ds.^2)/normsqS >= variance, 1) ; %// number of component to keep to reach a precision of 90% of the initial variance
    end
end
Nrank  = ceil(mean(Rank));
if(Nrank==1)
    Nrank =2;
end

