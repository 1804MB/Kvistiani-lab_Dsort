function [rez]=probability(rez)
Nrank = 3;
M_clust = length(unique(rez.st(:,end)));
h = waitbar(0, 'Computing joint probabilities...');
rez.PDF =cell(M_clust,1,1);
for i = 1:M_clust
    Channel = rez.Merge_cluster{i,5};
    Nchan = length(Channel);
    PC = rez.PC{i};
    if(length(PC)>500)
        for k = 1:Nchan
            Progress = (i/M_clust);
            waitbar(Progress)
            %extract the principal component of the cluster i for the largest
            %extract the principal component of the cluster i for the largest
            %projection
            score =  zeros(length(PC),Nrank);
            for irank =1:Nrank
                score(:,irank) = single(PC(k,irank,:));
            end
            
            % t-distribution mixture computation
            try                                       % to deal with batch processing , remove later. starts here
                [mu_t, S, nu] = fitt_fast(score);
            catch
                rm = rez.st(:,end)==i;
                rez.st(rm,:) = [];
                rez.best_channel_M(rm,:)=[];
                continue
            end                                      % ends here
            sigma_t = sqrt(diag(S));
            y_t = zeros(length(score),Nrank);
            % Compute the probability for a given spike to belong to cluster i
            for irank = 1:Nrank
                pd_t = makedist('tlocationscale','mu',mu_t(irank),'sigma',sigma_t(irank), 'nu', nu);
                y_t(:,irank) = pdf(pd_t,score(:,irank));
            end
            
            pdf_c = prod(y_t,2);
            pdf_c = pdf_c./max(pdf_c);
            rez.PDF{i,k,:} = pdf_c;
        end
    end
end

close(h)
