function plot_proba(rez,i)
     f=figure;
     channel  = rez.Merge_cluster{i,5};
     Nchan = length(channel);
     
     h=0.001;
     for k=1:Nchan
        PDF = rez.PDF{i,k};
        subplot(Nchan,1,k)        
        N = histcounts(PDF,0:h:1);
        bar(0:h:1-h,(N))
        title(['Channel:', num2str(channel(k))]);
        ylabel('Spike count')
        xlabel('Probability')
        xlim([-0.01 1.01])
        hold on
     end
         hold off
         movegui(f,'southwest');
