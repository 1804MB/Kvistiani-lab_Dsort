function plot_pc(rez,i)
    f=figure;
     channel  = rez.Merge_cluster{i,5};
     Nchan = length(channel);
     PC = rez.PC{i};
     ki =0;
     for k=1:Nchan
        subplot(Nchan,3,1+ki)  
        plot(squeeze(PC(k,1,:)),squeeze(PC(k,2,:)),'*')
        title(['Cluster:', num2str(i)]);
        ylabel('PC2')
        xlabel('PC1')
        xlim([min(squeeze(PC(k,1,:))) max(squeeze(PC(k,1,:)))])     
        ylim([min(squeeze(PC(k,2,:))) max(squeeze(PC(k,2,:)))])     
        hold on
        subplot(Nchan,3,2+ki)   
        plot(squeeze(PC(k,1,:)),squeeze(PC(k,3,:)),'*')
        ylabel('PC3')
        xlabel('PC1')
        xlim([min(squeeze(PC(k,1,:))) max(squeeze(PC(k,1,:)))])     
        ylim([min(squeeze(PC(k,3,:))) max(squeeze(PC(k,3,:)))])     
        title(['Channel:',num2str(channel(k))])
        hold on
        subplot(Nchan,3,3+ki)   
        plot(squeeze(PC(k,2,:)),squeeze(PC(k,3,:)),'*')
        ylabel('PC3')
        xlabel('PC2')
        xlim([min(squeeze(PC(k,2,:))) max(squeeze(PC(k,2,:)))])     
        ylim([min(squeeze(PC(k,3,:))) max(squeeze(PC(k,3,:)))])     
        hold on
        ki = ki +3;
     end
     movegui(f,'northwest');
         hold off
         
