function plot_pc_init_two(PC,Nchan,i,id1,id2)

                 figure;
             ki =0;
                   
            for k=1:Nchan
         
            subplot(Nchan,3,1+ki)  
            plot(squeeze(PC(k,1,id1)),squeeze(PC(k,2,id1)),'b.','MarkerSize',2)
            hold on;
            plot(squeeze(PC(k,1,id2)),squeeze(PC(k,2,id2)),'r.','MarkerSize',2)
            title(['Cluster:', num2str(i)]);
            ylabel('PC2')
            xlabel('PC1')
            xlim([min(squeeze(PC(k,1,:))) max(squeeze(PC(k,1,:)))])     
            ylim([min(squeeze(PC(k,2,:))) max(squeeze(PC(k,2,:)))])     
            hold on
            subplot(Nchan,3,2+ki)   
            plot(squeeze(PC(k,1,id1)),squeeze(PC(k,3,id1)),'b.','MarkerSize',2)
            hold on;
            plot(squeeze(PC(k,1,id2)),squeeze(PC(k,3,id2)),'r.','MarkerSize',2)
            ylabel('PC3')
            xlabel('PC1')
            xlim([min(squeeze(PC(k,1,:))) max(squeeze(PC(k,1,:)))])     
            ylim([min(squeeze(PC(k,3,:))) max(squeeze(PC(k,3,:)))])     
            title(['Channel:',num2str(k)])
            hold on
            subplot(Nchan,3,3+ki)   
            plot(squeeze(PC(k,2,id1)),squeeze(PC(k,3,id1)),'b.','MarkerSize',2)
            hold on;
            plot(squeeze(PC(k,2,id2)),squeeze(PC(k,3,id2)),'r.','MarkerSize',2)
            ylabel('PC3')
            xlabel('PC2')
            xlim([min(squeeze(PC(k,2,:))) max(squeeze(PC(k,2,:)))])     
            ylim([min(squeeze(PC(k,3,:))) max(squeeze(PC(k,3,:)))])     
            hold on
            ki = ki +3;
            end
            