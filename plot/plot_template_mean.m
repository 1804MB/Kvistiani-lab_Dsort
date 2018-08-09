function plot_template_mean(rez,i)


     temp = rez.Merge_cluster{i,1};
     nb_temp = length(rez.Merge_cluster{i,1});    
     Channel = rez.Merge_cluster{i,5};
     f=figure;
     for k = 1:length(Channel)
         text = [];
          subplot(ceil(length(Channel)/2),2,k)
          plot(squeeze(rez.M_template(:,Channel(k),i)), 'LineWidth',3);
          hold on;
          text{1} = 'mean';
          for jk=1:nb_temp
           plot(squeeze(rez.template(:,Channel(k), temp(jk))),'-.', 'LineWidth',1.5);
           text{1+jk} = ['template:',num2str(temp(jk))];
           hold on;
           
          end
           title([ 'Channel: ', num2str(Channel(k))])
           xlim([0 rez.ops.nt0]);
           xlabel('Sampling point')
           ylabel('Signal in mV')
          legend(text)
          movegui(f,'northeast');
          hold off
     end
