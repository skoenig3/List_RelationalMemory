
distancethreshold = [0  24 48 72 96  120  144 168 200 400;...%in pixels 24 pixels/dva
    24 48 72 96 120 144  168 200 400 1000];%in pixels 24 pixels/dva

nov_totalfixations = zeros(1,size(distancethreshold,2));
nov_allSalIOR = cell(1,size(distancethreshold,2));
rep_totalfixations = zeros(1,size(distancethreshold,2));
rep_allSalIOR = cell(1,size(distancethreshold,2));
for file = 1:length(cortex_files);
   load([data_dir cortex_files{file}(1:8) '-IOR'])
    for d = 1:length(distancethreshold);
        nov_allSalIOR{d} = [nov_allSalIOR{d}; nov_returnfixsal{d}];
        nov_totalfixations(d) = nov_totalfixations(d)+size(nov_returnfixsal{d},1);
        
        rep_allSalIOR{d} = [rep_allSalIOR{d}; rep_returnfixsal{d}];
        rep_totalfixations(d) = rep_totalfixations(d)+size(rep_returnfixsal{d},1);
    end
end

%%


nov_return_sal_pvalues = NaN(1,size(distancethreshold,2));
nov_return_sal_means = NaN(size(distancethreshold,2),2);
nov_return_sal_stds = NaN(size(distancethreshold,2),2);
nov_numreturns = NaN(1,size(distancethreshold,2));

rep_return_sal_pvalues = NaN(size(distancethreshold,2),3);
rep_return_sal_means = NaN(size(distancethreshold,2),2);
rep_return_sal_stds = NaN(size(distancethreshold,2),2);
rep_numreturns = NaN(1,size(distancethreshold,2));

labels = {};
for d = 1:size(distancethreshold,2);
    [~,p]=ttest2(nov_allSalIOR{d}(:,4),nov_allSalIOR{d}(:,8));
    nov_return_sal_pvalues(d) = p; %prior vs return
    nov_return_sal_means(d,1) = mean(nov_allSalIOR{d}(:,4));%average salience at prior
    nov_return_sal_means(d,2) = mean(nov_allSalIOR{d}(:,8));%average salience at return
    nov_return_sal_stds(d,1) = std(nov_allSalIOR{d}(:,4));%std of salience at prior
    nov_return_sal_stds(d,2) = std(nov_allSalIOR{d}(:,8));%std of salience at return
    nov_numreturns(d) = length(nov_allSalIOR{d}(:,8));
    
    [~,p]=ttest2(rep_allSalIOR{d}(:,4),rep_allSalIOR{d}(:,8));
    rep_return_sal_pvalues(d) = p; %prior vs return
    rep_return_sal_means(d,1) = mean(rep_allSalIOR{d}(:,4));%average salience at prior
    rep_return_sal_means(d,2) = mean(rep_allSalIOR{d}(:,8));%average salience at return
    rep_return_sal_stds(d,1) = std(rep_allSalIOR{d}(:,4));%std of salience at prior
    rep_return_sal_stds(d,2) = std(rep_allSalIOR{d}(:,8));%std of salience at return
    rep_numreturns(d) = length(rep_allSalIOR{d}(:,8));
    
    labels{d} = [num2str(round(distancethreshold(1,d)/24)) '-' num2str(round(distancethreshold(2,d)/24))];
end
%%
figure
subplot(1,2,1)
hold on
errorbar(nov_return_sal_means(:,2),nov_return_sal_stds(:,2)./sqrt(nov_numreturns'))
errorbar(nov_return_sal_means(:,1),nov_return_sal_stds(:,1)./sqrt(nov_numreturns'))
ylim([0.32 0.42])
for d = 1:length(nov_return_sal_pvalues)
    if nov_return_sal_pvalues(d) < 0.05/length(nov_return_sal_pvalues) %if the return is
        plot(d,0.4,'k*')
    end
end
xlim([0 length(distancethreshold)+1])
xlabel('Distance Tresholds (dva)')
ylabel('Salience (a.u.)')
set(gca,'XTick',1:length(distancethreshold))
set(gca,'XTickLabel',labels,'XTickLabelRotation',45)
hold off
box off
title('Novel Images')
legend('Return','Prior')


subplot(1,2,2)
hold on
errorbar(rep_return_sal_means(:,2),rep_return_sal_stds(:,2)./sqrt(rep_numreturns'))
errorbar(rep_return_sal_means(:,1),rep_return_sal_stds(:,1)./sqrt(rep_numreturns'))
ylim([0.25 0.45])
for d = 1:length(rep_return_sal_pvalues)
    if rep_return_sal_pvalues(d) < 0.05/length(rep_return_sal_pvalues) %if the return is
        plot(d,0.4,'k*')
    end
end
xlim([0 length(distancethreshold)+1])
xlabel('Distance Tresholds (dva)')
ylabel('Salience (a.u.)')
set(gca,'XTick',1:length(distancethreshold))
set(gca,'XTickLabel',labels,'XTickLabelRotation',45)
hold off
box off
title('Repeated Images')
subtitle(cortex_files{1}(1:2))

%%

nov_return_dur_pvalues = NaN(1,size(distancethreshold,2));
nov_return_dur_means = NaN(size(distancethreshold,2),2);
nov_return_dur_stds = NaN(size(distancethreshold,2),2);
nov_numreturns = NaN(1,size(distancethreshold,2));

rep_return_dur_pvalues = NaN(size(distancethreshold,2),3);
rep_return_dur_means = NaN(size(distancethreshold,2),2);
rep_return_dur_stds = NaN(size(distancethreshold,2),2);
rep_numreturns = NaN(1,size(distancethreshold,2));

labels = {};
for d = 1:size(distancethreshold,2);
    [~,p]=ttest2(nov_allSalIOR{d}(:,10),nov_allSalIOR{d}(:,11));
    nov_return_dur_pvalues(d) = p; %prior vs return
    nov_return_dur_means(d,1) = mean(nov_allSalIOR{d}(:,10));%average salience at prior
    nov_return_dur_means(d,2) = mean(nov_allSalIOR{d}(:,11));%average salience at return
    nov_return_dur_stds(d,1) = std(nov_allSalIOR{d}(:,10));%std of salience at prior
    nov_return_dur_stds(d,2) = std(nov_allSalIOR{d}(:,11));%std of salience at return
    nov_numreturns(d) = length(nov_allSalIOR{d}(:,11));
    
    [~,p]=ttest2(rep_allSalIOR{d}(:,10),rep_allSalIOR{d}(:,11));
    rep_return_dur_pvalues(d) = p; %prior vs return
    rep_return_dur_means(d,1) = mean(rep_allSalIOR{d}(:,10));%average salience at prior
    rep_return_dur_means(d,2) = mean(rep_allSalIOR{d}(:,11));%average salience at return
    rep_return_dur_stds(d,1) = std(rep_allSalIOR{d}(:,10));%std of salience at prior
    rep_return_dur_stds(d,2) = std(rep_allSalIOR{d}(:,11));%std of salience at return
    rep_numreturns(d) = length(rep_allSalIOR{d}(:,11));
    
    labels{d} = [num2str(round(distancethreshold(1,d)/24)) '-' num2str(round(distancethreshold(2,d)/24))];
end

figure
subplot(1,2,1)
hold on
errorbar(nov_return_dur_means(:,2),nov_return_dur_stds(:,2)./sqrt(nov_numreturns'))
errorbar(nov_return_dur_means(:,1),nov_return_dur_stds(:,1)./sqrt(nov_numreturns'))
for d = 1:length(nov_return_dur_pvalues)
    if nov_return_dur_pvalues(d) < 0.05/length(nov_return_dur_pvalues) %if the return is
        plot(d,225,'k*')
    end
end
xlim([0 length(distancethreshold)+1])
xlabel('Distance Tresholds (dva)')
ylabel('Salience (a.u.)')
set(gca,'XTick',1:length(distancethreshold))
set(gca,'XTickLabel',labels,'XTickLabelRotation',45)
hold off
box off
title('Novel Images')
legend('Return','Prior')
yl1 = ylim;

subplot(1,2,2)
hold on
errorbar(rep_return_dur_means(:,2),rep_return_dur_stds(:,2)./sqrt(rep_numreturns'))
errorbar(rep_return_dur_means(:,1),rep_return_dur_stds(:,1)./sqrt(rep_numreturns'))
for d = 1:length(rep_return_dur_pvalues)
    if rep_return_dur_pvalues(d) < 0.05/length(rep_return_dur_pvalues) %if the return is
        plot(d,225,'k*')
    end
end
xlim([0 length(distancethreshold)+1])
xlabel('Distance Tresholds (dva)')
ylabel('Salience (a.u.)')
set(gca,'XTick',1:length(distancethreshold))
set(gca,'XTickLabel',labels,'XTickLabelRotation',45)
hold off
box off
title('Repeated Images')
subtitle(cortex_files{1}(1:2))
yl2 = ylim;

ymin = min([yl1(1) yl2(1)]);
ymax = max([yl1(2) yl2(2)]);
subplot(1,2,1),ylim([ymin ymax])
subplot(1,2,2),ylim([ymin ymax])
%%
nov_mean_num_fixations_between = NaN(1,length(distancethreshold));
nov_std_num_fixations_between = NaN(1,length(distancethreshold));
nov_num_fixations = NaN(1,length(distancethreshold));

rep_mean_num_fixations_between = NaN(1,length(distancethreshold));
rep_std_num_fixations_between = NaN(1,length(distancethreshold));
rep_num_fixations = NaN(1,length(distancethreshold));

nov_mean_time_between = NaN(1,length(distancethreshold));
nov_std_time_between = NaN(1,length(distancethreshold));

rep_mean_time_between = NaN(1,length(distancethreshold));
rep_std_time_between = NaN(1,length(distancethreshold));
for d = 1:length(distancethreshold);
    nov_mean_num_fixations_between(d) =mean(nov_allSalIOR{d}(:,13)-nov_allSalIOR{d}(:,12));
    nov_std_num_fixations_between(d) = std(nov_allSalIOR{d}(:,13)-nov_allSalIOR{d}(:,12));
    nov_num_fixations = size(nov_allSalIOR{d},1);
    
    rep_mean_num_fixations_between(d) =mean(rep_allSalIOR{d}(:,13)-rep_allSalIOR{d}(:,12));
    rep_std_num_fixations_between(d) = std(rep_allSalIOR{d}(:,13)-rep_allSalIOR{d}(:,12));
    rep_num_fixations = size(rep_allSalIOR{d},1);
    
    nov_mean_time_between(d) =mean(nov_allSalIOR{d}(:,7)-nov_allSalIOR{d}(:,3));
    nov_std_time_between(d) = std(nov_allSalIOR{d}(:,7)-nov_allSalIOR{d}(:,3));
    
    rep_mean_time_between(d) =mean(rep_allSalIOR{d}(:,7)-rep_allSalIOR{d}(:,3));
    rep_std_time_between(d) = std(rep_allSalIOR{d}(:,7)-rep_allSalIOR{d}(:,3));
end
 

% figure
subplot(1,2,1)
hold on
errorbar(nov_mean_num_fixations_between,nov_std_num_fixations_between./sqrt(nov_num_fixations))
errorbar(rep_mean_num_fixations_between,rep_std_num_fixations_between./sqrt(rep_num_fixations))
hold off
xlim([0 length(distancethreshold)+1])
xlabel('Distance Tresholds (dva)')
ylabel('# of Fixations')
set(gca,'XTick',1:length(distancethreshold))
set(gca,'XTickLabel',labels,'XTickLabelRotation',45)
title('Number of Fixations Between Prior and Return Fixation')

subplot(1,2,2)
hold on
errorbar(nov_mean_time_between,nov_std_time_between./sqrt(nov_num_fixations))
errorbar(rep_mean_time_between,rep_std_time_between./sqrt(rep_num_fixations))
hold off
xlim([0 length(distancethreshold)+1])
xlabel('Distance Tresholds (dva)')
ylabel('Time (ms)')
set(gca,'XTick',1:length(distancethreshold))
set(gca,'XTickLabel',labels,'XTickLabelRotation',45)

title('Time Between Prior and Return Fixation')
%%

%---Vivian---%
pre_files = {'PW150416.4','PW150417.2','PW150420.2','PW150421.2',...
    'PW150422.2','PW150423.2','PW150424.2','PW150427.2',...
    'PW150428.2','PW150429.2','PW150430.2','PW150501.2',...
    'PW150504.2','PW150505.2','PW150506.2','PW150507.2',...
    'PW150511.2'};

post_files = {'PW160310.2','PW160311.2','PW160314.2','PW160315.2','PW160316.2',...
    'PW160317.2','PW160318.2','PW160321.2','PW160322.2','PW160323.2',...
    'PW160324.2','PW160325.2','PW160325.2','PW160329.2','PW160330.2'};
% 
% %---Red---%
pre_files = {'RR150423.2','RR150424.2','RR150427.2','RR150428.2','RR150429.2',...
    'RR150430.2','RR150501.2','RR150504.2','RR150505.2','RR150506.2',...
    'RR150507.2','RR150508.2','RR150511.2','RR150512.2','RR150513.2',...
    'RR150515.2','RR150518.2','RR150519.2','RR150520.2','RR150521.2'};

post_files = {'RR160324.2','RR160325.1','RR160328.2','RR160330.2',...
                'RR160331.2','RR160401.2','RR160405.2','RR160406.2',...
                'RR160407.2','RR160408.2','RR160411.2','RR160412.2',...
                'RR160413.2','RR160414.2','RR160415.2','RR160418.2',...
                'RR160419.2','RR160420.2','RR160421.2','RR160422.2'};


pre_nov_set_means_num_fix_time_between = NaN(2,length(pre_files));
post_nov_set_means_num_fix_time_between = NaN(2,length(post_files));
pre_rep_set_means_num_fix_time_between = NaN(2,length(pre_files));
post_rep_set_means_num_fix_time_between = NaN(2,length(post_files));

for file = 1:length(pre_files);
   load([data_dir cortex_files{file}(1:8) '-IOR.mat'])
   
   nov_data = [nov_returnfixsal{1}; nov_returnfixsal{2}]; %combine data for 0-2 dva part
   rep_data = [rep_returnfixsal{1}; rep_returnfixsal{2}]; %combine data for 0-2 dva part

   pre_nov_set_means_num_fix_time_between(1,file) = mean(nov_data(:,13)-nov_data(:,12));
   pre_nov_set_means_num_fix_time_between(2,file) = mean(nov_data(:,7)-nov_data(:,3));
   
   pre_rep_set_means_num_fix_time_between(1,file) = mean(rep_data(:,13)-rep_data(:,12));
   pre_rep_set_means_num_fix_time_between(2,file) = mean(rep_data(:,7)-rep_data(:,3));
end

for file = 1:length(post_files);
   load([data_dir cortex_files{file}(1:8) '-IOR'])
   
   nov_data = [nov_returnfixsal{1}; nov_returnfixsal{2}]; %combine data for 0-2 dva part
   rep_data = [rep_returnfixsal{1}; rep_returnfixsal{2}]; %combine data for 0-2 dva part

   post_nov_set_means_num_fix_time_between(1,file) = mean(nov_data(:,13)-nov_data(:,12));
   post_nov_set_means_num_fix_time_between(2,file) = mean(nov_data(:,7)-nov_data(:,3));
   
   post_rep_set_means_num_fix_time_between(1,file) = mean(rep_data(:,13)-rep_data(:,12));
   post_rep_set_means_num_fix_time_between(2,file) = mean(rep_data(:,7)-rep_data(:,3));
end
%%

figure
subplot(1,2,1)
bar([mean(pre_nov_set_means_num_fix_time_between(1,:))  mean(post_nov_set_means_num_fix_time_between(1,:)); ...
    mean(pre_rep_set_means_num_fix_time_between(1,:)) mean(post_rep_set_means_num_fix_time_between(1,:))])
hold on
errorb([mean(pre_nov_set_means_num_fix_time_between(1,:))  mean(post_nov_set_means_num_fix_time_between(1,:)); ...
    mean(pre_rep_set_means_num_fix_time_between(1,:)) mean(post_rep_set_means_num_fix_time_between(1,:))],...
    [std(pre_nov_set_means_num_fix_time_between(1,:))./sqrt(length(pre_files)), ...
    std(post_nov_set_means_num_fix_time_between(1,:))./sqrt(length(post_files));...
    std(pre_rep_set_means_num_fix_time_between(1,:))./sqrt(length(pre_files)), ...
    std(post_rep_set_means_num_fix_time_between(1,:))./sqrt(length(post_files))])
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
box off
ylabel('# of Fixations')
ylim([10 15])


subplot(1,2,2)
bar([mean(pre_nov_set_means_num_fix_time_between(2,:))  mean(post_nov_set_means_num_fix_time_between(2,:)); ...
    mean(pre_rep_set_means_num_fix_time_between(2,:)) mean(post_rep_set_means_num_fix_time_between(2,:))])
hold on
errorb([mean(pre_nov_set_means_num_fix_time_between(2,:))  mean(post_nov_set_means_num_fix_time_between(2,:)); ...
    mean(pre_rep_set_means_num_fix_time_between(2,:)) mean(post_rep_set_means_num_fix_time_between(2,:))],...
    [std(pre_nov_set_means_num_fix_time_between(2,:))./sqrt(length(pre_files)), ...
    std(post_nov_set_means_num_fix_time_between(2,:))./sqrt(length(post_files));...
    std(pre_rep_set_means_num_fix_time_between(2,:))./sqrt(length(pre_files)), ...
    std(post_rep_set_means_num_fix_time_between(2,:))./sqrt(length(post_files))])
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
box off
ylabel('Time (ms)')
ylim([2500 4000])
legend('Pre','Post')

subtitle([pre_files{1}(1:2) ': IOR, # of fixations and amount of time between prior and return fixations'])