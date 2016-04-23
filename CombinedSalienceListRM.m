function CombinedSalienceListRM(cortex_files,data_dir,figure_dir)
%written by Seth Konig 9/8/15
%code combines data for salience at fixation locations across image sets
%for the ListRM task. Data for individual sets comes from
%CalculateSalienceatFixationsListRM.m

all_novel = NaN(length(cortex_files),50);
all_repeat = NaN(length(cortex_files),50);
all_random = NaN(length(cortex_files),50);
% all_random = [];
median_novel = NaN(1,length(cortex_files));
median_repeat = NaN(1,length(cortex_files));

for file = 1:length(cortex_files)
    load([data_dir cortex_files{file}(1:8) '_' cortex_files{file}(end) '-Salience.mat'])
    
    all_novel(file,:) = nanmean(novel_sal);
    all_repeat(file,:) = nanmean(repeat_sal);
%     all_random = [all_random;random_sal];
	all_random(file,:) = nanmean(random_sal);
    
    temp = sum(~isnan(novel_sal),2);
    temp(temp == 0) = [];%if no fixations don't count
    median_novel(file) = median(temp);
    
    temp = sum(~isnan(repeat_sal),2);
    temp(temp == 0) = [];%if no fixations don't count
    median_repeat(file) = median(temp);
end

%trim down extra fixations
median_novel = median(median_novel);
median_repeat = median(median_repeat);
all_novel = all_novel(:,1:median_novel);
all_repeat = all_repeat(:,1:median_repeat);
all_random= all_random(:,1:median_novel);
[~,~,ci] = ztest(all_random(1:end),mean(all_novel(:,1)),nanstd(all_random(1:end)));

p_fix = NaN(1,size(all_repeat,2));
for fix = 1:size(all_repeat,2)
    [~,p_fix(fix)] = ttest2(all_novel(:,fix),all_repeat(:,fix));
end

figure
hold on
errorbar(mean(all_novel),std(all_novel)./sqrt(length(cortex_files)));
errorbar(mean(all_repeat),std(all_repeat)./sqrt(length(cortex_files)),'r');
errorbar(mean(all_random),std(all_random)./sqrt(length(cortex_files)),'g');
plot([0 median_novel+1],[ci(2) ci(2)],'k--')
hold off
legend('Novel','Repeat','Random Locations')
xlabel('Ordinal Fixation #')
ylabel('Normalized Salience')
title(cortex_files{1}(1:2));