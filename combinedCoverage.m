function combinedCoverage(data_dir,cortex_files)
% Written by Seth Konig
%function combines coverage esitmates across multipe sessions. Code
%importst -Coverage.mat data. 

all_nov_coverage = [];
all_rep_coverage = [];
for file = 1:length(cortex_files)
    load([data_dir cortex_files{file}(1:8) '_' cortex_files{file}(end) '-Coverage.mat']);
    all_nov_coverage(file) = 100*nanmean(nov_coverage);
    all_rep_coverage(file) = 100*nanmean(rep_coverage);
end

[~,p_coverage] = ttest2(all_nov_coverage,all_rep_coverage);

figure
hold on
bar([mean(all_nov_coverage) mean(all_rep_coverage)])
errorb([mean(all_nov_coverage) mean(all_rep_coverage)],...
    [std(all_nov_coverage) std(all_rep_coverage)]./sqrt(length(cortex_files)))
if p_coverage < 0.05
   plot(1.5,25,'k*') 
end
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('% Coverage')
ylim([20 25])
title([cortex_files{1}(1:2) ':  ListRM % Image Covered'])