%use from debug mode in combinedEyeMovementStats
%---Plot Fixation Durations By Fixation Number---%
nov_mean_fix_dur = all_novel_fix_dur;
rep_mean_fix_dur =all_repeat_fix_dur;
for monk = 1:4
    these_sets = find(which_monkey == monk);
    first_dur = nanmean(nov_mean_fix_dur(these_sets,1)); 
    nov_mean_fix_dur(these_sets,:) = nov_mean_fix_dur(these_sets,:)./first_dur;
    rep_mean_fix_dur(these_sets,:) = rep_mean_fix_dur(these_sets,:)./first_dur;
end
    

figure
hold on
errorbar(nanmean(nov_mean_fix_dur(:,1:22)),nanstd(nov_mean_fix_dur(:,1:22))...
    ./sqrt(sum(~isnan(nov_mean_fix_dur(:,1:22)))),'b')
errorbar(nanmean(rep_mean_fix_dur(:,1:22)),nanstd(rep_mean_fix_dur(:,1:22))...
    ./sqrt(sum(~isnan(rep_mean_fix_dur(:,1:22)))),'r')
hold off
xlabel('Ordinal Fixation #')
ylabel('Normalized Fixation Duration')
legend('Novel','Repeat')
title(cortex_files{1}(1:2))

%%
%---Plot Saccade Ampltiudes By Saccades Number---%
nov_mean_sac_amp = all_novel_sac_amp;
rep_mean_sac_amp =all_repeat_sac_amp;
for monk = 1:4
    these_sets = find(which_monkey == monk);
    first_dur = nanmean(nov_mean_sac_amp(these_sets,1)); 
    nov_mean_sac_amp(these_sets,:) = nov_mean_sac_amp(these_sets,:)./first_dur;
    rep_mean_sac_amp(these_sets,:) = rep_mean_sac_amp(these_sets,:)./first_dur;
end
    

figure
hold on
errorbar(nanmean(nov_mean_sac_amp(:,1:22)),nanstd(nov_mean_sac_amp(:,1:22))...
    ./sqrt(sum(~isnan(nov_mean_sac_amp(:,1:22)))),'b')
errorbar(nanmean(rep_mean_sac_amp(:,1:22)),nanstd(rep_mean_sac_amp(:,1:22))...
    ./sqrt(sum(~isnan(rep_mean_sac_amp(:,1:22)))),'r')
hold off
xlabel('Ordinal Saccade #')
ylabel('Normalized Saccade Amplitude')
legend('Novel','Repeat')
title(cortex_files{1}(1:2))

