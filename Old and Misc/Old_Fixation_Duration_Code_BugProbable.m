%---Calculate Pupil Data for Novel and Repeat Images---%
all_nov_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
all_rep_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
% then added 100 samples for 500 ms fixation on center cross

%---Calculate Fixation Durations and Saccade Ampliutdes for Novel and Repeat Images---%
all_nov_durs = NaN(length(cortex_files),50);%
all_rep_durs = NaN(length(cortex_files),50);%
all_nov_sacamps = NaN(length(cortex_files),50);%
all_rep_sacamps = NaN(length(cortex_files),50);%

%---Calculate the number of fixations per img---%
all_num_nov_fix = NaN(length(cortex_files),1);
all_num_rep_fix = NaN(length(cortex_files),1);

%all_pupil_vals = [];
for ls =1:length(cortex_files)
    load([data_dir,cortex_files{ls}(1:8) '_' cortex_files{ls}(end) '-fixation.mat']);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Determine Which Image Trials are "good"---%
    
    viewed = zeros(2,96);
    img_dur = zeros(2,96);
    for trial = 1:length(per);
        [r,c] = find(imgs == per(trial).cnd-1000);
        if viewed(r,c) ~= 0
            error('This should be the first viewing of this image, but it is not!')
        else
            viewed(r,c) = trial;
            img_dur(r,c) = per(trial).alltim(per(trial).allval == 24)-...
                per(trial).alltim(per(trial).allval == 23);
        end
    end
    
    %1st image trial was displayed with the wrong timing file, then when fixed
    %started the task over again. The first image therefore is not novel and thus
    %we will not be analyzing data for this image on this session.
    if strcmpi('PW150416_2',cortex_files{ls}(1:10))
        viewed(:,1) = 0;
        imgdur(:,1) = 0;
    end
    
    great_views = viewed;
    great_views(viewed & img_dur > 14000) = 0; %1000 ms buffer for accidental lookaways and cortex timing
    
    great_views(:,great_views(1,:) == 0 | great_views(2,:) == 0) = [];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Analyze the Data for the "good" trials---%
    nov_pupil = NaN(size(great_views,2),1500);
    rep_pupil = NaN(size(great_views,2),1500);
    
    novel_fixdurs = NaN(size(great_views,2),50);
    repeat_fixdurs = NaN(size(great_views,2),50);
    novel_sacamps = NaN(size(great_views,2),50);
    repeat_sacamps = NaN(size(great_views,2),50);
    
    nov_num_fix = NaN(size(great_views,2),1);
    rep_num_fix = NaN(size(great_views,2),1);
    for novrep = 1:size(great_views,2);
        
        
        %---Import Important Variables---%
        nov_allval = per(great_views(1,novrep)).allval;
        nov_alltim = per(great_views(1,novrep)).alltim;
        nov_img_on = nov_alltim(nov_allval == 23)-nov_alltim(nov_allval == 100);%image on relative to eye data start
        nov_img_off = nov_alltim(nov_allval == 24)-nov_alltim(nov_allval == 100);%image on relative to eye data start
        
        nov_fixations = fixationstats{great_views(1,novrep)}.fixations;
        nov_fixationtimes = fixationstats{great_views(1,novrep)}.fixationtimes;
        nov_saccadetimes = fixationstats{great_views(1,novrep)}.saccadetimes;
        nov_xy = fixationstats{great_views(1,novrep)}.XY;
        
        nov_trial_pupil = pupildata{great_views(1,novrep)};
        
        
        rep_allval = per(great_views(2,novrep)).allval;
        rep_alltim = per(great_views(2,novrep)).alltim;
        rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        
        rep_fixations = fixationstats{great_views(2,novrep)}.fixations;
        rep_fixationtimes = fixationstats{great_views(2,novrep)}.fixationtimes;
        rep_saccadetimes = fixationstats{great_views(2,novrep)}.saccadetimes;
        rep_xy = fixationstats{great_views(2,novrep)}.XY;
        
        rep_trial_pupil = pupildata{great_views(2,novrep)};
        
        %---Remove Fixations before and after image was on---%
        
        nov_pre_img_fixations = find(nov_fixationtimes(1,:) <= nov_img_on);
        nov_post_img_fixations = find(nov_fixationtimes(2,:) > nov_img_off);
        nov_pre_img_saccades = find(nov_saccadetimes(1,:) <= nov_img_on);
        nov_post_img_saccades = find(nov_saccadetimes(2,:) > nov_img_off);
        nov_fixations(:,nov_post_img_fixations) = [];
        nov_fixationtimes(:,nov_post_img_fixations) = [];
        nov_fixations(:,nov_pre_img_fixations) = [];
        nov_fixationtimes(:,nov_pre_img_fixations) = [];
        nov_saccadetimes(:,nov_post_img_saccades) = [];
        nov_saccadetimes(:,nov_pre_img_saccades) = [];
        
        rep_pre_img_fixations = find(rep_fixationtimes(1,:) <= rep_img_on);
        rep_post_img_fixations = find(rep_fixationtimes(2,:) > rep_img_off);
        rep_pre_img_saccades = find(rep_saccadetimes(1,:) <= rep_img_on);
        rep_post_img_saccades = find(rep_saccadetimes(2,:) > rep_img_off);
        rep_fixations(:,rep_post_img_fixations) = [];
        rep_fixationtimes(:,rep_post_img_fixations) = [];
        rep_fixations(:,rep_pre_img_fixations) = [];
        rep_fixationtimes(:,rep_pre_img_fixations) = [];
        rep_saccadetimes(:,rep_post_img_saccades) = [];
        rep_saccadetimes(:,rep_pre_img_saccades) = [];
        
        
        %---Remove Central Fixations that continued to Occur once image was on---%
%         first_out = find((nov_fixations(1,:) > fixwin | nov_fixations(1,:) < -fixwin) | ...
%             (nov_fixations(2,:) > fixwin | nov_fixations(2,:) < -fixwin));
%         first_out = first_out(1);
%         %         if first_out > 2
%         %             disp('First out > 2')
%         %         end
%         nov_fixations(:,1:first_out-1) = [];
%         nov_fixationtimes(:,1:first_out-1) = [];
%         
%         first_out = find((rep_fixations(1,:) > fixwin | rep_fixations(1,:) < -fixwin) | ...
%             (rep_fixations(2,:) > fixwin | rep_fixations(2,:) < -fixwin));
%         first_out = first_out(1);
%         %         if first_out > 2
%         %             disp('First out > 2')
%         %         end
%         rep_fixations(:,1:first_out-1) = [];
%         rep_fixationtimes(:,1:first_out-1) = [];
        
        %---Calculate the Number of Fixations per image--%
        nov_num_fix(novrep) =size(nov_fixations,2);
        rep_num_fix(novrep) =size(rep_fixations,2);
        
        %probably shouldn't happen but just in case
        if size(nov_fixationtimes,2) > 50;
            nov_fixationtimes(:,51:end) = [];
        end
        if size(rep_fixationtimes,2) > 50;
            rep_fixationtimes(:,51:end) = [];
        end
        
        nov_fixationtimes(:,1) = [];
        rep_fixationtimes(:,2) = [];
        %---Calculate Fixation Durations---%
        novfixdurs = (nov_fixationtimes(2,:) - nov_fixationtimes(1,:))+ 1;
        novel_fixdurs(novrep,1:length(novfixdurs)) = novfixdurs;
        
        repfixdurs = (rep_fixationtimes(2,:) - rep_fixationtimes(1,:))+ 1;
        repeat_fixdurs(novrep,1:length(repfixdurs)) = repfixdurs;
        
        %---Calculate Saccade Amplitudes---%
        for s = 1:size(nov_saccadetimes,2);
            sacx = nov_xy(1,nov_saccadetimes(2,s))-nov_xy(1,nov_saccadetimes(1,s));
            sacy = nov_xy(2,nov_saccadetimes(2,s))-nov_xy(2,nov_saccadetimes(1,s));
            novel_sacamps(novrep,s) =sqrt(sacx^2+sacy^2);
        end
        for s = 1:size(rep_saccadetimes,2);
            sacx = rep_xy(1,rep_saccadetimes(2,s))-rep_xy(1,rep_saccadetimes(1,s));
            sacy = rep_xy(2,rep_saccadetimes(2,s))-rep_xy(2,rep_saccadetimes(1,s));
            repeat_sacamps(novrep,s) =sqrt(sacx^2+sacy^2);
        end
        
        %---Calculate Pupil Diameter For Image Duration---%
        %all_pupil_vals = [all_pupil_vals;nov_trial_pupil;rep_trial_pupil];
        
        nov_img_on = round(nov_img_on/5);
        rep_img_on = round(rep_img_on/5);
        
        nov_trial_pupil = nov_trial_pupil(nov_img_on-100:nov_img_on+1399); %only want 1st 7 seconds+ 500 ms fixation on center cross
        rep_trial_pupil = rep_trial_pupil(rep_img_on-100:rep_img_on+1399);%only want 1st 7 seconds + 500 ms fixation on center cross
        
        %remove blinks artifacts? estimate from histogram should replace
        %with blink detection
        nov_trial_pupil(nov_trial_pupil < -1250) = NaN;
        rep_trial_pupil(rep_trial_pupil < -1250) = NaN;
        
        %normalize pupil data on every trial but starting value. May need
        %to change to import pupil data before image turns on
        nov_trial_pupil = nov_trial_pupil./abs(mean(nov_trial_pupil(81:100)))+2;
        rep_trial_pupil = rep_trial_pupil./abs(mean(rep_trial_pupil(81:100)))+2;
        
        nov_pupil(novrep,:) = nov_trial_pupil;
        rep_pupil(novrep,:) = rep_trial_pupil;
        
    end
    
    %---Average across all novel/rep images within a session---%
    %only take up to the median number of fixations for novel and repeat
    %presetnations sepperately
    nov_median_num_fix = floor(median(sum(~isnan(novel_fixdurs'))));
    nov_median_num_sac = floor(median(sum(~isnan(novel_sacamps'))));
    all_nov_durs(ls,1:nov_median_num_fix) = nanmean(novel_fixdurs(:,1:nov_median_num_fix));
    all_num_nov_fix(ls) = mean(nov_num_fix);
    all_nov_sacamps(ls,1:nov_median_num_sac) = nanmean(novel_sacamps(:,1:nov_median_num_sac));
    all_nov_pupil(ls,:) = nanmean(nov_pupil);
    
    rep_median_num_fix = floor(median(sum(~isnan(repeat_fixdurs'))));
    rep_median_num_sac = floor(median(sum(~isnan(repeat_sacamps'))));
    all_rep_durs(ls,1:rep_median_num_fix) = nanmean(repeat_fixdurs(:,1:rep_median_num_fix));
    all_num_rep_fix(ls) = mean(rep_num_fix);
    all_rep_sacamps(ls,1:rep_median_num_sac) = nanmean(repeat_sacamps(:,1:rep_median_num_sac));
    all_rep_pupil(ls,:) = nanmean(rep_pupil);
    
end
%%
%%%%%%%%%%%%%%%%%
%---Plot Data---%

%---Plot Number of Fixations for Novel vs repat images---%
[~,p] = ttest2(all_num_nov_fix,all_num_rep_fix);
figure
hold on
bar([mean(all_num_nov_fix) mean(all_num_rep_fix)])
errorb([mean(all_num_nov_fix) mean(all_num_rep_fix)],...
    [std(all_num_nov_fix) std(all_num_rep_fix)]...
    ./[sqrt(length(all_num_nov_fix)) sqrt(length(all_num_rep_fix))])
if p < 0.05
    plot(1.5,25,'*k','markersize',5)
end
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Number of fixiations')
ylim([22 25])

%---Plot Fixation Durations By Fixation Number---%
nov_median_num_fix = ceil(median(sum(~isnan(all_nov_durs'))));
rep_median_num_fix = ceil(median(sum(~isnan(all_rep_durs'))));
all_nov_durs = all_nov_durs(:,1:nov_median_num_fix);
all_rep_durs = all_rep_durs(:,1:rep_median_num_fix);

figure
hold on
errorbar(nanmean(all_nov_durs),nanstd(all_nov_durs)...
    ./sqrt(sum(~isnan(all_nov_durs))),'b')
errorbar(nanmean(all_rep_durs),nanstd(all_rep_durs)...
    ./sqrt(sum(~isnan(all_rep_durs))),'r')
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
legend('Novel','Repeat')

%---Plot Saccade Ampltiudes By Saccades Number---%
nov_median_num_sac = ceil(median(sum(~isnan(all_nov_sacamps'))));
rep_median_num_sac = ceil(median(sum(~isnan(all_rep_sacamps'))));
all_nov_sacamps = all_nov_sacamps(:,1:nov_median_num_sac);
all_rep_sacamps = all_rep_sacamps(:,1:rep_median_num_sac);

figure
hold on
errorbar(nanmean(all_nov_sacamps),nanstd(all_nov_sacamps)...
    ./sqrt(sum(~isnan(all_nov_sacamps))),'b')
errorbar(nanmean(all_rep_sacamps),nanstd(all_rep_sacamps)...
    ./sqrt(sum(~isnan(all_rep_sacamps))),'r')
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
legend('Novel','Repeat')

%---Plot Pupil Diameter Over Time---%
t = 1:5:7500;
figure
hold on
dofill(t,all_nov_pupil/1000,'blue',1,60)
dofill(t,all_rep_pupil/1000,'red',1,60)
yl = ylim;
plot([500 500],[0.75 1.1],'--k')
hold off
ylim(yl)
xlabel('Time from Image Onset (ms)')
ylabel('Normalize Pupil Data (a.u.)')
set(gca,'Xtick',0:500:7500)
set(gca,'XtickLabel',num2cell([0:500:7500]-500))
legend('Novel Images','Repeat Images')
title('Normalized Pupil Horizontal Diameter for Images Displayed < 10 seconds')