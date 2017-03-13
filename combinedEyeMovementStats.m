function combinedEyeMovementStats(data_dir,cortex_files)
%function combines viewing behavior across multiple sessions
% written by Seth Koenig
% 1) calculate fixation durations
% 2) calculate saccade amplitudes
% 3) calculate pupil diameter 

first_fix = cell(1,length(cortex_files));
first_fix_location = cell(2,length(cortex_files));
all_blinks = cell(1,2);

all_novel_fix_durt = NaN(length(cortex_files),7000); %fixation durations
all_repeat_fix_durt = NaN(length(cortex_files),7000); %fixation durations


all_nov_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
all_rep_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
all_novel_fix_dur = NaN(length(cortex_files),50); %fixation durations
all_repeat_fix_dur = NaN(length(cortex_files),50); %fixation durations
all_novel_sac_amp = NaN(length(cortex_files),50); %saccade amplitudes
all_repeat_sac_amp = NaN(length(cortex_files),50); %saccade amplitudes
nov_time_out = NaN(length(cortex_files),96); %time when first left
rep_time_out = NaN(length(cortex_files),96); %time when first left
nov_first_dur = NaN(length(cortex_files),96);
rep_first_dur = NaN(length(cortex_files),96);
which_monkey = [];
for file = 1:length(cortex_files)
    load([data_dir cortex_files{file}(1:8) '_' cortex_files{file}(end) '-fixation.mat'])
    if strcmpi( cortex_files{file}(1:2),'MF')
          which_monkey = [which_monkey 1];
    elseif strcmpi( cortex_files{file}(1:2),'PW')
        which_monkey = [which_monkey 2];
    elseif strcmpi( cortex_files{file}(1:2),'TO')
        which_monkey = [which_monkey 3];
    elseif strcmpi( cortex_files{file}(1:2),'RR')
        which_monkey = [which_monkey 4];
    end
    
    first_fix{file} = NaN(2,96);
    first_fix_location{1,file} = NaN(2,96);
    first_fix_location{2,file} = NaN(2,96);
      
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
    imgnum = 1:96;
    
    %1st image trial was displayed with the wrong timing file, then when fixed
    %started the task over again. The first image therefore is not novel and thus
    %we will not be analyzing data for this image on this session.
    if strcmpi('PW150416_2',cortex_files{file}(1:10))
        viewed(:,1) = 0;
        imgdur(:,1) = 0;
    end
        
    %Let's not even look at images in which the monkey looked away more than the
    %image was displayed for. Also should take care of EOG overflow trials
    viewed(viewed & img_dur > 14000) = 0; %1000 ms buffer for accidental lookaways and cortex timing
    imgnum(viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    viewed(:,viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    
    nov_trial_by_trial_pupil = NaN(96,1500); %each trials pupil data
    rep_trial_by_trial_pupil = NaN(96,1500); %each trials pupil data
    nov_trial_by_trial_fixdurs = NaN(96,50); %each trials fixation durations
    rep_trial_by_trial_fixdurs = NaN(96,50); %each trials fixation durations
    nov_trial_by_trial_sacamps = NaN(96,50); %each trials saccade amplitudes
    rep_trial_by_trial_sacamps = NaN(96,50); %each trials fixation amplitudes
    
    nov_trial_by_trial_fixdurts = NaN(96,7000); %each trials fixation durations
    rep_trial_by_trial_fixdurts = NaN(96,7000); %each trials fixation durations
    
    for novrep = 1:size(viewed,2);
        
        
        %---Grab Important Vairables---%
        %for novel images
        nov_allval = per(viewed(1,novrep)).allval;
        nov_alltim = per(viewed(1,novrep)).alltim;
        nov_img_on = nov_alltim(nov_allval == 23)-nov_alltim(nov_allval == 100);%image on relative to eye data start
        nov_img_off = nov_alltim(nov_allval == 24)-nov_alltim(nov_allval == 100);%image on relative to eye data start
        
        nov_x = fixationstats{viewed(1,novrep)}.XY(1,nov_img_on:nov_img_off);
        nov_y = fixationstats{viewed(1,novrep)}.XY(2,nov_img_on:nov_img_off);
        nov_pupil = pupildata{viewed(1,novrep)}(round(nov_img_on/5):round(nov_img_off/5));
        
        nov_fix = fixationstats{viewed(1,novrep)}.fixations;
        nov_fixtimes = fixationstats{viewed(1,novrep)}.fixationtimes;
        nov_sactimes = fixationstats{viewed(1,novrep)}.saccadetimes;
        
        pre_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
        post_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
        pre_img_saccades = find(nov_sactimes(1,:) <= nov_img_on);
        post_img_saccades = find(nov_sactimes(2,:) > nov_img_off);
        
        nov_fix(:,post_img_fix) = [];
        nov_fix(:,pre_img_fix) = [];
        nov_fixtimes(:,post_img_fix) = [];
        nov_fixtimes(:,pre_img_fix ) = [];
        nov_sactimes(:,post_img_saccades) = [];
        nov_sactimes(:,pre_img_saccades) = [];
        
        
        %for repeat images
        rep_allval = per(viewed(2,novrep)).allval;
        rep_alltim = per(viewed(2,novrep)).alltim;
        rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        
        rep_x = fixationstats{viewed(2,novrep)}.XY(1,rep_img_on:rep_img_off);
        rep_y = fixationstats{viewed(2,novrep)}.XY(2,rep_img_on:rep_img_off);
        rep_pupil = pupildata{viewed(2,novrep)}(round(rep_img_on/5):round(rep_img_off/5));
        
        rep_fix = fixationstats{viewed(2,novrep)}.fixations;
        rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
        rep_sactimes = fixationstats{viewed(2,novrep)}.saccadetimes;
        
        pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        pre_img_saccades = find(rep_sactimes(1,:) <= rep_img_on);
        post_img_saccades = find(rep_sactimes(2,:) > rep_img_off);
        
        rep_fix(:,post_img_fix) = [];
        rep_fix(:,pre_img_fix) = [];
        rep_fixtimes(:,post_img_fix) = [];
        rep_fixtimes(:,pre_img_fix) = [];
        rep_sactimes(:,post_img_saccades) = [];
        rep_sactimes(:,pre_img_saccades) = [];
        
        %---Remove Central Fixations that continued to Occur once image was on---%
        %         nov_fix(:,1) = [];
        %         nov_fixtimes(:,1) = [];
        %         rep_fix(:,1) = [];
        %         rep_fixtimes(:,1) = [];
        
        %         first_out = find((nov_fix(1,:) > imageX/2+fixwin*24 | nov_fix(1,:) <  imageX/2-24*fixwin) | ...
        %             (nov_fix(2,:) > imageY/2+fixwin*24 | nov_fix(2,:) < imageY/2-fixwin*24));
        %         nov_fix(:,1:first_out-1) = [];
        %         nov_fixtimes(:,1:first_out-1) = [];
        %
        %         first_out = find((rep_fix(1,:) > imageX/2+fixwin*24 | rep_fix(1,:) <  imageX/2-24*fixwin) | ...
        %             (rep_fix(2,:) > imageY/2+fixwin*24 | rep_fix(2,:) < imageY/2-fixwin*24));
        %         rep_fix(:,1:first_out-1) = [];
        %         rep_fixtimes(:,1:first_out-1) = [];
        
        %remove eye movements that start within the 1st 150 ms of the image
        %turning on
%         nov_fix(:,nov_fixtimes(1,:) < 150+nov_img_on) = [];
%         rep_fix(:,rep_fixtimes(1,:) < 150+rep_img_on) = [];
%         nov_fixtimes(:,nov_fixtimes(1,:) < 150+nov_img_on) = [];
%         rep_fixtimes(:,rep_fixtimes(1,:) < 150+rep_img_on) = [];
%         nov_sactimes(:,nov_sactimes(1,:) < 150+nov_img_on) = [];
%         rep_sactimes(:,rep_sactimes(1,:) < 150+rep_img_on) = [];
        
        %---Find When the monkey Blinked,when monkey looked away---%
        %pupil values at 0 diameter
        
        [nov_blink_ind,nov_nan_ind,nov_time_out(file,novrep)] = ...
            findbad_eye_ind(nov_pupil,nov_x,1000);
        [rep_blink_ind,rep_nan_ind,rep_time_out(file,novrep)] = ...
            findbad_eye_ind(rep_pupil,rep_x,1000);
        
        if ~isempty(nov_blink_ind)
            all_blinks{1} = [all_blinks{1} nov_blink_ind(1)];
        end
        if ~isempty(rep_blink_ind)
            all_blinks{2} = [all_blinks{2} rep_blink_ind(1)];
        end
        
        
        %---Get Pupil Data over time---%
        nov_pupil = pupildata{viewed(1,novrep)}((round(nov_img_on/5)-100):round(nov_img_off/5)); %grab 500 ms before img on too
        rep_pupil = pupildata{viewed(2,novrep)}((round(rep_img_on/5)-100):round(rep_img_off/5)); %grab 500 ms before img on too
        
        %remove blinks from pupil data
        if ~isempty(nov_blink_ind)
            for b = 1:size(nov_blink_ind,1)
                ind = nov_blink_ind(b,:);
                ind(ind == 0) = [];
                nov_pupil(ind+100) = NaN;
            end
        end
        if ~isempty(rep_blink_ind)
            for b = 1:size(rep_blink_ind,1)
                ind = rep_blink_ind(b,:);
                ind(ind == 0) = [];
                rep_pupil(ind+100) = NaN;
            end
        end
        
        %remove pupil data when not looking at picture
        if ~isempty(nov_nan_ind)
            for b = 1:size(nov_nan_ind,1)
                ind = nov_nan_ind(b,:);
                ind(ind == 0) = [];
                nov_pupil(round(ind(1)/5):round(ind(end)/5)) = NaN;
            end
        end
        if ~isempty(rep_nan_ind)
            for b = 1:size(rep_nan_ind,1)
                ind = rep_nan_ind(b,:);
                ind(ind == 0) = [];
                rep_pupil(round(ind(1)/5):round(ind(end)/5)) = NaN;
            end
        end
        
        %remove data after looked away too much
        if ~isnan(nov_time_out(file,novrep))
            nov_pupil(round(nov_time_out(file,novrep)/5):end) = NaN;
        end
        if ~isnan(rep_time_out(file,novrep))
            rep_pupil(round(rep_time_out(file,novrep)/5):end) = NaN;
        end
        
        %normalize to 100 ms of prestimulus levels
        nov_pupil = nov_pupil./abs(mean(nov_pupil(81:100)))+2;
        rep_pupil = rep_pupil./abs(mean(rep_pupil(81:100)))+2;
        
        %only look at the 1st 7 seconds the images is up
        nov_pupil = nov_pupil(1:1500);
        rep_pupil = rep_pupil(1:1500);
        
        nov_trial_by_trial_pupil(novrep,:) = nov_pupil;
        rep_trial_by_trial_pupil(novrep,:) = rep_pupil;
        
        %---Calculate Fixation Durations and Saccade Amplitudes---%
        
        if ~isnan(nov_time_out(file,novrep))
            post_attention_fix = find(nov_fixtimes(1,:) > nov_time_out(file,novrep));
            nov_fixtimes(:,post_attention_fix) = [];
        end
        
        if ~isnan(nov_time_out(file,novrep))
            post_attention_fix = find(rep_fixtimes(1,:) > rep_time_out(file,novrep));
            rep_fixtimes(:,post_attention_fix) = [];
        end
        
        %fixation duration over time
        for f = 1:size(nov_fixtimes,2)
            dur =  nov_fixtimes(2,f)-nov_fixtimes(1,f)+1;
            nov_trial_by_trial_fixdurts(novrep,[nov_fixtimes(1,f):nov_fixtimes(2,f)]-nov_img_on) = dur; %each trials fixation durations
        end
        for f = 1:size(rep_fixtimes,2)
            dur =  rep_fixtimes(2,f)-rep_fixtimes(1,f)+1;
            rep_trial_by_trial_fixdurts(novrep,[rep_fixtimes(1,f):rep_fixtimes(2,f)]-rep_img_on) = dur; %each trials fixation durations
        end
        
        %fixation duration by ordinal fixation number
        nov_trial_by_trial_fixdurs(novrep,1:length(nov_fixtimes)) = nov_fixtimes(2,:)-nov_fixtimes(1,:)+1;
        rep_trial_by_trial_fixdurs(novrep,1:length(rep_fixtimes)) = rep_fixtimes(2,:)-rep_fixtimes(1,:)+1;
        
        if ~isempty(nov_fixtimes) && ~isempty(rep_fixtimes)
            first_fix{file}(1,imgnum(novrep)) = nov_fixtimes(2,1)-nov_fixtimes(1,1)+1;
            first_fix{file}(2,imgnum(novrep)) = rep_fixtimes(2,1)-rep_fixtimes(1,1)+1;
            
            first_fix_location{1,file}(:,imgnum(novrep)) = nov_fix(:,1);
            first_fix_location{2,file}(:,imgnum(novrep)) = rep_fix(:,1);
        end
        
        nov_x = fixationstats{viewed(1,novrep)}.XY(1,:);
        nov_y = fixationstats{viewed(1,novrep)}.XY(2,:);
        
        for sac = 1:size(nov_sactimes,2)
            sacx = nov_x(nov_sactimes(2,sac)) - nov_x(nov_sactimes(1,sac));
            sacy = nov_y(nov_sactimes(2,sac)) - nov_y(nov_sactimes(1,sac));
            nov_trial_by_trial_sacamps(novrep,sac) = sqrt(sacx^2+sacy^2);
        end
        
        rep_x = fixationstats{viewed(2,novrep)}.XY(1,:);
        rep_y = fixationstats{viewed(2,novrep)}.XY(2,:);
        for sac = 1:size(rep_sactimes,2)
            sacx = rep_x(rep_sactimes(2,sac)) - rep_x(rep_sactimes(1,sac));
            sacy = rep_y(rep_sactimes(2,sac)) - rep_y(rep_sactimes(1,sac));
            rep_trial_by_trial_sacamps(novrep,sac) = sqrt(sacx^2+sacy^2);
        end
    end
    
    all_nov_pupil(file,:) = nanmean(nov_trial_by_trial_pupil);
    all_rep_pupil(file,:) = nanmean(rep_trial_by_trial_pupil);
    
    nov_first_dur(file,:) = nov_trial_by_trial_fixdurs(:,1)';
    rep_first_dur(file,:) = rep_trial_by_trial_fixdurs(:,1)';
    
    nov_median_num_fix = sum(~isnan(nov_trial_by_trial_fixdurs'));
    nov_median_num_fix(nov_median_num_fix == 0) = [];
    nov_median_num_fix = round(median(nov_median_num_fix));
    all_novel_fix_dur(file,1:nov_median_num_fix) = ...
        nanmean(nov_trial_by_trial_fixdurs(:,1:nov_median_num_fix));
    
    rep_median_num_fix = sum(~isnan(rep_trial_by_trial_fixdurs'));
    rep_median_num_fix(rep_median_num_fix == 0) = [];
    rep_median_num_fix = round(median(rep_median_num_fix));
    all_repeat_fix_dur(file,1:rep_median_num_fix) = ...
        nanmean(rep_trial_by_trial_fixdurs(:,1:rep_median_num_fix));
    
    
    nov_median_num_fix = sum(~isnan(nov_trial_by_trial_sacamps'));
    nov_median_num_fix(nov_median_num_fix == 0) = [];
    nov_median_num_fix = round(median(nov_median_num_fix));
    all_novel_sac_amp(file,1:nov_median_num_fix) = ...
        nanmean(nov_trial_by_trial_sacamps(:,1:nov_median_num_fix));
    
    rep_median_num_fix = sum(~isnan(rep_trial_by_trial_sacamps'));
    rep_median_num_fix(rep_median_num_fix == 0) = [];
    rep_median_num_fix = round(median(rep_median_num_fix));
    all_repeat_sac_amp(file,1:rep_median_num_fix) = ...
        nanmean(rep_trial_by_trial_sacamps(:,1:rep_median_num_fix));
    
    all_novel_fix_durt(file,:) = nanmean(nov_trial_by_trial_fixdurts(:,1:7000)); %fixation durations
    all_repeat_fix_durt(file,:) = nanmean(rep_trial_by_trial_fixdurts(:,1:7000)); %fixation durations
end
%%
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
title([cortex_files{1}(1:2) ' : Normalized Pupil Horizontal Diameter'])

%%
%---Plot Fixation durations over time---%
t = 1:7000;
figure
hold on
dofill(t,all_novel_fix_durt(:,1:end)/1000,'blue',1,180)
dofill(t,all_repeat_fix_durt(:,1:end)/1000,'red',1,180)
hold off
xlabel('Time from Image Onset (ms)')
ylabel('Fixation Duration (ms)')
set(gca,'Xtick',0:500:7000)
set(gca,'XtickLabel',num2cell([0:500:7000]))
legend('Novel Images','Repeat Images')
title([cortex_files{1}(1:2) ' : Fixation Duration over time'])
%%
%---Plot Fixation Durations By Fixation Number---%
nov_median_num_fix = ceil(median(sum(~isnan(all_novel_fix_dur'))));
rep_median_num_fix = ceil(median(sum(~isnan(all_repeat_fix_dur'))));
all_novel_fix_dur = all_novel_fix_dur(:,1:nov_median_num_fix);
all_repeat_fix_dur = all_repeat_fix_dur(:,1:rep_median_num_fix);

figure
hold on
errorbar(nanmean(all_novel_fix_dur),nanstd(all_novel_fix_dur)...
    ./sqrt(sum(~isnan(all_novel_fix_dur))),'b')
errorbar(nanmean(all_repeat_fix_dur),nanstd(all_repeat_fix_dur)...
    ./sqrt(sum(~isnan(all_repeat_fix_dur))),'r')
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
legend('Novel','Repeat')
title(cortex_files{1}(1:2))
%%
%---Plot Saccade Ampltiudes By Saccades Number---%
nov_median_num_sac = ceil(median(sum(~isnan(all_novel_sac_amp'))));
rep_median_num_sac = ceil(median(sum(~isnan(all_repeat_sac_amp'))));
all_novel_sac_amp = all_novel_sac_amp(:,1:nov_median_num_sac)/24;
all_repeat_sac_amp = all_repeat_sac_amp(:,1:rep_median_num_sac)/24;

figure
hold on
errorbar(nanmean(all_novel_sac_amp),nanstd(all_novel_sac_amp)...
    ./sqrt(sum(~isnan(all_novel_sac_amp))),'b')
errorbar(nanmean(all_repeat_sac_amp),nanstd(all_repeat_sac_amp)...
    ./sqrt(sum(~isnan(all_repeat_sac_amp))),'r')
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
legend('Novel','Repeat')
title(cortex_files{1}(1:2))