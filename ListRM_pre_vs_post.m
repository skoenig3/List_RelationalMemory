%written by Seth Konig 6/8/16
clear,clc
fixwin = 3.5; %size of the fixation window/2. Was a width of 7 dva
imageX = 800;
imageY = 600;

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Cortex Data\';
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Image Sets\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Figures\';

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Eye movement Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Pre-Lesion Data---%
pre_nov_pupil = NaN(length(pre_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
pre_rep_pupil = NaN(length(pre_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
pre_novel_fix_dur = NaN(length(pre_files),50); %fixation durations
pre_repeat_fix_dur = NaN(length(pre_files),50); %fixation durations
pre_novel_sac_amp = NaN(length(pre_files),50); %saccade amplitudes
pre_repeat_sac_amp = NaN(length(pre_files),50); %saccade amplitudes
pre_fix_medians = NaN(2,length(pre_files));
pre_sac_medians = NaN(2,length(pre_files));
for file = 1:length(pre_files)
    load([data_dir pre_files{file}(1:8) '_' pre_files{file}(end) '-fixation.mat'])
    
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
    if strcmpi('PW150416_2',pre_files{file}(1:10))
        viewed(:,1) = 0;
        imgdur(:,1) = 0;
    end
    
    %Let's not even look at images in which the monkey looked away more than the
    %image was displayed for. Also should take care of EOG overflow trials
    viewed(viewed & img_dur > 12000) = 0; %1000 ms buffer for accidental lookaways and cortex timing
    imgnum(viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    viewed(:,viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    
    nov_trial_by_trial_pupil = NaN(96,1500); %each trials pupil data
    rep_trial_by_trial_pupil = NaN(96,1500); %each trials pupil data
    nov_trial_by_trial_fixdurs = NaN(96,50); %each trials fixation durations
    rep_trial_by_trial_fixdurs = NaN(96,50); %each trials fixation durations
    nov_trial_by_trial_sacamps = NaN(96,50); %each trials saccade amplitudes
    rep_trial_by_trial_sacamps = NaN(96,50); %each trials fixation amplitudes
    
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
        
        prep_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
        post_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
        prep_img_saccades = find(nov_sactimes(1,:) <= nov_img_on);
        post_img_saccades = find(nov_sactimes(2,:) > nov_img_off);
        
        nov_fix(:,post_img_fix) = [];
        nov_fix(:,prep_img_fix) = [];
        nov_fixtimes(:,post_img_fix) = [];
        nov_fixtimes(:,prep_img_fix ) = [];
        nov_sactimes(:,post_img_saccades) = [];
        nov_sactimes(:,prep_img_saccades) = [];
        
        
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
        
        prep_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        prep_img_saccades = find(rep_sactimes(1,:) <= rep_img_on);
        post_img_saccades = find(rep_sactimes(2,:) > rep_img_off);
        
        rep_fix(:,post_img_fix) = [];
        rep_fix(:,prep_img_fix) = [];
        rep_fixtimes(:,post_img_fix) = [];
        rep_fixtimes(:,prep_img_fix) = [];
        rep_sactimes(:,post_img_saccades) = [];
        rep_sactimes(:,prep_img_saccades) = [];
        
        %---Find When the monkey Blinked,when monkey looked away---%
        %pupil values at 0 diameter
        
        [nov_blink_ind,nov_nan_ind,nov_time_out(file,novrep)] = ...
            findbad_eye_ind(nov_pupil,nov_x,1000);
        [rep_blink_ind,rep_nan_ind,rep_time_out(file,novrep)] = ...
            findbad_eye_ind(rep_pupil,rep_x,1000);
        
        
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
    
    pre_nov_pupil(file,:) = nanmean(nov_trial_by_trial_pupil);
    pre_rep_pupil(file,:) = nanmean(rep_trial_by_trial_pupil);
    
    nov_median_num_fix = sum(~isnan(nov_trial_by_trial_fixdurs'));
    nov_median_num_fix(nov_median_num_fix == 0) = [];
    nov_median_num_fix = round(median(nov_median_num_fix));
    pre_fix_medians(1,file) = nov_median_num_fix;
    pre_novel_fix_dur(file,1:nov_median_num_fix) = ...
        nanmean(nov_trial_by_trial_fixdurs(:,1:nov_median_num_fix));
    
    rep_median_num_fix = sum(~isnan(rep_trial_by_trial_fixdurs'));
    rep_median_num_fix(rep_median_num_fix == 0) = [];
    rep_median_num_fix = round(median(rep_median_num_fix));
    pre_fix_medians(2,file) = rep_median_num_fix;
    pre_repeat_fix_dur(file,1:rep_median_num_fix) = ...
        nanmean(rep_trial_by_trial_fixdurs(:,1:rep_median_num_fix));
    
    
    nov_median_num_fix = sum(~isnan(nov_trial_by_trial_sacamps'));
    nov_median_num_fix(nov_median_num_fix == 0) = [];
    nov_median_num_fix = round(median(nov_median_num_fix));
    pre_sac_medians(1,file) = nov_median_num_fix;
    pre_novel_sac_amp(file,1:nov_median_num_fix) = ...
        nanmean(nov_trial_by_trial_sacamps(:,1:nov_median_num_fix));
    
    rep_median_num_fix = sum(~isnan(rep_trial_by_trial_sacamps'));
    rep_median_num_fix(rep_median_num_fix == 0) = [];
    rep_median_num_fix = round(median(rep_median_num_fix));
    pre_sac_medians(2,file) = rep_median_num_fix;
    pre_repeat_sac_amp(file,1:rep_median_num_fix) = ...
        nanmean(rep_trial_by_trial_sacamps(:,1:rep_median_num_fix));
end

%---Post-Lesion Data---%
post_nov_pupil = NaN(length(post_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
post_rep_pupil = NaN(length(post_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
post_novel_fix_dur = NaN(length(post_files),50); %fixation durations
post_repeat_fix_dur = NaN(length(post_files),50); %fixation durations
post_novel_sac_amp = NaN(length(post_files),50); %saccade amplitudes
post_repeat_sac_amp = NaN(length(post_files),50); %saccade amplitudes
post_fix_medians = NaN(2,length(post_files));
post_sac_medians = NaN(2,length(post_files));
for file = 1:length(post_files)
    load([data_dir post_files{file}(1:8) '_' post_files{file}(end) '-fixation.mat'])
    
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
    if strcmpi('PW150416_2',post_files{file}(1:10))
        viewed(:,1) = 0;
        imgdur(:,1) = 0;
    end
    
    %Let's not even look at images in which the monkey looked away more than the
    %image was displayed for. Also should take care of EOG overflow trials
    viewed(viewed & img_dur > 12000) = 0; %1000 ms buffer for accidental lookaways and cortex timing
    imgnum(viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    viewed(:,viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    
    nov_trial_by_trial_pupil = NaN(96,1500); %each trials pupil data
    rep_trial_by_trial_pupil = NaN(96,1500); %each trials pupil data
    nov_trial_by_trial_fixdurs = NaN(96,50); %each trials fixation durations
    rep_trial_by_trial_fixdurs = NaN(96,50); %each trials fixation durations
    nov_trial_by_trial_sacamps = NaN(96,50); %each trials saccade amplitudes
    rep_trial_by_trial_sacamps = NaN(96,50); %each trials fixation amplitudes
    
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
        
        prep_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
        post_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
        prep_img_saccades = find(nov_sactimes(1,:) <= nov_img_on);
        post_img_saccades = find(nov_sactimes(2,:) > nov_img_off);
        
        nov_fix(:,post_img_fix) = [];
        nov_fix(:,prep_img_fix) = [];
        nov_fixtimes(:,post_img_fix) = [];
        nov_fixtimes(:,prep_img_fix ) = [];
        nov_sactimes(:,post_img_saccades) = [];
        nov_sactimes(:,prep_img_saccades) = [];
        
        
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
        
        prep_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        prep_img_saccades = find(rep_sactimes(1,:) <= rep_img_on);
        post_img_saccades = find(rep_sactimes(2,:) > rep_img_off);
        
        rep_fix(:,post_img_fix) = [];
        rep_fix(:,prep_img_fix) = [];
        rep_fixtimes(:,post_img_fix) = [];
        rep_fixtimes(:,prep_img_fix) = [];
        rep_sactimes(:,post_img_saccades) = [];
        rep_sactimes(:,prep_img_saccades) = [];
        
        %---Find When the monkey Blinked,when monkey looked away---%
        %pupil values at 0 diameter
        
        [nov_blink_ind,nov_nan_ind,nov_time_out(file,novrep)] = ...
            findbad_eye_ind(nov_pupil,nov_x,1000);
        [rep_blink_ind,rep_nan_ind,rep_time_out(file,novrep)] = ...
            findbad_eye_ind(rep_pupil,rep_x,1000);
        
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
    
      post_nov_pupil(file,:) = nanmean(nov_trial_by_trial_pupil);
    post_rep_pupil(file,:) = nanmean(rep_trial_by_trial_pupil);
    
    nov_median_num_fix = sum(~isnan(nov_trial_by_trial_fixdurs'));
    nov_median_num_fix(nov_median_num_fix == 0) = [];
    nov_median_num_fix = round(median(nov_median_num_fix));
    post_fix_medians(1,file) = nov_median_num_fix;
    post_novel_fix_dur(file,1:nov_median_num_fix) = ...
        nanmean(nov_trial_by_trial_fixdurs(:,1:nov_median_num_fix));
    
    rep_median_num_fix = sum(~isnan(rep_trial_by_trial_fixdurs'));
    rep_median_num_fix(rep_median_num_fix == 0) = [];
    rep_median_num_fix = round(median(rep_median_num_fix));
    post_fix_medians(2,file) = rep_median_num_fix;
    post_repeat_fix_dur(file,1:rep_median_num_fix) = ...
        nanmean(rep_trial_by_trial_fixdurs(:,1:rep_median_num_fix));
    
    
    nov_median_num_fix = sum(~isnan(nov_trial_by_trial_sacamps'));
    nov_median_num_fix(nov_median_num_fix == 0) = [];
    nov_median_num_fix = round(median(nov_median_num_fix));
    post_sac_medians(1,file) = nov_median_num_fix;
    post_novel_sac_amp(file,1:nov_median_num_fix) = ...
        nanmean(nov_trial_by_trial_sacamps(:,1:nov_median_num_fix));
    
    rep_median_num_fix = sum(~isnan(rep_trial_by_trial_sacamps'));
    rep_median_num_fix(rep_median_num_fix == 0) = [];
    rep_median_num_fix = round(median(rep_median_num_fix));
    post_sac_medians(2,file) = rep_median_num_fix;
    post_repeat_sac_amp(file,1:rep_median_num_fix) = ...
        nanmean(rep_trial_by_trial_sacamps(:,1:rep_median_num_fix));
end
%%
%---Get the median number of fixations/saccades---%
%for novel images
pre_nov_median_num_fix = median(pre_fix_medians(1,:));
pre_nov_median_num_sac = median(pre_sac_medians(1,:));
post_nov_median_num_fix = median(post_fix_medians(1,:));
post_nov_median_num_sac = median(post_sac_medians(1,:));

pre_novel_fix_dur = pre_novel_fix_dur(:,1:pre_nov_median_num_fix);
pre_novel_sac_amp = pre_novel_sac_amp(:,1:pre_nov_median_num_sac);
post_novel_fix_dur = post_novel_fix_dur(:,1:post_nov_median_num_fix);
post_novel_sac_amp = post_novel_sac_amp(:,1:post_nov_median_num_sac);

%for repeated images
pre_rep_median_num_fix = median(pre_fix_medians(2,:));
pre_rep_median_num_sac = median(pre_sac_medians(2,:));
post_rep_median_num_fix = median(post_fix_medians(2,:));
post_rep_median_num_sac = median(post_sac_medians(2,:));

pre_repeat_fix_dur = pre_repeat_fix_dur(:,1:pre_nov_median_num_fix);
pre_repeat_sac_amp = pre_repeat_sac_amp(:,1:pre_nov_median_num_sac);
post_repeat_fix_dur = post_repeat_fix_dur(:,1:post_nov_median_num_fix);
post_repeat_sac_amp = post_repeat_sac_amp(:,1:post_nov_median_num_sac);
%%

%---Absolute Plots for Eye movements---%figure
subplot(1,3,1)
hold on
errorbar(mean(pre_novel_fix_dur),std(pre_novel_fix_dur)./sqrt(size(pre_novel_fix_dur,1)),'b')
errorbar(mean(pre_repeat_fix_dur),std(pre_repeat_fix_dur)./sqrt(size(pre_repeat_fix_dur,1)),'r')
errorbar(mean(post_novel_fix_dur),std(post_novel_fix_dur)./sqrt(size(post_novel_fix_dur,1)),'g')
errorbar(mean(post_repeat_fix_dur),std(post_repeat_fix_dur)./sqrt(size(post_repeat_fix_dur,1)),'m')
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
xlim([0 max([pre_nov_median_num_fix,pre_rep_median_num_fix,post_nov_median_num_fix,post_rep_median_num_fix])+1])
axis square
title('Fixation Durations')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','SouthEast')

subplot(1,3,2)
hold on
errorbar(mean(pre_novel_sac_amp)/24,std(pre_novel_sac_amp)/24./sqrt(size(pre_novel_sac_amp,1)),'b')
errorbar(mean(pre_repeat_sac_amp)/24,std(pre_repeat_sac_amp)/24./sqrt(size(pre_repeat_sac_amp,1)),'r')
errorbar(mean(post_novel_sac_amp)/24,std(post_novel_sac_amp)/24./sqrt(size(post_novel_sac_amp,1)),'g')
errorbar(mean(post_repeat_sac_amp)/24,std(post_repeat_sac_amp)/24./sqrt(size(post_repeat_sac_amp,1)),'m')
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
xlim([0 max([pre_nov_median_num_sac,pre_rep_median_num_sac,post_nov_median_num_sac,post_rep_median_num_sac])+1])
axis square
title('Saccade Durations')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEast')

t = 1:5:7500;
subplot(1,3,3)
hold on
dofill(t,pre_nov_pupil/1000,'b',1,60)
dofill(t,pre_rep_pupil/1000,'r',1,60)
dofill(t,post_nov_pupil/1000,'g',1,60)
dofill(t,post_rep_pupil/1000,'m',1,60)
yl = ylim;
plot([500 500],[yl(1) yl(2)],'--k')
hold off
ylim(yl);
axis square
xlabel('Time from Image Onset (ms)')
ylabel('Normalize Pupil Data (a.u.)')
set(gca,'Xtick',0:500:7500)
set(gca,'XtickLabel',num2cell([0:500:7500]-500))
title('Pupil Diameter')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEast')

subtitle([pre_files{1}(1:2) ': ListRM Pre vs Post Absolute Eye Movements'])
save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_ListRM_Pre_vs_Post_Absolute_Eye_Movements'])
%%
%---Relative Plots for Eye movements---%
%normalize everything to novel condition...first fixation/saccade and all of pupil
pre_repeat_fix_dur = pre_repeat_fix_dur./mean(pre_novel_fix_dur(:,1));
pre_novel_fix_dur = pre_novel_fix_dur./mean(pre_novel_fix_dur(:,1));
pre_repeat_sac_amp = pre_repeat_sac_amp./mean(pre_novel_sac_amp(:,1));
pre_novel_sac_amp = pre_novel_sac_amp./mean(pre_novel_sac_amp(:,1));

post_repeat_fix_dur = post_repeat_fix_dur./mean(post_novel_fix_dur(:,1));
post_novel_fix_dur = post_novel_fix_dur./mean(post_novel_fix_dur(:,1));
post_repeat_sac_amp = post_repeat_sac_amp./mean(post_novel_sac_amp(:,1));
post_novel_sac_amp = post_novel_sac_amp./mean(post_novel_sac_amp(:,1));

figure
subplot(2,3,1)
hold on
errorbar(mean(pre_novel_fix_dur),std(pre_novel_fix_dur)./sqrt(size(pre_novel_fix_dur,1)),'b')
errorbar(mean(pre_repeat_fix_dur),std(pre_repeat_fix_dur)./sqrt(size(pre_repeat_fix_dur,1)),'r')
errorbar(mean(post_novel_fix_dur),std(post_novel_fix_dur)./sqrt(size(post_novel_fix_dur,1)),'g')
errorbar(mean(post_repeat_fix_dur),std(post_repeat_fix_dur)./sqrt(size(post_repeat_fix_dur,1)),'m')
hold off
xlabel('Ordinal Fixation #')
ylabel('Normalized to 1st Novel Fixation')
xlim([0 max([pre_nov_median_num_fix,pre_rep_median_num_fix,post_nov_median_num_fix,post_rep_median_num_fix])+1])
axis square
title('Fixation Durations')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','SouthEast')

subplot(2,3,2)
hold on
errorbar(mean(pre_novel_sac_amp),std(pre_novel_sac_amp)./sqrt(size(pre_novel_sac_amp,1)),'b')
errorbar(mean(pre_repeat_sac_amp),std(pre_repeat_sac_amp)./sqrt(size(pre_repeat_sac_amp,1)),'r')
errorbar(mean(post_novel_sac_amp),std(post_novel_sac_amp)./sqrt(size(post_novel_sac_amp,1)),'g')
errorbar(mean(post_repeat_sac_amp),std(post_repeat_sac_amp)./sqrt(size(post_repeat_sac_amp,1)),'m')
hold off
xlabel('Ordinal Saccade #')
ylabel('Normalized to 1st Novel Saccade')
xlim([0 max([pre_nov_median_num_sac,pre_rep_median_num_sac,post_nov_median_num_sac,post_rep_median_num_sac])+1])
axis square
title('Saccade Durations')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEast')

t = 1:5:7500;
subplot(2,3,3)
hold on
dofill(t,pre_rep_pupil./pre_nov_pupil/1000,'k',1,60)
dofill(t,post_rep_pupil./post_nov_pupil/1000,'c',1,60)
yl = ylim;
plot([0 7500],[1 1],'k')
plot([500 500],[yl(1) yl(2)],'--k')
hold off
xlim([0 7500])
ylim(yl);
axis square
xlabel('Time from Image Onset (ms)')
ylabel('Normalize Pupil Data by Novel Condition')
set(gca,'Xtick',0:500:7500)
set(gca,'XtickLabel',num2cell([0:500:7500]-500))
title('Pupil Diameter')
legend('Pre','Post','Location','NorthEast')


subplot(2,3,4)
hold on
errorbar(mean(pre_repeat_fix_dur./pre_novel_fix_dur),std(pre_repeat_fix_dur./pre_novel_fix_dur)./sqrt(size(pre_repeat_fix_dur,1)),'k')
errorbar(mean(post_repeat_fix_dur./post_novel_fix_dur),std(post_repeat_fix_dur./post_novel_fix_dur)./sqrt(size(post_repeat_fix_dur,1)),'c')
hold off
axis square
xlabel('Ordinal Fixation #')
ylabel('Normalized to Novel Fixation Durations')
legend('Pre','Post','Location','NorthEast')

subplot(2,3,5)
hold on
errorbar(mean(pre_repeat_sac_amp./pre_novel_sac_amp),std(pre_repeat_sac_amp./pre_novel_sac_amp)./sqrt(size(pre_repeat_sac_amp,1)),'k')
errorbar(mean(post_repeat_sac_amp./post_novel_sac_amp),std(post_repeat_sac_amp./post_novel_sac_amp)./sqrt(size(post_repeat_sac_amp,1)),'c')
hold off
axis square
xlabel('Ordinal Fixation #')
ylabel('Normalized to Novel Saccade Amplitudes')
legend('Pre','Post','Location','NorthEast')

subtitle([pre_files{1}(1:2) ': ListRM Pre vs Post Relative Eye Movements'])
%save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_ListRM_Pre_vs_Post_Relative_Eye_Movements'])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Spatial Coverage Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pre_coverage = NaN(2,length(pre_files));
for file = 1:length(pre_files)
       load([data_dir pre_files{file}(1:8) '_' pre_files{file}(end) '-Coverage.mat'])
       pre_coverage(1,file) = 100*nanmean(nov_coverage);
       pre_coverage(2,file) = 100*nanmean(rep_coverage);
end

post_coverage = NaN(2,length(post_files));
for file = 1:length(post_files)
       load([data_dir post_files{file}(1:8) '_' post_files{file}(end) '-Coverage.mat'])
       post_coverage(1,file) = 100*nanmean(nov_coverage);
       post_coverage(2,file) = 100*nanmean(rep_coverage);
end

figure
subplot(1,2,1)
hold on
bar([mean(pre_coverage');mean(post_coverage')]','grouped');
errorb([mean(pre_coverage');mean(post_coverage')]',...
    [std(pre_coverage')./sqrt(size(pre_coverage,2)); std(post_coverage')./sqrt(size(post_coverage,2))]')

[~,p] = ttest2(pre_coverage(1,:),post_coverage(1,:));
if p < 0.05
    plot(1,mean(pre_coverage(1,:))+1,'k*')
end
[~,p] = ttest2(pre_coverage(2,:),post_coverage(2,:));
if p < 0.05
    plot(2,mean(pre_coverage(2,:))+1,'k*')
end
hold off
legend('Pre','Post')
set(gca,'Xtick',[1,2]);
set(gca,'XtickLabel',{'Novel','Repeat'})
legend('Pre','Post')
ylabel('% Coverage')

means = [mean(pre_coverage');mean(post_coverage')];
stds = [std(pre_coverage');std(post_coverage')];

ylim([min(means(1:end))-max(stds(1:end)) max(means(1:end))+max(stds(1:end))])
title('Absolute')
axis square


%normalize to novel condition
pre_coverage(2,:) = pre_coverage(2,:)./mean(pre_coverage(1,:));
pre_coverage(1,:) = pre_coverage(1,:)./mean(pre_coverage(1,:));
post_coverage(2,:) = post_coverage(2,:)./mean(post_coverage(1,:));
post_coverage(1,:) = post_coverage(1,:)./mean(post_coverage(1,:));

subplot(1,2,2)
hold on
bar([mean(pre_coverage');mean(post_coverage')]','grouped');
errorb([mean(pre_coverage');mean(post_coverage')]',...
    [std(pre_coverage')./sqrt(size(pre_coverage,2)); std(post_coverage')./sqrt(size(post_coverage,2))]')

[~,p] = ttest2(pre_coverage(2,:),post_coverage(2,:));
if p < 0.05
    plot(2,mean(pre_coverage(2,:))+1,'k*')
end
hold off
legend('Pre','Post')
set(gca,'Xtick',[1,2]);
set(gca,'XtickLabel',{'Novel','Repeat'})
legend('Pre','Post')
ylabel('Relative Coverage')
title('Relative')
ylim([0.95 1.05])
axis square

subtitle('Viewing Coverage')
save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_ListRM_Pre_vs_Post_Image_Coverage'])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Spatial Similariy/KL Divergence Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pre_similarity = NaN(1,length(pre_files));
pre_similarity_random_chance =  NaN(1,length(pre_files));
pre_similarity_shuffled_chance =  NaN(1,length(pre_files));
pre_KL = NaN(length(pre_files),6);
pre_KL_random_chance = NaN(length(pre_files),6);
pre_KL_shuffled_chance = NaN(length(pre_files),6);
for file = 1:length(pre_files)
       load([data_dir pre_files{file}(1:8) '_' pre_files{file}(end) '-KLdivergence.mat'])
      
       pre_similarity(file) = nanmean(ObservedSimilarity); 
       pre_similarity_random_chance(file) = nanmean(ShuffledSimilarity);
       pre_similarity_shuffled_chance(file) =  nanmean(ShuffledSimilarity2);
       
       pre_KL(file,:) = nanmean(KLnorm);
       pre_KL_random_chance(file,:) = nanmean(KLshuff);
       pre_KL_shuffled_chance(file,:) = nanmean(KLshuff2);
end

post_similarity = NaN(1,length(post_files));
post_similarity_random_chance =  NaN(1,length(post_files));
post_similarity_shuffled_chance =  NaN(1,length(post_files));
post_KL = NaN(length(post_files),6);
post_KL_random_chance = NaN(length(post_files),6);
post_KL_shuffled_chance = NaN(length(post_files),6);
for file = 1:length(post_files)
       load([data_dir post_files{file}(1:8) '_' post_files{file}(end) '-KLdivergence.mat'])
      
       post_similarity(file) = nanmean(ObservedSimilarity); 
       post_similarity_random_chance(file) = nanmean(ShuffledSimilarity);
       post_similarity_shuffled_chance(file) =  nanmean(ShuffledSimilarity2);
       
       post_KL(file,:) = nanmean(KLnorm);
       post_KL_random_chance(file,:) = nanmean(KLshuff);
       post_KL_shuffled_chance(file,:) = nanmean(KLshuff2);
end


figure
subplot(4,4,[1 2 5 6])
hold on
errorbar(mean(pre_KL(:,1:5)),std(pre_KL(:,1:5))./sqrt(size(pre_KL,1)),'b')
errorbar(mean(post_KL(:,1:5)),std(post_KL(:,1:5))./sqrt(size(post_KL,1)),'r')
plot(mean(pre_KL_random_chance(:,1:5)),'k')
plot(mean(post_KL_random_chance(:,1:5)),'--k')
plot(mean(pre_KL_shuffled_chance(:,1:5)),'g')
plot(mean(post_KL_shuffled_chance(:,1:5)),'m')
pvals = [];
for g = 1:5
   [~,p] = ttest2(pre_KL(:,g),post_KL(:,g));
   pvals(g) = p;
   if p < 0.05
       plot(g,mean(pre_KL(:,g))+std(pre_KL(:,g)),'*k')
   end
end
hold off
legend('Observed Pre','Observed Post','Pre Chance','Post Chance','Pre Shuffled','Post Shuffled','Location','SouthEast')
set(gca,'Xtick',[1:5])
set(gca,'XtickLabel',{'1-5','6-10','11-15','15-20','21-25'})
% xlabel('Fixation #s')
ylabel('KL Divergence (bits)')
title('Absolute')

subplot(4,4,[3 7])
hold on
bar(1,mean(pre_KL(:,6)))
errorbar(1,mean(pre_KL(:,6)),std(pre_KL(:,6))./sqrt(size(pre_KL,1)),'k','linewidth',3)
bar(2,mean(post_KL(:,6)),'r')
errorbar(2,mean(post_KL(:,6)),std(post_KL(:,6))./sqrt(size(post_KL,1)),'k','linewidth',3)
errorbar(1,mean(pre_KL_shuffled_chance(:,6)),std(pre_KL_shuffled_chance(:,6))./sqrt(size(pre_KL_shuffled_chance,1)),'g','linewidth',3)
errorbar(2,mean(post_KL_shuffled_chance(:,6)),std(post_KL_shuffled_chance(:,6))./sqrt(size(post_KL_shuffled_chance,1)),'m','linewidth',3)
[~,p] = ttest2(pre_KL(:,6),post_KL(:,6));
if p < 0.05
   plot(1.5,mean(pre_KL(:,6)+1),'k*')
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Pre','Post'})
ylabel('KL Divergence (bits)')
title('Absolute: All Fixations')

pre_KL = pre_KL./pre_KL_shuffled_chance;
post_KL = post_KL./post_KL_shuffled_chance;

subplot(4,4,[9 10 13 14])
hold on
errorbar(mean(pre_KL(:,1:5)),std(pre_KL(:,1:5))./sqrt(size(pre_KL,1)),'b')
errorbar(mean(post_KL(:,1:5)),std(post_KL(:,1:5))./sqrt(size(post_KL,1)),'r')
pvals = [];
for g = 1:5
   [~,p] = ttest2(pre_KL(:,g),post_KL(:,g));
   pvals(g) =  p;
   if p < 0.05
       plot(g,mean(pre_KL(:,g))+0.05,'k*')
   end
end
hold off
legend('Pre','Post','Location','SouthEast')
set(gca,'Xtick',[1:5])
set(gca,'XtickLabel',{'1-5','6-10','11-15','15-20','21-25'})
xlabel('Fixation #s')
ylabel('KL Divergence Relative to Shuffled Measures')
title('Relative')

subplot(4,4,[11 15])
hold on
bar(1,mean(pre_KL(:,6)))
errorbar(1,mean(pre_KL(:,6)),std(pre_KL(:,6))./sqrt(size(pre_KL,1)),'k','linewidth',3)
bar(2,mean(post_KL(:,6)),'r')
errorbar(2,mean(post_KL(:,6)),std(post_KL(:,6))./sqrt(size(post_KL,1)),'k','linewidth',3)
[~,p] = ttest2(pre_KL(:,6),post_KL(:,6));
if p < 0.05
   plot(1.5,mean(pre_KL(:,6)+0.05),'k*')
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Pre','Post'})
ylabel('KL Divergence Relative to Shuffled Measures')
title('Relative: All Fixations')

subplot(4,4,[4 8])
hold on
bar(1,mean(pre_similarity),'b')
errorbar(1,mean(pre_similarity),std(pre_similarity)./sqrt(size(pre_similarity,2)),'k','linewidth',3)
bar(2,mean(post_similarity),'r')
errorbar(2,mean(post_similarity),std(post_similarity)./sqrt(size(post_similarity,2)),'k','linewidth',3)
errorbar(1,mean(pre_similarity_shuffled_chance),std(pre_similarity_shuffled_chance)./sqrt(size(pre_similarity_shuffled_chance,2)),'g','linewidth',3)
errorbar(2,mean(post_similarity_shuffled_chance),std(post_similarity_shuffled_chance)./sqrt(size(post_similarity_shuffled_chance,2)),'m','linewidth',3)
[~,p] = ttest2(pre_similarity,post_similarity);
if p < 0.05
   plot(1.5,mean(pre_similarity)+100,'*k')
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Pre','Post'})
ylabel('Coverage Similarity')
title('Absolute')

pre_similarity = pre_similarity./pre_similarity_shuffled_chance;
post_similarity = post_similarity./post_similarity_shuffled_chance;

subplot(4,4,[12 16])
hold on
bar(1,mean(pre_similarity),'b')
errorbar(1,mean(pre_similarity),std(pre_similarity)./sqrt(size(pre_similarity,2)),'k','linewidth',3)
bar(2,mean(post_similarity),'r')
errorbar(2,mean(post_similarity),std(post_similarity)./sqrt(size(post_similarity,2)),'k','linewidth',3)
[~,p] = ttest2(pre_similarity,post_similarity);
if p < 0.05
   plot(1.5,mean(pre_similarity)+0.1,'*k')
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Pre','Post'})
ylabel('Relative Coverage Similarity')
title('Relative')

subtitle([pre_files{1}(1:2) ': Fixation Location Similarity'])
save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_ListRM_Pre_vs_Post_SpatialSimilarity'])
