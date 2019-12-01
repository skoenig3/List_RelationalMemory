%written by Seth Konig 6/8/16
clar
fixwin = 3.5; %size of the fixation window/2. Was a width of 7 dva
imageX = 800;
imageY = 600;

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Cortex Data\';
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Image Sets\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Figures\';

%---Vivian---%
% pre_files = {'PW150416.4','PW150417.2','PW150420.2','PW150421.2',...
%     'PW150422.2','PW150423.2','PW150424.2','PW150427.2',...
%     'PW150428.2','PW150429.2','PW150430.2','PW150501.2',...
%     'PW150504.2','PW150505.2','PW150506.2','PW150507.2',...
%     'PW150511.2'};
% 
% post_files = {'PW160310.2','PW160311.2','PW160314.2','PW160315.2','PW160316.2',...
%     'PW160317.2','PW160318.2','PW160321.2','PW160322.2','PW160323.2',...
%     'PW160324.2','PW160325.2','PW160325.2','PW160329.2','PW160330.2'};
%
% % %---Red---%
% pre_files = {'RR150423.2','RR150424.2','RR150427.2','RR150428.2','RR150429.2',...
%     'RR150430.2','RR150501.2','RR150504.2','RR150505.2','RR150506.2',...
%     'RR150507.2','RR150508.2','RR150511.2','RR150512.2','RR150513.2',...
%     'RR150515.2','RR150518.2','RR150519.2','RR150520.2','RR150521.2'};
%
% post_files = {'RR160324.2','RR160325.1','RR160328.2','RR160330.2',...
%                 'RR160331.2','RR160401.2','RR160405.2','RR160406.2',...
%                 'RR160407.2','RR160408.2','RR160411.2','RR160412.2',...
%                 'RR160413.2','RR160414.2','RR160415.2','RR160418.2',...
%                 'RR160419.2','RR160420.2','RR160421.2','RR160422.2'};
%
%---Tobii---%
% pre_files = {'TO150513.2','TO150514.2','TO150515.2','TO150518.2','TO150519.2',...
%     'TO150520.2','TO150521.2','TO150522.2','TO150526.2','TO150527.2',...
%     'TO150528.2','TO150529.2','TO150601.2','TO150602.2','TO150603.2',...
%     'TO150604.2','TO150605.2'};
%
% post_files = {'TO170404.2','TO170405.2','TO170406.2','TO170407.2',...
%     'TO170410.2','TO170411.2','TO170412.2','TO170414.2',...
%     'TO170419.2','TO170421.2','TO170425.2','TO170426.2',...
%     'TO170427.2','TO170428.2','TO170501.2'};

%---Manfred---%
% pre_files = {'MF170111.2','MF170112.2','MF170117.2','MF170118.2',...
%                 'MF170119.2','MF170120.2','MF170123.2','MF170126.2',...
%                 'MF170130.2','MF170131.2','MF170201.2','MF170202.2',...
%                 'MF170209.2','MF170210.2','MF170213.2'};
% 
% post_files = {'MF180425.2','MF180426.2','MF180427.2','MF180430.2',...
%                 'MF180501.2','MF180502.2','MF180503.2','MF180504.2',...
%                 'MF180507.2','MF180508.2','MF180509.2','MF180511.2',...
%                 'MF180514.2','MF180515.2','MF180516.2','MF180517.2',...
%                 'MF180518.2','MF180521.2'};


%---Combined---%
pre_files = {'PW150416.4','PW150417.2','PW150420.2','PW150421.2',...
    'PW150422.2','PW150423.2','PW150424.2','PW150427.2',...
    'PW150428.2','PW150429.2','PW150430.2','PW150501.2',...
    'PW150504.2','PW150505.2','PW150506.2','PW150507.2',...
    'PW150511.2',...
    'TO150513.2','TO150514.2','TO150515.2','TO150518.2','TO150519.2',...
    'TO150520.2','TO150521.2','TO150522.2','TO150526.2','TO150527.2',...
    'TO150528.2','TO150529.2','TO150601.2','TO150602.2','TO150603.2',...
    'TO150604.2','TO150605.2',...
    'RR150423.2','RR150424.2','RR150427.2','RR150428.2','RR150429.2',...
    'RR150430.2','RR150501.2','RR150504.2','RR150505.2','RR150506.2',...
    'RR150507.2','RR150508.2','RR150511.2','RR150512.2','RR150513.2',...
    'RR150515.2','RR150518.2','RR150519.2','RR150520.2','RR150521.2',...
    'MF170111.2','MF170112.2','MF170117.2','MF170118.2',...
    'MF170119.2','MF170120.2','MF170123.2','MF170126.2',...
    'MF170130.2','MF170131.2','MF170201.2','MF170202.2',...
    'MF170209.2','MF170210.2','MF170213.2'};

post_files = {'PW160310.2','PW160311.2','PW160314.2','PW160315.2','PW160316.2',...
    'PW160317.2','PW160318.2','PW160321.2','PW160322.2','PW160323.2',...
    'PW160324.2','PW160325.2','PW160325.2','PW160329.2','PW160330.2',...
    'TO170404.2','TO170405.2','TO170406.2','TO170407.2',...
    'TO170410.2','TO170411.2','TO170412.2','TO170414.2',...
    'TO170419.2','TO170421.2','TO170425.2','TO170426.2',...
    'TO170427.2','TO170428.2','TO170501.2',...
    'RR160324.2','RR160325.1','RR160328.2','RR160330.2',...
    'RR160331.2','RR160401.2','RR160405.2','RR160406.2',...
    'RR160407.2','RR160408.2','RR160411.2','RR160412.2',...
    'RR160413.2','RR160414.2','RR160415.2','RR160418.2',...
    'RR160419.2','RR160420.2','RR160421.2','RR160422.2',...
    'MF180425.2','MF180426.2','MF180427.2','MF180430.2',...
    'MF180501.2','MF180502.2','MF180503.2','MF180504.2',...
    'MF180507.2','MF180508.2','MF180509.2','MF180511.2',...
    'MF180514.2','MF180515.2','MF180516.2','MF180517.2',...
    'MF180518.2','MF180521.2'};


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
        
        if nov_img_off-nov_img_on > 10500
            continue
        end
        
        nov_x = fixationstats{viewed(1,novrep)}.XY(1,nov_img_on:nov_img_off);
        nov_y = fixationstats{viewed(1,novrep)}.XY(2,nov_img_on:nov_img_off);
        nov_pupil = pupildata{viewed(1,novrep)}(round(nov_img_on/5):round(nov_img_off/5));
        
        nov_fix = fixationstats{viewed(1,novrep)}.fixations;
        nov_fixtimes = fixationstats{viewed(1,novrep)}.fixationtimes;
        nov_sactimes = fixationstats{viewed(1,novrep)}.saccadetimes;
        
        prep_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
        postp_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
        prep_img_saccades = find(nov_sactimes(1,:) <= nov_img_on);
        postp_img_saccades = find(nov_sactimes(2,:) > nov_img_off);
        
        nov_fix(:,postp_img_fix) = [];
        nov_fix(:,prep_img_fix) = [];
        nov_fixtimes(:,postp_img_fix) = [];
        nov_fixtimes(:,prep_img_fix ) = [];
        nov_sactimes(:,postp_img_saccades) = [];
        nov_sactimes(:,prep_img_saccades) = [];
        
        
        %for repeat images
        rep_allval = per(viewed(2,novrep)).allval;
        rep_alltim = per(viewed(2,novrep)).alltim;
        rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        
        if rep_img_off-rep_img_on > 10050
            continue
        end
        
        rep_x = fixationstats{viewed(2,novrep)}.XY(1,rep_img_on:rep_img_off);
        rep_y = fixationstats{viewed(2,novrep)}.XY(2,rep_img_on:rep_img_off);
        rep_pupil = pupildata{viewed(2,novrep)}(round(rep_img_on/5):round(rep_img_off/5));
        
        rep_fix = fixationstats{viewed(2,novrep)}.fixations;
        rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
        rep_sactimes = fixationstats{viewed(2,novrep)}.saccadetimes;
        
        prep_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        postp_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        prep_img_saccades = find(rep_sactimes(1,:) <= rep_img_on);
        postp_img_saccades = find(rep_sactimes(2,:) > rep_img_off);
        
        rep_fix(:,postp_img_fix) = [];
        rep_fix(:,prep_img_fix) = [];
        rep_fixtimes(:,postp_img_fix) = [];
        rep_fixtimes(:,prep_img_fix) = [];
        rep_sactimes(:,postp_img_saccades) = [];
        rep_sactimes(:,prep_img_saccades) = [];
        
        %---Analyze pupil data---%
        %find when monkey blinked
        [nov_blink_ind,nov_nan_ind,nov_time_out] = findbad_eye_ind(nov_pupil,nov_x,1000);
        [rep_blink_ind,rep_nan_ind,rep_time_out] = findbad_eye_ind(rep_pupil,rep_x,1000);
        
        nov_pupil = pupildata{viewed(1,novrep)}((round(nov_img_on/5)-100):round(nov_img_off/5));
        rep_pupil = pupildata{viewed(2,novrep)}((round(rep_img_on/5)-100):round(rep_img_off/5));
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
        
        %normalize to 100 ms of prestimulus levels
        nov_pupil = nov_pupil+2000;
        rep_pupil = rep_pupil+2000;
        nov_pupil = nov_pupil./abs(mean(nov_pupil(81:100)));
        rep_pupil = rep_pupil./abs(mean(rep_pupil(81:100)));
        
        %only look at the 1st 5 seconds the images is up
        if length(nov_pupil) >= 1500
            nov_pupil = nov_pupil(1:1500);
        else
            nov_pupil = [nov_pupil; NaN(1500-length(nov_pupil),1)];
        end
        if length(rep_pupil) >= 1500
            rep_pupil = rep_pupil(1:1500);
        else
            rep_pupil = [rep_pupil; NaN(1500-length(rep_pupil),1)];
        end
        
        nov_trial_by_trial_pupil(novrep,:) = nov_pupil;
        rep_trial_by_trial_pupil(novrep,:) = rep_pupil;
        
        %---Calculate Fixation Durations and Saccade Amplitudes---%
        %fixation duration by ordinal fixation number
        nov_trial_by_trial_fixdurs(novrep,1:length(nov_fixtimes)) = nov_fixtimes(2,:)-nov_fixtimes(1,:)+1;
        rep_trial_by_trial_fixdurs(novrep,1:length(rep_fixtimes)) = rep_fixtimes(2,:)-rep_fixtimes(1,:)+1;
        
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
    
    nov_median_num_fix = sum(~isnan(nov_trial_by_trial_fixdurs'));
    nov_median_num_fix(nov_median_num_fix == 0) = [];
    nov_median_num_fix = round(median(nov_median_num_fix));
    pre_novel_fix_dur(file,1:nov_median_num_fix) = ...
        nanmean(nov_trial_by_trial_fixdurs(:,1:nov_median_num_fix));
    
    rep_median_num_fix = sum(~isnan(rep_trial_by_trial_fixdurs'));
    rep_median_num_fix(rep_median_num_fix == 0) = [];
    rep_median_num_fix = round(median(rep_median_num_fix));
    pre_repeat_fix_dur(file,1:rep_median_num_fix) = ...
        nanmean(rep_trial_by_trial_fixdurs(:,1:rep_median_num_fix));
    
    nov_median_num_sac = sum(~isnan(nov_trial_by_trial_sacamps'));
    nov_median_num_sac(nov_median_num_sac == 0) = [];
    nov_median_num_sac = round(median(nov_median_num_sac));
    pre_novel_sac_amp(file,1:nov_median_num_sac) = ...
        nanmean(nov_trial_by_trial_sacamps(:,1:nov_median_num_sac));
    
    rep_median_num_sac = sum(~isnan(rep_trial_by_trial_sacamps'));
    rep_median_num_sac(rep_median_num_sac == 0) = [];
    rep_median_num_sac = round(median(rep_median_num_sac));
    pre_repeat_sac_amp(file,1:rep_median_num_sac) = ...
        nanmean(rep_trial_by_trial_sacamps(:,1:rep_median_num_sac));
    
    pre_nov_pupil(file,:) = nanmean(nov_trial_by_trial_pupil);
    pre_rep_pupil(file,:) = nanmean(rep_trial_by_trial_pupil);
    
end

%---Post-Lesion Data---%
post_nov_pupil = NaN(length(post_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
post_rep_pupil = NaN(length(post_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
post_novel_fix_dur = NaN(length(post_files),50); %fixation durations
post_repeat_fix_dur = NaN(length(post_files),50); %fixation durations
post_novel_sac_amp = NaN(length(post_files),50); %saccade amplitudes
post_repeat_sac_amp = NaN(length(post_files),50); %saccade amplitudes
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
        
        if nov_img_off-nov_img_on > 10500
            continue
        end
        
        nov_x = fixationstats{viewed(1,novrep)}.XY(1,nov_img_on:nov_img_off);
        nov_y = fixationstats{viewed(1,novrep)}.XY(2,nov_img_on:nov_img_off);
        nov_pupil = pupildata{viewed(1,novrep)}(round(nov_img_on/5):round(nov_img_off/5));
        
        nov_fix = fixationstats{viewed(1,novrep)}.fixations;
        nov_fixtimes = fixationstats{viewed(1,novrep)}.fixationtimes;
        nov_sactimes = fixationstats{viewed(1,novrep)}.saccadetimes;
        
        prep_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
        postp_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
        prep_img_saccades = find(nov_sactimes(1,:) <= nov_img_on);
        postp_img_saccades = find(nov_sactimes(2,:) > nov_img_off);
        
        nov_fix(:,postp_img_fix) = [];
        nov_fix(:,prep_img_fix) = [];
        nov_fixtimes(:,postp_img_fix) = [];
        nov_fixtimes(:,prep_img_fix ) = [];
        nov_sactimes(:,postp_img_saccades) = [];
        nov_sactimes(:,prep_img_saccades) = [];
        
        
        %for repeat images
        rep_allval = per(viewed(2,novrep)).allval;
        rep_alltim = per(viewed(2,novrep)).alltim;
        rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        
        if rep_img_off-rep_img_on > 10050
            continue
        end
        
        rep_x = fixationstats{viewed(2,novrep)}.XY(1,rep_img_on:rep_img_off);
        rep_y = fixationstats{viewed(2,novrep)}.XY(2,rep_img_on:rep_img_off);
        rep_pupil = pupildata{viewed(2,novrep)}(round(rep_img_on/5):round(rep_img_off/5));
        
        rep_fix = fixationstats{viewed(2,novrep)}.fixations;
        rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
        rep_sactimes = fixationstats{viewed(2,novrep)}.saccadetimes;
        
        prep_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        postp_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        prep_img_saccades = find(rep_sactimes(1,:) <= rep_img_on);
        postp_img_saccades = find(rep_sactimes(2,:) > rep_img_off);
        
        rep_fix(:,postp_img_fix) = [];
        rep_fix(:,prep_img_fix) = [];
        rep_fixtimes(:,postp_img_fix) = [];
        rep_fixtimes(:,prep_img_fix) = [];
        rep_sactimes(:,postp_img_saccades) = [];
        rep_sactimes(:,prep_img_saccades) = [];
        
        %---Analyze pupil data---%
        %find when monkey blinked
        [nov_blink_ind,nov_nan_ind,nov_time_out] = findbad_eye_ind(nov_pupil,nov_x,1000);
        [rep_blink_ind,rep_nan_ind,rep_time_out] = findbad_eye_ind(rep_pupil,rep_x,1000);
        
        nov_pupil = pupildata{viewed(1,novrep)}((round(nov_img_on/5)-100):round(nov_img_off/5));
        rep_pupil = pupildata{viewed(2,novrep)}((round(rep_img_on/5)-100):round(rep_img_off/5));
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
        
        %normalize to 100 ms of prestimulus levels
        nov_pupil = nov_pupil+2000;
        rep_pupil = rep_pupil+2000;
        nov_pupil = nov_pupil./abs(mean(nov_pupil(81:100)));
        rep_pupil = rep_pupil./abs(mean(rep_pupil(81:100)));
        
        %only look at the 1st 5 seconds the images is up
        if length(nov_pupil) >= 1500
            nov_pupil = nov_pupil(1:1500);
        else
            nov_pupil = [nov_pupil; NaN(1500-length(nov_pupil),1)];
        end
        if length(rep_pupil) >= 1500
            rep_pupil = rep_pupil(1:1500);
        else
            rep_pupil = [rep_pupil; NaN(1500-length(rep_pupil),1)];
        end
        
        nov_trial_by_trial_pupil(novrep,:) = nov_pupil;
        rep_trial_by_trial_pupil(novrep,:) = rep_pupil;
        
        %---Calculate Fixation Durations and Saccade Amplitudes---%
        %fixation duration by ordinal fixation number
        nov_trial_by_trial_fixdurs(novrep,1:length(nov_fixtimes)) = nov_fixtimes(2,:)-nov_fixtimes(1,:)+1;
        rep_trial_by_trial_fixdurs(novrep,1:length(rep_fixtimes)) = rep_fixtimes(2,:)-rep_fixtimes(1,:)+1;
        
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
    
    nov_median_num_fix = sum(~isnan(nov_trial_by_trial_fixdurs'));
    nov_median_num_fix(nov_median_num_fix == 0) = [];
    nov_median_num_fix = round(median(nov_median_num_fix));
    post_novel_fix_dur(file,1:nov_median_num_fix) = ...
        nanmean(nov_trial_by_trial_fixdurs(:,1:nov_median_num_fix));
    
    rep_median_num_fix = sum(~isnan(rep_trial_by_trial_fixdurs'));
    rep_median_num_fix(rep_median_num_fix == 0) = [];
    rep_median_num_fix = round(median(rep_median_num_fix));
    post_repeat_fix_dur(file,1:rep_median_num_fix) = ...
        nanmean(rep_trial_by_trial_fixdurs(:,1:rep_median_num_fix));
    
    nov_median_num_sac = sum(~isnan(nov_trial_by_trial_sacamps'));
    nov_median_num_sac(nov_median_num_sac == 0) = [];
    nov_median_num_sac = round(median(nov_median_num_sac));
    post_novel_sac_amp(file,1:nov_median_num_sac) = ...
        nanmean(nov_trial_by_trial_sacamps(:,1:nov_median_num_sac));
    
    rep_median_num_sac = sum(~isnan(rep_trial_by_trial_sacamps'));
    rep_median_num_sac(rep_median_num_sac == 0) = [];
    rep_median_num_sac = round(median(rep_median_num_sac));
    post_repeat_sac_amp(file,1:rep_median_num_sac) = ...
        nanmean(rep_trial_by_trial_sacamps(:,1:rep_median_num_sac));
    
    post_nov_pupil(file,:) = nanmean(nov_trial_by_trial_pupil);
    post_rep_pupil(file,:) = nanmean(rep_trial_by_trial_pupil);
end


%%

npre = length(pre_files);
npost = length(post_files);
%%
median_pre_nov_fix_count = median(sum(~isnan(pre_novel_fix_dur')));
median_pre_rep_fix_count = median(sum(~isnan(pre_repeat_fix_dur')));
median_post_nov_fix_count = median(sum(~isnan(post_novel_fix_dur')));
median_post_rep_fix_count = median(sum(~isnan(post_repeat_fix_dur')));
%%

pre_nov = mean(pre_novel_fix_dur(:,3:15),2);
pre_rep = mean(pre_repeat_fix_dur(:,3:15),2);
post_nov = mean(post_novel_fix_dur(:,3:15),2);
post_rep = mean(post_repeat_fix_dur(:,3:15),2);

ids = [[ones(npre,1); 2*ones(npre,1); ones(npost,1); 2*ones(npost,1)] ...
    [ones(npre,1); ones(npre,1); 2*ones(npost,1); 2*ones(npost,1)]];

vals = [pre_nov; pre_rep; post_nov; post_rep];
[P_ANOVA_fixdurs] = anovan(vals,ids,'model','interaction','varnames',{'Novelty','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

figure

subplot(2,2,1)
hold all
errorbar(nanmean(pre_novel_fix_dur(:,1:median_pre_nov_fix_count)),nanstd(pre_novel_fix_dur(:,1:median_pre_nov_fix_count))./sqrt(npre))
errorbar(nanmean(pre_repeat_fix_dur(:,1:median_pre_rep_fix_count)),nanstd(pre_repeat_fix_dur(:,1:median_pre_rep_fix_count))./sqrt(npre))
errorbar(nanmean(post_novel_fix_dur(:,1:median_post_nov_fix_count)),nanstd(post_novel_fix_dur(:,1:median_post_nov_fix_count))./sqrt(npost))
errorbar(nanmean(post_repeat_fix_dur(:,1:median_post_rep_fix_count)),nanstd(post_repeat_fix_dur(:,1:median_post_rep_fix_count))./sqrt(npost))
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEastOutside')
xlim([0 median_post_nov_fix_count+1])
title('Raw Fixation Durations')

[~,p_nov] = ttest2(pre_nov,post_nov);
[~,p_rep] = ttest2(pre_rep,post_rep);

subplot(2,2,3)
errorb([[mean(pre_nov);mean(post_nov)] [mean(pre_rep);mean(post_rep)]]',...
    [[std(pre_nov)./sqrt(npre) ;std(post_nov)./sqrt(npost)] ...
    [std(pre_rep)./sqrt(npre) ;std(post_rep)./sqrt(npost)]]')
hold on
if p_nov < 0.05
    plot(1,mean(pre_nov)+2*std(pre_nov),'k*')
end
if p_rep < 0.05
    plot(2,mean(pre_rep)+2*std(pre_rep),'k*')
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Fixation Duration (ms)')
title(sprintf(['ANOVA Ordinal Fixations 3-15: \n p_{lesion} = ' num2str(P_ANOVA_fixdurs(2),2) ...
    ' , p_{Novelty} = ' num2str(P_ANOVA_fixdurs(1),2) ', p_{inter} = ' num2str(P_ANOVA_fixdurs(3),2)]))
axis square
box off
legend('Pre','Post','Location','NorthEastOutside')


norm_pre_nov_fixdur = pre_novel_fix_dur(:,1:median_pre_nov_fix_count);
norm_pre_rep_fixdur = pre_repeat_fix_dur(:,1:median_pre_rep_fix_count);
norm_pre_rep_fixdur = norm_pre_rep_fixdur./nanmean(norm_pre_nov_fixdur(:,1));
norm_pre_nov_fixdur = norm_pre_nov_fixdur./nanmean(norm_pre_nov_fixdur(:,1));
norm_post_rep_fixdur = post_repeat_fix_dur(:,1:median_post_rep_fix_count);
norm_post_nov_fixdur = post_novel_fix_dur(:,1:median_post_nov_fix_count);
norm_post_rep_fixdur = norm_post_rep_fixdur./nanmean(norm_post_nov_fixdur(:,1));
norm_post_nov_fixdur = norm_post_nov_fixdur./nanmean(norm_post_nov_fixdur(:,1));


norm_pre_nov = mean(norm_pre_nov_fixdur(:,3:15),2);
norm_pre_rep = mean(norm_pre_rep_fixdur(:,3:15),2);
norm_post_nov = mean(norm_post_nov_fixdur(:,3:15),2);
norm_post_rep = mean(norm_post_rep_fixdur(:,3:15),2);

vals = [norm_pre_nov; norm_pre_rep; norm_post_nov; norm_post_rep];
[P_ANOVA_fixdurs] = anovan(vals,ids,'model','interaction','varnames',{'Novelty','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

subplot(2,2,2)
hold all
errorbar(nanmean(norm_pre_nov_fixdur),nanstd(norm_pre_nov_fixdur)./sqrt(npre))
errorbar(nanmean(norm_pre_rep_fixdur),nanstd(norm_pre_rep_fixdur)./sqrt(npre))
errorbar(nanmean(norm_post_nov_fixdur),nanstd(norm_post_nov_fixdur)./sqrt(npost))
errorbar(nanmean(norm_post_rep_fixdur),nanstd(norm_post_rep_fixdur)./sqrt(npost))
hold off
xlabel('Ordinal Fixation #')
ylabel('Relative Fixation Duration ')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEastOutside')
xlim([0 median_post_nov_fix_count+1])
title('Normalized Fixation Durations')


[~,norm_p_nov] = ttest2(norm_pre_nov,norm_post_nov);
[~,norm_p_rep] = ttest2(norm_pre_rep,norm_post_rep);

subplot(2,2,4)
errorb([[mean(norm_pre_nov);mean(norm_post_nov)] [mean(norm_pre_rep);mean(norm_post_rep)]]',...
    [[std(norm_pre_nov)./sqrt(npre) ;std(norm_post_nov)./sqrt(npost)] ...
    [std(norm_pre_rep)./sqrt(npre) ;std(norm_post_rep)./sqrt(npost)]]')
hold on
if norm_p_nov < 0.05
    plot(1,mean(norm_pre_nov)+2*std(norm_pre_nov),'k*')
end
if norm_p_rep < 0.05
    plot(2,mean(norm_pre_rep)+2*std(norm_pre_rep),'k*')
end
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Relative Fixation Duration ')
title(sprintf(['ANOVA Ordinal Fixations 3-15: \n p_{lesion} = ' num2str(P_ANOVA_fixdurs(2),2) ...
    ' , p_{Novelty} = ' num2str(P_ANOVA_fixdurs(1),2) ', p_{inter} = ' num2str(P_ANOVA_fixdurs(3),2)]))
axis square
box off
legend('Pre','Post','Location','NorthEastOutside')
%%subtitle([pre_files{1}(1:2) ': Pre vs. Post ListRM'])

%%
median_pre_nov_sac_count = median(sum(~isnan(pre_novel_sac_amp')));
median_pre_rep_sac_count = median(sum(~isnan(pre_repeat_sac_amp')));
median_post_nov_sac_count = median(sum(~isnan(post_novel_sac_amp')));
median_post_rep_sac_count = median(sum(~isnan(post_repeat_sac_amp')));


npre = size(pre_novel_sac_amp,1);
npost = size(post_novel_sac_amp,1);

pre_nov = nanmean(pre_novel_sac_amp(:,2:5),2)/24;
pre_rep = nanmean(pre_repeat_sac_amp(:,2:5),2)/24;
post_nov = nanmean(post_novel_sac_amp(:,2:5),2)/24;
post_rep = nanmean(post_repeat_sac_amp(:,2:5),2)/24;

ids = [[ones(npre,1); 2*ones(npre,1); ones(npost,1); 2*ones(npost,1)] ...
    [ones(npre,1); ones(npre,1); 2*ones(npost,1); 2*ones(npost,1)]];

vals = [pre_nov; pre_rep; post_nov; post_rep];
[P_ANOVA_sacamps] = anovan(vals,ids,'model','interaction','varnames',{'Novelty','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

figure

subplot(2,2,1)
hold all
errorbar(nanmean(pre_novel_sac_amp(:,1:median_pre_nov_sac_count))/24,nanstd(pre_novel_sac_amp(:,1:median_pre_nov_sac_count))/24./sqrt(npre))
errorbar(nanmean(pre_repeat_sac_amp(:,1:median_pre_rep_sac_count))/24,nanstd(pre_repeat_sac_amp(:,1:median_pre_rep_sac_count))/24./sqrt(npre))
errorbar(nanmean(post_novel_sac_amp(:,1:median_post_nov_sac_count))/24,nanstd(post_novel_sac_amp(:,1:median_post_nov_sac_count))/24./sqrt(npost))
errorbar(nanmean(post_repeat_sac_amp(:,1:median_post_rep_sac_count))/24,nanstd(post_repeat_sac_amp(:,1:median_post_rep_sac_count))/24./sqrt(npost))
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEastOutside')
xlim([0 median_post_nov_sac_count+1])
title('Raw Saccade Amplitudes')

[~,p_nov] = ttest2(pre_nov,post_nov);
[~,p_rep] = ttest2(pre_rep,post_rep);

subplot(2,2,3)
errorb([[mean(pre_nov);mean(post_nov)] [mean(pre_rep);mean(post_rep)]]',...
    [[std(pre_nov)./sqrt(npre) ;std(post_nov)./sqrt(npost)] ...
    [std(pre_rep)./sqrt(npre) ;std(post_rep)./sqrt(npost)]]')
hold on
if p_nov < 0.05 
    plot(1,mean(pre_nov)+2*std(pre_nov),'k*')
end
if p_rep < 0.05
    plot(2,mean(pre_rep)+2*std(pre_rep),'k*')
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Saccade Amplitude (dva)')
title(sprintf(['ANOVA Saccades 2-5: \n p_{lesion} = ' num2str(P_ANOVA_sacamps(2),2) ...
    ' , p_{Novelty} = ' num2str(P_ANOVA_sacamps(1),2) ', p_{inter} = ' num2str(P_ANOVA_sacamps(3),2)]))
axis square
box off
legend('Pre','Post','Location','NorthEastOutside')


norm_pre_nov_sacamp = pre_novel_sac_amp(:,1:median_pre_nov_sac_count);
norm_pre_rep_sacamp = pre_repeat_sac_amp(:,1:median_pre_rep_sac_count);
norm_pre_rep_sacamp = norm_pre_rep_sacamp./nanmean(norm_pre_rep_sacamp(:,1));
norm_pre_nov_sacamp = norm_pre_nov_sacamp./nanmean(norm_pre_nov_sacamp(:,1));
norm_post_nov_sacamp = post_novel_sac_amp(:,1:median_post_nov_sac_count);
norm_post_rep_sacamp = post_repeat_sac_amp(:,1:median_post_rep_sac_count);
norm_post_rep_sacamp = norm_post_rep_sacamp./nanmean(norm_post_rep_sacamp(:,1));
norm_post_nov_sacamp = norm_post_nov_sacamp./nanmean(norm_post_nov_sacamp(:,1));


norm_pre_nov = mean(norm_pre_nov_sacamp(:,3:15),2);
norm_pre_rep = mean(norm_pre_rep_sacamp(:,3:15),2);
norm_post_nov = mean(norm_post_nov_sacamp(:,3:15),2);
norm_post_rep = mean(norm_post_rep_sacamp(:,3:15),2);

vals = [norm_pre_nov; norm_pre_rep; norm_post_nov; norm_post_rep];
[P_ANOVA_sacamps] = anovan(vals,ids,'model','interaction','varnames',{'Novelty','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

subplot(2,2,2)
hold all
errorbar(nanmean(norm_pre_nov_sacamp),nanstd(norm_pre_nov_sacamp)./sqrt(npre))
errorbar(nanmean(norm_pre_rep_sacamp),nanstd(norm_pre_rep_sacamp)./sqrt(npre))
errorbar(nanmean(norm_post_nov_sacamp),nanstd(norm_post_nov_sacamp)./sqrt(npost))
errorbar(nanmean(norm_post_rep_sacamp),nanstd(norm_post_rep_sacamp)./sqrt(npost))
hold off
xlabel('Ordinal Saccade #')
ylabel('Relative Saccade Amplitude')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEastOutside')
xlim([0 median_post_nov_sac_count+1])
title('Normalized Saccade Amplitudes')


[~,norm_p_nov] = ttest2(norm_pre_nov,norm_post_nov);
[~,norm_p_rep] = ttest2(norm_pre_rep,norm_post_rep);

subplot(2,2,4)
errorb([[mean(norm_pre_nov);mean(norm_post_nov)] [mean(norm_pre_rep);mean(norm_post_rep)]]',...
    [[std(norm_pre_nov)./sqrt(npre) ;std(norm_post_nov)./sqrt(npost)] ...
    [std(norm_pre_rep)./sqrt(npre) ;std(norm_post_rep)./sqrt(npost)]]')
hold on
if norm_p_nov < 0.05
    plot(1,mean(norm_pre_nov)+2*std(norm_pre_nov),'k*')
end
if norm_p_rep < 0.05
    plot(2,mean(norm_pre_rep)+2*std(norm_pre_rep),'k*')
end
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Relative Saccade Amplitude')
title(sprintf(['ANOVA Saccades 2-5: \n p_{lesion} = ' num2str(P_ANOVA_sacamps(2),2) ...
    ' , p_{Novelty} = ' num2str(P_ANOVA_sacamps(1),2) ', p_{inter} = ' num2str(P_ANOVA_sacamps(3),2)]))
axis square
box off
legend('Pre','Post','Location','NorthEastOutside')

%subtitle([pre_files{1}(1:2) ': Pre vs. Post ListRM'])

%%
%---Plot Pupil Diameter Over Time---%
t = 1:5:7500;
figure
hold on
dofill(t,pre_nov_pupil/1000,'r',1,60)
dofill(t,pre_rep_pupil/1000,'b',1,60)
dofill(t,post_nov_pupil/1000,'m',1,60)
dofill(t,post_rep_pupil/1000,'g',1,60)
yl = ylim;
plot([500 500],[yl(1) yl(2)],'--k')
hold off
ylim(yl)
xlim([0 7500])
xlabel('Time from Image Onset (ms)')
ylabel('Normalize Pupil Data (a.u.)')
set(gca,'Xtick',0:500:7500)
set(gca,'XtickLabel',num2cell([0:500:7500]-500))
title('Pupil Diameter')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat')


subtitle([pre_files{1}(1:2) ': Pre vs. Post ListRM'])
%% Calculate Salience at Fixation Locations..No Stats

pre_novel_sal = NaN(length(pre_files),50);
pre_repeat_sal = NaN(length(pre_files),50);
for file = 1:length(pre_files)
    load([data_dir pre_files{file}(1:8) '_' pre_files{file}(end) '-Salience.mat'])
    
    median_novel = sum(~isnan(novel_sal),2);
    median_novel(median_novel == 0) = [];%if no fixations don't count
    median_novel = round(median(median_novel));
    
    median_repeat = sum(~isnan(repeat_sal),2);
    median_repeat(median_repeat == 0) = [];%if no fixations don't count
    median_repeat = round(median(median_repeat));
    
    pre_novel_sal(file,1:median_novel) = nanmean(novel_sal(:,1:median_novel));
    pre_repeat_sal(file,1:median_repeat) = nanmean(repeat_sal(:,1:median_repeat));
end
%trim down extra fixations
pre_median_novel = round(median(sum(~isnan(pre_novel_sal),2)));
pre_median_repeat = round(median(sum(~isnan(pre_repeat_sal),2)));

post_novel_sal = NaN(length(post_files),50);
post_repeat_sal = NaN(length(post_files),50);
for file = 1:length(post_files)
    load([data_dir post_files{file}(1:8) '_' post_files{file}(end) '-Salience.mat'])
    
    median_novel = sum(~isnan(novel_sal),2);
    median_novel(median_novel == 0) = [];%if no fixations don't count
    median_novel = round(median(median_novel));
    
    median_repeat = sum(~isnan(repeat_sal),2);
    median_repeat(median_repeat == 0) = [];%if no fixations don't count
    median_repeat = round(median(median_repeat));
    
    post_novel_sal(file,1:median_novel) = nanmean(novel_sal(:,1:median_novel));
    post_repeat_sal(file,1:median_repeat) = nanmean(repeat_sal(:,1:median_repeat));
end
%trim down extra fixations
post_median_novel = round(median(sum(~isnan(post_novel_sal),2)));
post_median_repeat = round(median(sum(~isnan(post_repeat_sal),2)));

figure
hold all
errorbar(mean(pre_novel_sal(:,1:pre_median_novel)),std(pre_novel_sal(:,1:pre_median_novel))./sqrt(length(pre_files)));
errorbar(mean(pre_repeat_sal(:,1:pre_median_repeat)),std(pre_repeat_sal(:,1:pre_median_repeat))./sqrt(length(pre_files)));
errorbar(mean(post_novel_sal(:,1:post_median_novel)),std(post_novel_sal(:,1:post_median_novel))./sqrt(length(post_files)));
errorbar(mean(post_repeat_sal(:,1:post_median_repeat)),std(post_repeat_sal(:,1:post_median_repeat))./sqrt(length(post_files)));
hold off
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat')
xlabel('Ordinal Fixation #')
ylabel('Normalized Salience')
title('Salience @ Fixation Locations')
subtitle([pre_files{1}(1:2) ': Pre vs. Post ListRM'])
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
errorb([mean(pre_coverage');mean(post_coverage')]',...
    [std(pre_coverage')./sqrt(size(pre_coverage,2)); std(post_coverage')./sqrt(size(post_coverage,2))]')
hold on
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
errorb([mean(pre_coverage');mean(post_coverage')]',...
    [std(pre_coverage')./sqrt(size(pre_coverage,2)); std(post_coverage')./sqrt(size(post_coverage,2))]')
legend('Pre','Post')
set(gca,'Xtick',[1,2]);
set(gca,'XtickLabel',{'Novel','Repeat'})
legend('Pre','Post')
ylabel('Relative Coverage')
title('Relative')
ylim([0.95 1.05])
axis square

subtitle('Viewing Coverage')
%save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_ListRM_Pre_vs_Post_Image_Coverage'])
