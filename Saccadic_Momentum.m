clear,clc
fixwin = 3.5; %size of the fixation window/2. Was a width of 7 dva
imageX = 800;
imageY = 600;

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Cortex Data\';
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Image Sets\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Figures\';


pre_files = {'PW150416.4','PW150417.2','PW150420.2','PW150421.2',...
    'PW150422.2','PW150423.2','PW150424.2','PW150427.2',...
    'PW150428.2','PW150429.2','PW150430.2','PW150501.2',...
    'PW150504.2','PW150505.2','PW150506.2','PW150507.2',...
    'PW150511.2','TO150513.2','TO150514.2','TO150515.2','TO150518.2','TO150519.2',...
    'TO150520.2','TO150521.2','TO150522.2','TO150526.2','TO150527.2',...
    'TO150528.2','TO150529.2','TO150601.2','TO150602.2','TO150603.2',...
    'TO150604.2','TO150605.2','RR150423.2','RR150424.2','RR150427.2','RR150428.2','RR150429.2',...
    'RR150430.2','RR150501.2','RR150504.2','RR150505.2','RR150506.2',...
    'RR150507.2','RR150508.2','RR150511.2','RR150512.2','RR150513.2',...
    'RR150515.2','RR150518.2','RR150519.2','RR150520.2','RR150521.2',...
    'MF170111.2','MF170112.2','MF170117.2','MF170118.2',...
    'MF170119.2','MF170120.2','MF170123.2','MF170126.2',...
    'MF170130.2','MF170131.2','MF170201.2','MF170202.2',...
    'MF170209.2','MF170210.2','MF170213.2'};


nov_saccade_change_angles = NaN(7500,50);
nov_saccade_change_amps = NaN(7500,50);
nov_fixation_durations = NaN(7500,50);
nov_saccade_angles = NaN(7500,50);

rep_saccade_change_angles = NaN(7500,50);
rep_saccade_change_amps = NaN(7500,50);
rep_fixation_durations = NaN(7500,50);
rep_saccade_angles = NaN(7500,50);

img_count = 1;
for file = 1:length(pre_files);
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
    
    %1st image trial was displayed with the wrong timing file, then when fixed
    %started the task over again. The first image therefore is not novel and thus
    %we will not be analyzing data for this image on this session.
    if strcmpi('PW150416_2',pre_files{file}(1:10))
        viewed(:,1) = 0;
        imgdur(:,1) = 0;
    end
    
    %Let's not even look at images in which the monkey looked away more than the
    %image was displayed for. Also should take care of EOG overflow trials
    viewed(viewed & img_dur > 14000) = 0; %1000 ms buffer for accidental lookaways and cortex timing
    viewed(:,viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    
    for novrep = 1:size(viewed,2);
        %---Grab Important Vairables---%
        %for novel images
        nov_allval = per(viewed(1,novrep)).allval;
        nov_alltim = per(viewed(1,novrep)).alltim;
        nov_img_on = nov_alltim(nov_allval == 23)-nov_alltim(nov_allval == 100);%image on relative to eye data start
        nov_img_off = nov_alltim(nov_allval == 24)-nov_alltim(nov_allval == 100);%image on relative to eye data start
        
        nov_xy = fixationstats{viewed(1,novrep)}.XY;
        nov_xy = nov_xy(:,nov_img_on:nov_img_off);
        if (nov_img_off-nov_img_on) > 10500
            continue
        end
        
        nov_fix = fixationstats{viewed(1,novrep)}.fixations;
        nov_fixtimes = fixationstats{viewed(1,novrep)}.fixationtimes;
        nov_sactimes = fixationstats{viewed(1,novrep)}.saccadetimes;
        
        pre_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
        post_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
        nov_fix(:,post_img_fix) = [];
        nov_fix(:,pre_img_fix) = [];
        nov_fixtimes(:,post_img_fix) = [];
        nov_fixtimes(:,pre_img_fix ) = [];
        pre_img_sac = find(nov_sactimes(1,:) <= nov_img_on);
        post_img_sac = find(nov_sactimes(2,:) > nov_img_off);
        nov_sactimes(:,post_img_sac) = [];
        nov_sactimes(:,pre_img_sac ) = [];
        
        
        for f = 2:size(nov_fix,2)
            
            prior_sac = find(nov_sactimes(2,:) == nov_fixtimes(1,f)-1);%prior saccade
            post_sac =  find(nov_sactimes(1,:) == nov_fixtimes(2,f)+1);%next saccade
            if isempty(prior_sac) || isempty(post_sac)
                continue
            end
            prior_fix = find(nov_sactimes(1,prior_sac)-1 == nov_fixtimes(2,:));
            post_fix = find(nov_sactimes(2,post_sac)+1 == nov_fixtimes(1,:));
            if isempty(prior_fix) || isempty(post_fix)
                continue
            end
            
            this_fix = nov_fix(:,f);
            prior_fix = nov_fix(:,prior_fix);
            post_fix = nov_fix(:,post_fix);
            
            prior_amp = sqrt(sum((prior_fix-this_fix).^2));
            post_amp = sqrt(sum((post_fix-this_fix).^2));
            
            if prior_amp < 48 || post_amp < 48 %too small
                continue
            end
            
            pre_sac_angle = atan2d(this_fix(2)-prior_fix(2),this_fix(1)-prior_fix(1));
            post_sac_angle = atan2d(post_fix(2)-this_fix(2),post_fix(1)-this_fix(1));
            
            nov_saccade_change_angles(img_count,f) = post_sac_angle-pre_sac_angle;
            nov_saccade_change_amps(img_count,f) = post_amp-prior_amp;
            nov_fixation_durations(img_count,f) = diff(nov_fixtimes(:,f))+1;
            nov_saccade_angles(img_count,f) = pre_sac_angle;%just do prior for now
        end
        
        
        %---Grab Important Vairables---%
        %for novel images
        rep_allval = per(viewed(2,novrep)).allval;
        rep_alltim = per(viewed(2,novrep)).alltim;
        rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        
        rep_xy = fixationstats{viewed(2,novrep)}.XY;
        rep_xy = rep_xy(:,rep_img_on:rep_img_off);
        if (rep_img_off-rep_img_on) > 10500
            continue
        end
        
        rep_fix = fixationstats{viewed(2,novrep)}.fixations;
        rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
        rep_sactimes = fixationstats{viewed(2,novrep)}.saccadetimes;
        
        pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        rep_fix(:,post_img_fix) = [];
        rep_fix(:,pre_img_fix) = [];
        rep_fixtimes(:,post_img_fix) = [];
        rep_fixtimes(:,pre_img_fix ) = [];
        pre_img_sac = find(rep_sactimes(1,:) <= rep_img_on);
        post_img_sac = find(rep_sactimes(2,:) > rep_img_off);
        rep_sactimes(:,post_img_sac) = [];
        rep_sactimes(:,pre_img_sac ) = [];
        
        
        for f = 2:size(rep_fix,2)
            
            prior_sac = find(rep_sactimes(2,:) == rep_fixtimes(1,f)-1);%prior saccade
            post_sac =  find(rep_sactimes(1,:) == rep_fixtimes(2,f)+1);%next saccade
            if isempty(prior_sac) || isempty(post_sac)
                continue
            end
            prior_fix = find(rep_sactimes(1,prior_sac)-1 == rep_fixtimes(2,:));
            post_fix = find(rep_sactimes(2,post_sac)+1 == rep_fixtimes(1,:));
            if isempty(prior_fix) || isempty(post_fix)
                continue
            end
            
            this_fix = rep_fix(:,f);
            prior_fix = rep_fix(:,prior_fix);
            post_fix = rep_fix(:,post_fix);
            
            prior_amp = sqrt(sum((prior_fix-this_fix).^2));
            post_amp = sqrt(sum((post_fix-this_fix).^2));
            
            if prior_amp < 48 || post_amp < 48 %too small
                continue
            end
            
            pre_sac_angle = atan2d(this_fix(2)-prior_fix(2),this_fix(1)-prior_fix(1));
            post_sac_angle = atan2d(post_fix(2)-this_fix(2),post_fix(1)-this_fix(1));
            
            rep_saccade_change_angles(img_count,f) = post_sac_angle-pre_sac_angle;
            rep_saccade_change_amps(img_count,f) = post_amp-prior_amp;
            rep_fixation_durations(img_count,f) = diff(rep_fixtimes(:,f))+1;
            rep_saccade_angles(img_count,f) = pre_sac_angle;%just do prior for now
        end
        
        img_count = img_count +1;
    end
end
%%
nov_saccade_change_amps = nov_saccade_change_amps/24;%convert to dva
rep_saccade_change_amps = rep_saccade_change_amps/24;%convert to dva

%%
hist_angle_bin = 3;
degrees = [-180:hist_angle_bin:180];

nov_change_amps = nov_saccade_change_amps;
nov_change_amps(isnan(nov_change_amps)) = [];

rep_change_amps = rep_saccade_change_amps;
rep_change_amps(isnan(rep_change_amps)) = [];

nov_sac_angles = nov_saccade_angles;
nov_sac_angles(isnan(nov_sac_angles)) = [];

rep_sac_angles = rep_saccade_angles;
rep_sac_angles(isnan(rep_sac_angles)) = [];

nov_fixdurs = nov_fixation_durations;
nov_fixdurs(isnan(nov_fixdurs)) = [];

rep_fixdurs = rep_fixation_durations;
rep_fixdurs(isnan(rep_fixdurs)) = [];

nov_angles = nov_saccade_change_angles(:);
nov_angles(isnan(nov_angles)) = [];
nov_angles(nov_angles < 0) = nov_angles(nov_angles < 0)+360;
nov_angles = nov_angles-180;

rep_angles = rep_saccade_change_angles(:);
rep_angles(isnan(rep_angles)) = [];
rep_angles(rep_angles < 0) = rep_angles(rep_angles < 0)+360;
rep_angles = rep_angles-180;


%%
nov_hist_angle_dist = NaN(1,length(degrees)-1);
rep_hist_angle_dist = NaN(1,length(degrees)-1);
nov_hist_saccade_angles =  NaN(1,length(degrees)-1);
rep_hist_saccade_angles =  NaN(1,length(degrees)-1);

for d = 2:length(degrees)
    these_dirs = nov_angles >= degrees(d-1) & nov_angles < degrees(d);
    nov_hist_angle_dist(d-1) = sum(these_dirs);
    
    these_dirs = nov_sac_angles >= degrees(d-1) & nov_sac_angles < degrees(d);
    nov_hist_saccade_angles(d-1) = sum(these_dirs);
    
    these_dirs = rep_angles >= degrees(d-1) & rep_angles < degrees(d);
    rep_hist_angle_dist(d-1) = sum(these_dirs);
    
    these_dirs = rep_sac_angles >= degrees(d-1) & rep_sac_angles < degrees(d);
    rep_hist_saccade_angles(d-1) = sum(these_dirs);
    
    
end
nov_hist_angle_dist = nov_hist_angle_dist/sum(nov_hist_angle_dist);
rep_hist_angle_dist = rep_hist_angle_dist/sum(rep_hist_angle_dist);

nov_hist_saccade_angles = nov_hist_saccade_angles/sum(nov_hist_saccade_angles);
rep_hist_saccade_angles = rep_hist_saccade_angles/sum(rep_hist_saccade_angles);
%%

hist_angle_bin = 8;
degrees2 = [0:hist_angle_bin:180];

nov_fix_dur_direction = NaN(1,length(degrees2)-1);
rep_fix_dur_direction = NaN(1,length(degrees2)-1);
for d = 2:length(degrees2)
    these_dirs = abs(nov_angles) >= degrees2(d-1) & abs(nov_angles) < degrees2(d);
    nov_fix_dur_direction(d-1) = mean(nov_fixdurs(these_dirs));
        
    these_dirs = abs(rep_angles) >= degrees2(d-1) & abs(rep_angles) < degrees2(d);
    rep_fix_dur_direction(d-1) = mean(rep_fixdurs(these_dirs));

end
%%
amp_bin = 2;
amp_bins = [-6:amp_bin:8];

nov_amp_dir_dur = NaN(1,length(amp_bins));
rep_amp_dir_dur = NaN(1,length(amp_bins));

these_dirs = abs(nov_angles) >165 & abs(nov_angles) <= 180;
these_fix_durs = nov_fixdurs(these_dirs);
these_sac_amps = nov_change_amps(these_dirs);

for apb = 1:length(amp_bins)-1
    these_amps = these_sac_amps < amp_bins(apb+1) & these_sac_amps >= amp_bins(apb);
    nov_amp_dir_dur(apb) = mean(these_fix_durs(these_amps));
end


these_dirs = abs(rep_angles) >165 & abs(rep_angles) <= 180;
these_fix_durs = rep_fixdurs(these_dirs);
these_sac_amps = rep_change_amps(these_dirs);

for apb = 1:length(amp_bins)-1
    these_amps = these_sac_amps < amp_bins(apb+1) & these_sac_amps >= amp_bins(apb);
    rep_amp_dir_dur(apb) = mean(these_fix_durs(these_amps));
end
    

%%
figure

subplot(2,2,1)
polarplot(degrees*pi/180,[nov_hist_angle_dist nov_hist_angle_dist(1)])
hold on
polarplot(degrees*pi/180,[rep_hist_angle_dist rep_hist_angle_dist(1)])
hold off
title('Distribution of \DeltaSaccade Angle')

subplot(2,2,3)
polarplot(degrees*pi/180,[nov_hist_saccade_angles nov_hist_saccade_angles(1)])
hold on
polarplot(degrees*pi/180,[rep_hist_saccade_angles rep_hist_saccade_angles(1)])
hold off
title('Distribution of Saccade Angles')

subplot(2,2,2)
plot(degrees2(1:end-1),nov_fix_dur_direction)
hold on
plot(degrees2(1:end-1),rep_fix_dur_direction-(rep_fix_dur_direction(1)-nov_fix_dur_direction(1)))
hold off
xlabel('\DeltaSaccade Angle')
ylabel('Fixation Duration (ms)')
xlim([0 180])
box off
title('Saccadic Momentum')
ylim([160 210])


subplot(2,2,4)
plot(amp_bins(1:end-1),nov_amp_dir_dur(1:end-1))
hold on
plot(amp_bins(1:end-1),rep_amp_dir_dur(1:end-1)-(mean(rep_amp_dir_dur(1:end-1))-mean(nov_amp_dir_dur(1:end-1))))
hold off
ylabel('Fixation Duration')
xlabel('\DeltaSaccade Amplitude')
title('IOR for \DeltaSaccade Angle = 180 +/- 30')
box off
xlim([-6 6])

subtitle('Novel (Blue) & Repeat (Red)')