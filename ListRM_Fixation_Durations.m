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
    
end
