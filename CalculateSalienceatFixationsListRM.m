function CalculateSalienceatFixationsListRM(cortex_file,data_dir,image_dir,imageX,imageY)
%written by Seth Konig 9/8/15
%code calculates salience at fixation locations for ListRM image sets


novel_sal = NaN(96,50); %Salience at fixations location for novel images
repeat_sal = NaN(96,50);%Salience at fixations location for repeat images
random_sal = NaN(96,50); %Salience at random locations

load([data_dir cortex_file(1:8) '_' cortex_file(end) '-fixation.mat'])
setnum = item_set(7:8);
image_dir2 = [image_dir 'LRM' setnum '\'];

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
if strcmpi('PW150416_2',cortex_file(1:10))
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
    pre_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
    post_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
    
    nov_fix(:,post_img_fix) = [];
    nov_fix(:,pre_img_fix) = [];
    nov_fixtimes(:,post_img_fix) = [];
    nov_fixtimes(:,pre_img_fix ) = [];
    
    
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
    pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
    post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
    
    rep_fix(:,post_img_fix) = [];
    rep_fix(:,pre_img_fix) = [];
    rep_fixtimes(:,post_img_fix) = [];
    rep_fixtimes(:,pre_img_fix) = [];
    
    
    %---Find When the monkey Blinked,when monkey looked away---%
    %pupil values at 0 diameter
    [~,~,nov_time_out] = findbad_eye_ind(nov_pupil,nov_x,1000);
    [~,~,rep_time_out] = findbad_eye_ind(rep_pupil,rep_x,1000);
    
    
    %---Calculate Fixation Durations and Saccade Amplitudes---%
    if ~isnan(nov_time_out)
        post_attention_fix = find(nov_fixtimes(1,:) > nov_time_out);
        nov_fixtimes(:,post_attention_fix) = [];
        nov_fix(:,post_attention_fix) = [];
    end
    
    if ~isnan(nov_time_out)
        post_attention_fix = find(rep_fixtimes(1,:) > rep_time_out);
        rep_fixtimes(:,post_attention_fix) = [];
        rep_fix(:,post_attention_fix) = [];
    end
    
    
    %---Remove fixations if more than 50
    if size(nov_fix,2) > 50
        nov_fix = nov_fix(:,1:50);
    end
    if size(rep_fix,2) > 50
        rep_fix = rep_fix(:,1:50);
    end
    
    %---Load saliencemap---% 
    if imgnum(novrep) < 10
        load([image_dir2 'S' setnum 'I0' num2str(imgnum(novrep)) '-saliencemap.mat'],'fullmap')
    else
        load([image_dir2 'S' setnum 'I' num2str(imgnum(novrep)) '-saliencemap.mat'],'fullmap')
    end
    
    %---Calculate Observed Salience at fixation locations---%
    for fix = 1:size(nov_fix,2) %for novel
        fixx = round(nov_fix(1,fix));
        fixx(fixx < 1) = 1;
        fixx(fixx > imageX) = imageX;
        fixy = imageY-round(nov_fix(2,fix)); %flip since salience map is a matrix
        fixy(fixy < 1) = 1;
        fixy(fixy > imageY) = imageY;
        novel_sal(imgnum(novrep),fix) = fullmap(fixy,fixx);
    end
    
    for fix = 1:size(rep_fix,2) %for repeat
        fixx = round(rep_fix(1,fix));
        fixx(fixx < 1) = 1;
        fixx(fixx > imageX) = imageX;
        fixy = imageY-round(rep_fix(2,fix)); %flip since salience map is a matrix
        fixy(fixy < 1) = 1;
        fixy(fixy > imageY) = imageY;
        repeat_sal(imgnum(novrep),fix) = fullmap(fixy,fixx);
    end
    
    avg_num_fix = (size(nov_fix,2)+size(rep_fix,2))/2;
    for fix = 1:avg_num_fix %for random locations
        x = randi(imageX-1)+1;
        y = randi(imageY-1)+1;
        random_sal(imgnum(novrep),fix) = fullmap(y,x);
    end
    
end

save([data_dir cortex_file(1:8) '_' cortex_file(end) '-Salience.mat'],...
    'novel_sal','repeat_sal','random_sal')