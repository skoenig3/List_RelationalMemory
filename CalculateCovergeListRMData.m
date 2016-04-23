function CalculateCovergeListRMData(data_dir,fixationfile)
%written by Seth Konig 7/22/15. Code determines the percentage of the image
%viewed/attended to by the monkey assuming an area of attention with a
%radius of 2 dva (IOR_area).

load([data_dir fixationfile]);

imageX = 800;
imageY = 600;
fixwin = 3.5;
imgsize = imageX*imageY;

IOR_area = 48;
[rr,cc] = meshgrid(1:imageX,1:imageY);

nov_coverage = NaN(1,96); %percent coverage by image
rep_coverage = NaN(1,96); %percent coverage by image

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
if strcmpi('PW150416_2',fixationfile(1:10))
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
    
    %     nov_x = fixationstats{viewed(1,novrep)}.XY(1,nov_img_on:nov_img_off);
    %     nov_y = fixationstats{viewed(1,novrep)}.XY(2,nov_img_on:nov_img_off);
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
    
    %     rep_x = fixationstats{viewed(2,novrep)}.XY(1,rep_img_on:rep_img_off);
    %     rep_y = fixationstats{viewed(2,novrep)}.XY(2,rep_img_on:rep_img_off);
    rep_fix = fixationstats{viewed(2,novrep)}.fixations;
    rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
    
    pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
    post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
    rep_fix(:,post_img_fix) = [];
    rep_fix(:,pre_img_fix) = [];
    
    %---Remove Central Fixations that continued to Occur once image was on---%
    first_out = find((nov_fix(1,1) > imageX/2+fixwin*24 | nov_fix(1,1) <  imageX/2-24*fixwin) | ...
        (nov_fix(2,1) > imageY/2+fixwin*24 | nov_fix(2,1) < imageY/2-fixwin*24));
    nov_fix(:,1) = [];
    nov_fixtimes(:,1) = [];
    
    first_out = find((rep_fix(1,1) > imageX/2+fixwin*24 | rep_fix(1,1) <  imageX/2-24*fixwin) | ...
        (rep_fix(2,1) > imageY/2+fixwin*24 | rep_fix(2,1) < imageY/2-fixwin*24));
    rep_fix(:,1) = [];
    rep_fixtimes(:,1) = [];
    
    
    nov_img = zeros(imageY,imageX);
    for i = 1:size(nov_fix,2)
        C = sqrt((rr-nov_fix(1,i)).^2+(cc-nov_fix(2,i)).^2)<=IOR_area;
        nov_img(C) = 1;
    end
    
    rep_img = zeros(imageY,imageX);
    for i = 1:size(rep_fix,2)
        C = sqrt((rr-rep_fix(1,i)).^2+(cc-rep_fix(2,i)).^2)<=IOR_area;
        rep_img(C) = 1;
    end
    
    %     figure
    %     subplot(1,2,1)
    %     hold on
    %     imagesc(nov_img);
    %     plot(nov_x,nov_y);
    %     hold off
    %     axis off
    %
    %     subplot(1,2,2)
    %     hold on
    %     imagesc(rep_img);
    %     plot(rep_x,rep_y);
    %     hold off
    %     axis off
    
    nov_coverage(novrep) = sum(sum(nov_img))/imgsize;
    rep_coverage(novrep) = sum(sum(rep_img))/imgsize;
end

save([data_dir fixationfile(1:8) '_' fixationfile(10) '-Coverage.mat'],'nov_coverage','rep_coverage')
