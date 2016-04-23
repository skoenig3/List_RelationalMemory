function DetermineListRMSpatialSimilarity(data_dir,fixationfile)
%Written 8/26/15 by Seth Konig. Code based on
%DetermineListSQCortexKLdivergence.m and CalculateCoverageListRMData.m
%Code calculates KLdivergence values aka similarity in fixation locations
%during novel and familiar presentation of an image.

imageX = 800;
imageY = 600;

IOR_area = 48;
[rr,cc] = meshgrid(1:imageX,1:imageY);

load([data_dir fixationfile]);

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
if strcmpi('PW150416_2',fixationfile(1:10))
    viewed(:,1) = 0;
    imgdur(:,1) = 0;
end


%Let's not even look at images in which the monkey looked away more than the
%image was displayed for. Also should take care of EOG overflow trials
viewed(viewed & img_dur > 14000) = 0; %1000 ms buffer for accidental lookaways and cortex timing
imgnum(viewed(1,:) == 0 | viewed(2,:) == 0) = [];
viewed(:,viewed(1,:) == 0 | viewed(2,:) == 0) = [];

%get random novel/repeats pairs but pairs can't be for the same image
ind = 1:size(viewed,2);
for trys = 1:100;
    randind = 1:size(viewed,2);
    while any(randind == ind)
        randind = ind(randperm(size(viewed,2)));
    end
end

viewedrand = [viewed(1,:);viewed(randind)];
imgnumrand = [imgnum;imgnum(randind)];

KLnorm = NaN(96,6); %Observed Data
KLshuff = NaN(96,6);%Random Shuffled Data
KLshuff2 = NaN(96,6); %Chance levels calculated by Random pairing of novel and repeat imagse.

ObservedSimilarity = NaN(1,96); %Observed
ShuffledSimilarity = NaN(1,96); %Random Shuffled
ShuffledSimilarity2 = NaN(1,96); %Chance levels calculated by Random pairing of novel and repeat imagse.

similiarity_num_fix = NaN(2,96);%number of fixations used for novel (row1) and repeat (row2)
similiarity_num_fix2 = NaN(2,96);%number of fixations used for novel (row1) and repeat (row2) for ShuffledSimilarity2

for novrep = 1:size(viewed,2);
    
    %---Grab Eye Data---%
    %for novel image
    nov_allval = per(viewed(1,novrep)).allval;
    nov_alltim = per(viewed(1,novrep)).alltim;
    nov_img_on = nov_alltim(nov_allval == 23)-nov_alltim(nov_allval == 100);%image on relative to eye data start
    nov_img_off = nov_alltim(nov_allval == 24)-nov_alltim(nov_allval == 100);%image on relative to eye data start
    
    nov_fixations = fixationstats{viewed(1,novrep)}.fixations;
    nov_fixtimes = fixationstats{viewed(1,novrep)}.fixationtimes;
    
    nov_fixations(:,nov_fixtimes(2,:) > nov_img_off) = []; %post image fixations
    nov_fixations(:,nov_fixtimes(1,:) <= nov_img_on) = []; %pre image fixations
    nov_fixations = round(nov_fixations);
    
    %make sure fixations are within image borders
    nov_fixations(1,nov_fixations(1,:) < 1) = 1;
    nov_fixations(1,nov_fixations(1,:) > imageX) = imageX;
    nov_fixations(2,nov_fixations(2,:) < 1) = 1;
    nov_fixations(2,nov_fixations(2,:) > imageY) = imageY;
    
    
    %for repeat image
    rep_allval = per(viewed(2,novrep)).allval;
    rep_alltim = per(viewed(2,novrep)).alltim;
    rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
    rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
    
    rep_fixations = fixationstats{viewed(2,novrep)}.fixations;
    rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
    
    rep_fixations(:,rep_fixtimes(2,:) > rep_img_off) = []; %post image fixations
    rep_fixations(:,rep_fixtimes(1,:) <= rep_img_on) = []; %pre image fixations
    rep_fixations = round(rep_fixations);
    
    %make sure fixations are within image borders
    rep_fixations(1,rep_fixations(1,:) < 1) = 1;
    rep_fixations(1,rep_fixations(1,:) > imageX) = imageX;
    rep_fixations(2,rep_fixations(2,:) < 1) = 1;
    rep_fixations(2,rep_fixations(2,:) > imageY) = imageY;
    
    
    %for random repeat image paired with novel image
    rep2_allval = per(viewedrand(2,novrep)).allval;
    rep2_alltim = per(viewedrand(2,novrep)).alltim;
    rep2_img_on = rep2_alltim(rep2_allval == 23)-rep2_alltim(rep2_allval == 100);%image on relative to eye data start
    rep2_img_off = rep2_alltim(rep2_allval == 24)-rep2_alltim(rep2_allval == 100);%image on relative to eye data start
    
    rep2_fixations = fixationstats{viewedrand(2,novrep)}.fixations;
    rep2_fixtimes = fixationstats{viewedrand(2,novrep)}.fixationtimes;
    
    rep2_fixations(:,rep2_fixtimes(2,:) > rep2_img_off) = []; %post image fixations
    rep2_fixations(:,rep2_fixtimes(1,:) <= rep2_img_on) = []; %pre image fixations
    rep2_fixations = round(rep2_fixations);
    
    %make sure fixations are within image borders
    rep2_fixations(1,rep2_fixations(1,:) < 1) = 1;
    rep2_fixations(1,rep2_fixations(1,:) > imageX) = imageX;
    rep2_fixations(2,rep2_fixations(2,:) < 1) = 1;
    rep2_fixations(2,rep2_fixations(2,:) > imageY) = imageY;
    
    %generate 100 random x and y positions for random fixations
    %use same location for KL divergence and Overlap measures
    shuff_nov_fixations_x = randi(imageX-1,100)+1;
    shuff_nov_fixations_y = randi(imageY-1,100)+1;
    shuff_rep_fixations_x = randi(imageX-1,100)+1;
    shuff_rep_fixations_y = randi(imageY-1,100)+1;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Overlap Similarity Measures---%
    
    similiarity_num_fix(1,imgnum(novrep)) = size(nov_fixations,2); %number of fixations used for novel image
    similiarity_num_fix(2,imgnum(novrep)) = size(rep_fixations,2); %number of fixations used for repeat image
    similiarity_num_fix2(1,imgnumrand(1,novrep)) = size(nov_fixations,2); %number of fixations used for novel image
    similiarity_num_fix2(1,imgnumrand(2,novrep)) = size(rep2_fixations,2); %number of fixations used for random repeat image
    %these also apply to KL divergence
    
    nov_img = zeros(imageY,imageX);
    nov_img_shuff = zeros(imageY,imageX);
    for i = 1:size(nov_fixations,2)
        C = sqrt((rr-nov_fixations(1,i)).^2+(cc-nov_fixations(2,i)).^2)<=IOR_area;
        nov_img(C) = 1;
        
        C = sqrt((rr-shuff_nov_fixations_x(i)).^2+(cc-shuff_nov_fixations_y(i)).^2)<=IOR_area;
        nov_img_shuff(C) = 1;
    end
    
    rep_img = zeros(imageY,imageX);
    rep_img_shuff = zeros(imageY,imageX);
    for i = 1:size(rep_fixations,2)
        C = sqrt((rr-rep_fixations(1,i)).^2+(cc-rep_fixations(2,i)).^2)<=IOR_area;
        rep_img(C) = 1;
        
        C = sqrt((rr-shuff_rep_fixations_x(i)).^2+(cc-shuff_rep_fixations_y(i)).^2)<=IOR_area;
        rep_img_shuff(C) = 1;
    end
    
    rep2_img = zeros(imageY,imageX);
    rep2_img_shuff = zeros(imageY,imageX);
    for i = 1:size(rep2_fixations,2)
        C = sqrt((rr-rep2_fixations(1,i)).^2+(cc-rep2_fixations(2,i)).^2)<=IOR_area;
        rep2_img(C) = 1;
    end
    
    total_fix = size(nov_fixations,2)+size(rep_fixations,2);
    ObservedSimilarity(imgnum(novrep))= sum(sum(nov_img == 1 & rep_img == 1))/total_fix;
    ShuffledSimilarity(imgnum(novrep))= sum(sum(nov_img_shuff == 1 & rep_img_shuff == 1))/total_fix;
    
    total_fix2 = size(nov_fixations,2)+size(rep2_fixations,2);
    ShuffledSimilarity2(imgnumrand(1,novrep))= sum(sum(nov_img == 1 & rep2_img == 1))/total_fix2;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---KL Divergence Measures---%
    
    %storage of fixations locations
    novelfixations = cell(1,6);
    repeatfixations = cell(1,6);
    shuffled_novelfixations = cell(1,6);
    shuffled_repeatfixations = cell(1,6);
    for i = 1:6; %fill with a matrix of zeros
        novelfixations{i} = zeros(imageY,imageX);
        repeatfixations{i}=zeros(imageY,imageX);
        shuffled_novelfixations{i} = zeros(imageY,imageX);
        shuffled_repeatfixations{i} = zeros(imageY,imageX);
    end
    
    %we want to take the same number of fixations from the novel and repeat trials
    maxfixations = min(size(nov_fixations,2),size(rep_fixations,2));
    maxfixations(maxfixations > 25) = 25; %don't care if there are more than 20 fixations
    %since puting into groups of five remove the remainder of #/5,
    maxfixations = maxfixations-rem(maxfixations,5);
    if maxfixations >=5 %want at least 5 fixations
        for fixation = 1:maxfixations;
            %for all fixations 1-all
            
            nov_fix_x = nov_fixations(1,fixation);%horizontal fixation position for novel presentation
            nov_fix_y = nov_fixations(2,fixation);%vertical fixation position for novel presentation
            rep_fix_x = rep_fixations(1,fixation);%horizonal fixation position for repeat presentation
            rep_fix_y = rep_fixations(2,fixation);%vertical fixation position for repeat presentation
            
            shuff_nov_fix_x = shuff_nov_fixations_x(fixation);
            shuff_nov_fix_y = shuff_nov_fixations_y(fixation);
            shuff_rep_fix_x = shuff_rep_fixations_x(fixation);
            shuff_rep_fix_y = shuff_rep_fixations_y(fixation);
            
            novelfixations{6}(nov_fix_y,nov_fix_x) = novelfixations{6}(nov_fix_y,nov_fix_x)+1;
            repeatfixations{6}(rep_fix_y,nov_fix_x) = repeatfixations{6}(rep_fix_y,nov_fix_x)+1;
            shuffled_novelfixations{6}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{6}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
            shuffled_repeatfixations{6}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{6}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
            
            
            %put fixations in their appropriate PDF. Mark matrix with a
            %1 where there was a fixation
            if fixation <= 5
                novelfixations{1}(nov_fix_y,nov_fix_x) = novelfixations{1}(nov_fix_y,nov_fix_x)+1;
                repeatfixations{1}(rep_fix_y,nov_fix_x) = repeatfixations{1}(rep_fix_y,nov_fix_x)+1;
                shuffled_novelfixations{1}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{1}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                shuffled_repeatfixations{1}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{1}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
            elseif fixation <= 10
                novelfixations{2}(nov_fix_y,nov_fix_x) = novelfixations{2}(nov_fix_y,nov_fix_x)+1;
                repeatfixations{2}(rep_fix_y,nov_fix_x) = repeatfixations{2}(rep_fix_y,nov_fix_x)+1;
                shuffled_novelfixations{2}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{2}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                shuffled_repeatfixations{2}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{2}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
            elseif fixation <=15
                novelfixations{3}(nov_fix_y,nov_fix_x) = novelfixations{3}(nov_fix_y,nov_fix_x)+1;
                repeatfixations{3}(rep_fix_y,nov_fix_x) = repeatfixations{3}(rep_fix_y,nov_fix_x)+1;
                shuffled_novelfixations{3}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{3}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                shuffled_repeatfixations{3}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{3}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
            elseif fixation <=20
                novelfixations{4}(nov_fix_y,nov_fix_x) = novelfixations{4}(nov_fix_y,nov_fix_x)+1;
                repeatfixations{4}(rep_fix_y,nov_fix_x) = repeatfixations{4}(rep_fix_y,nov_fix_x)+1;
                shuffled_novelfixations{4}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{4}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                shuffled_repeatfixations{4}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{4}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
            elseif fixation <=25
                novelfixations{5}(nov_fix_y,nov_fix_x) = novelfixations{5}(nov_fix_y,nov_fix_x)+1;
                repeatfixations{5}(rep_fix_y,nov_fix_x) = repeatfixations{5}(rep_fix_y,nov_fix_x)+1;
                shuffled_novelfixations{5}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{5}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                shuffled_repeatfixations{5}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{5}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
            end
        end
    end
    
    for i = 1:6
        if i == 6
            if sum(sum(novelfixations{i})) >= 5 &&  sum(sum(repeatfixations{i})) >= 5
                Distance = KL_Divergence(novelfixations{i},repeatfixations{i});
                KLnorm(imgnum(novrep),i) = Distance;
                Distance = KL_Divergence(shuffled_novelfixations{i},shuffled_repeatfixations{i});
                KLshuff(imgnum(novrep),i) = Distance;
            else
                KLnorm(imgnum(novrep),i) = NaN;
                KLshuff(imgnum(novrep),i) = NaN;
            end
        else
            if sum(sum(novelfixations{i})) == 5 &&  sum(sum(repeatfixations{i})) == 5
                Distance = KL_Divergence(novelfixations{i},repeatfixations{i});
                KLnorm(imgnum(novrep),i) = Distance;
                Distance = KL_Divergence(shuffled_novelfixations{i},shuffled_repeatfixations{i});
                KLshuff(imgnum(novrep),i) = Distance;
            else
                KLnorm(imgnum(novrep),i) = NaN;
                KLshuff(imgnum(novrep),i) = NaN;
            end
        end
    end
    
    %now do KL divergence for random pairs of novel and repeat images
    
    %storage of fixations locations
    novelfixations = cell(1,6);
    repeatfixations = cell(1,6);
    for i = 1:6; %fill with a matrix of zeros
        novelfixations{i} = zeros(imageY,imageX);
        repeatfixations{i}=zeros(imageY,imageX);
    end
    
    
    %we want to take the same number of fixations from the novel and repeat trials
    maxfixations = min(size(nov_fixations,2),size(rep2_fixations,2));
    maxfixations(maxfixations > 25) = 25; %don't care if there are more than 25 fixations
    %since puting into groups of five remove the remainder of #/5,
    maxfixations = maxfixations-rem(maxfixations,5);
    if maxfixations >=5 %want at least 5 fixations
        for fixation = 1:maxfixations;
            %for all fixations 1-all
            
            nov_fix_x = nov_fixations(1,fixation);%horizontal fixation position for novel presentation
            nov_fix_y = nov_fixations(2,fixation);%vertical fixation position for novel presentation
            rep2_fix_x = rep2_fixations(1,fixation);%horizonal fixation position for repeat presentation
            rep2_fix_y = rep2_fixations(2,fixation);%vertical fixation position for repeat presentation
            
            novelfixations{6}(nov_fix_y,nov_fix_x) = novelfixations{6}(nov_fix_y,nov_fix_x)+1;
            repeatfixations{6}(rep2_fix_y,nov_fix_x) = repeatfixations{6}(rep2_fix_y,nov_fix_x)+1;

            %put fixations in their appropriate PDF. Mark matrix with a
            %1 where there was a fixation
            if fixation <= 5
                novelfixations{1}(nov_fix_y,nov_fix_x) = novelfixations{1}(nov_fix_y,nov_fix_x)+1;
                repeatfixations{1}(rep2_fix_y,nov_fix_x) = repeatfixations{1}(rep2_fix_y,nov_fix_x)+1;
            elseif fixation <= 10
                novelfixations{2}(nov_fix_y,nov_fix_x) = novelfixations{2}(nov_fix_y,nov_fix_x)+1;
                repeatfixations{2}(rep2_fix_y,nov_fix_x) = repeatfixations{2}(rep2_fix_y,nov_fix_x)+1;
            elseif fixation <=15
                novelfixations{3}(nov_fix_y,nov_fix_x) = novelfixations{3}(nov_fix_y,nov_fix_x)+1;
                repeatfixations{3}(rep2_fix_y,nov_fix_x) = repeatfixations{3}(rep2_fix_y,nov_fix_x)+1;
            elseif fixation <=20
                novelfixations{4}(nov_fix_y,nov_fix_x) = novelfixations{4}(nov_fix_y,nov_fix_x)+1;
                repeatfixations{4}(rep2_fix_y,nov_fix_x) = repeatfixations{4}(rep2_fix_y,nov_fix_x)+1;
            elseif fixation <=25
                novelfixations{5}(nov_fix_y,nov_fix_x) = novelfixations{5}(nov_fix_y,nov_fix_x)+1;
                repeatfixations{5}(rep2_fix_y,nov_fix_x) = repeatfixations{5}(rep2_fix_y,nov_fix_x)+1;
            end
        end
    end

    for i = 1:6
        if i == 6
            if sum(sum(novelfixations{i})) >= 5 &&  sum(sum(repeatfixations{i})) >= 5
                Distance = KL_Divergence(novelfixations{i},repeatfixations{i});
                KLshuff2(imgnumrand(1,novrep),i) = Distance;
            else
                KLshuff2(imgnumrand(1,novrep),i) = NaN;
            end
        else
            if sum(sum(novelfixations{i})) == 5 &&  sum(sum(repeatfixations{i})) == 5
                Distance = KL_Divergence(novelfixations{i},repeatfixations{i});
                KLshuff2(imgnumrand(1,novrep),i) = Distance;
            else
                KLshuff2(imgnumrand(1,novrep),i) = NaN;
            end
        end
    end
    
    
end

save([data_dir fixationfile(1:end-13) '-KLdivergence.mat'],'KLshuff','KLshuff2',...
    'KLnorm','ObservedSimilarity','ShuffledSimilarity','ShuffledSimilarity2',...
    'similiarity_num_fix','similiarity_num_fix','imgnum','imgnumrand')

end
