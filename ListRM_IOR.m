function ListRM_IOR(data_dir,cortex_file,image_dir,figure_dir)
% modified from SalienceIOR Seth Konig 4/7/2017

% function determines the rate of return fixations, the time between return
% fixations, time within trial of return, and salience at returned
% location.

% Inputs:
%   cortex_file: Fixations extracted from cortex e.g. MP120606_1-fixation.mat
%   ImageX,ImageY: x and y dimensions of images
%   Pairings: take all pairing or only closest unique pairing between first
%   and second fixation
%Outputs:
%   A .mat file named [cortex_file(1:end-13) '-IOR']
%   containg saccade statistics. See variable statvariablenames for
%   detailed explanation of variables in .mat file.


imageX = 800;
imageY = 600;
distancethreshold = [0  24 48 72 96  120  144 168 200 400;...%in pixels 24 pixels/dva
    24 48 72 96 120 144  168 200 400 1000];%in pixels 24 pixels/dva

min_middle_distance=10*24;%so at least 1 5 dva saccade away from area
min_num_fix = 10;

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
    img_dur(:,1) = 0;
end

%fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal fixdist fix1dur fix2dur fixnum1 fixnum2 trialnum
nov_returnfixsal = cell(1,size(distancethreshold,2));
rep_returnfixsal = cell(1,size(distancethreshold,2));
nov_count = ones(1,size(distancethreshold,2));
rep_count = ones(1,size(distancethreshold,2));
for d = 1:size(distancethreshold,2)
    nov_returnfixsal{d} = NaN(500,14);
    rep_returnfixsal{d} = NaN(5000,14);
end

for novrep = 1:size(viewed,2);
    if any(viewed(:,novrep) == 0)
        continue
    end
    
    %---Grab Important Vairables---%
    %for novel images
    nov_allval = per(viewed(1,novrep)).allval;
    nov_alltim = per(viewed(1,novrep)).alltim;
    nov_img_on = nov_alltim(nov_allval == 23)-nov_alltim(nov_allval == 100);%image on relative to eye data start
    nov_img_off = nov_alltim(nov_allval == 24)-nov_alltim(nov_allval == 100);%image on relative to eye data start
    
%     nov_x = fixationstats{viewed(1,novrep)}.XY(1,nov_img_on:nov_img_off);
%     nov_y = fixationstats{viewed(1,novrep)}.XY(2,nov_img_on:nov_img_off);
%     nov_pupil = pupildata{viewed(1,novrep)}(round(nov_img_on/5):round(nov_img_off/5));
    
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
%     rep_pupil = pupildata{viewed(2,novrep)}(round(rep_img_on/5):round(rep_img_off/5));
    
    rep_fix = fixationstats{viewed(2,novrep)}.fixations;
    rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
    pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
    post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
    
    rep_fix(:,post_img_fix) = [];
    rep_fix(:,pre_img_fix) = [];
    rep_fixtimes(:,post_img_fix) = [];
    rep_fixtimes(:,pre_img_fix) = [];
    
    
%     %---Find When the monkey Blinked,when monkey looked away---%
%     %pupil values at 0 diameter
%     [~,~,nov_time_out] = findbad_eye_ind(nov_pupil,nov_x,1000);
%     [~,~,rep_time_out] = findbad_eye_ind(rep_pupil,rep_x,1000);
%     
%     
%     %---Calculate Fixation Durations and Saccade Amplitudes---%
%     if ~isnan(nov_time_out)
%         post_attention_fix = find(nov_fixtimes(1,:) > nov_time_out);
%         nov_fixtimes(:,post_attention_fix) = [];
%         nov_fix(:,post_attention_fix) = [];
%     end
%     
%     if ~isnan(nov_time_out)
%         post_attention_fix = find(rep_fixtimes(1,:) > rep_time_out);
%         rep_fixtimes(:,post_attention_fix) = [];
%         rep_fix(:,post_attention_fix) = [];
%     end
    
    
    %---Remove fixations if more than 50
    if size(nov_fix,2) > 50
        nov_fix = nov_fix(:,1:50);
    end
    if size(rep_fix,2) > 50
        rep_fix = rep_fix(:,1:50);
    end
    
    %---Load fullmap---%
    if imgnum(novrep) < 10
        load([image_dir2 'S' setnum 'I0' num2str(imgnum(novrep)) '-saliencemap.mat'],'fullmap')
    else
        load([image_dir2 'S' setnum 'I' num2str(imgnum(novrep)) '-saliencemap.mat'],'fullmap')
    end
    
    %----Determine Novel Fixation Pairs---%
    if size(nov_fix,2) > min_num_fix
        N=size(nov_fix,2);
        [x,y]=meshgrid(1:N);
        i=find(ones(N)-eye(N)); %forms nov_pairs except for self-pairing
        i=[x(i), y(i)];
        i(i(:,1) > i(:,2),:) = []; %repeat nov_pairs
        i(i(:,1)+1 == i(:,2),:) = []; %removes consecutive in time nov_pairs
        dist =sqrt((nov_fix(1,i(:,1))-nov_fix(1,i(:,2))).^2 +...
            (nov_fix(2,i(:,1))-nov_fix(2,i(:,2))).^2);
        wnov_count = 1;
        nov_pairs = NaN(ceil(size(nov_fix,2)/2),3);
        while ~isempty(i);
            [minn,mind] = min(dist);
            minn = minn(1); %if multiple just take the first
            mind = mind(1); %if multiple just take the first
            
            middlefixes = i(mind,1)+1:i(mind,2)-1;
            middist = sqrt((nov_fix(1,i(mind,1))-nov_fix(1,middlefixes)).^2+...
                (nov_fix(2,i(mind,2))-nov_fix(2,middlefixes)).^2);
            
            if any(middist >= min_middle_distance) %so at least 1 larger saccade away from area
                tind = find(minn > distancethreshold(1,:) & minn <= distancethreshold(2,:));
                if ~isempty(tind);
                    nov_pairs(wnov_count,:) = [i(mind,:) tind];
                end
                %remove all instances of fixations with prior or return number
                [rmvind1,~] = find(i(:,1) == i(mind,1));
                [rmvind2,~] = find(i(:,2) == i(mind,1));
                [rmvind3,~] = find(i(:,1) == i(mind,2));
                [rmvind4,~] = find(i(:,2) == i(mind,2));
                rmvind = [rmvind1; rmvind2; rmvind3; rmvind4];
                rmvind = unique(rmvind);
                i(rmvind,:) = [];
                dist(rmvind) = [];
                wnov_count = wnov_count+1;
            else %then remove pair
                dist(mind) = [];
                i(mind,:) = [];
            end
        end
        nov_pairs(isnan(nov_pairs(:,1)),:) = [];
        
        for ip = 1:size(nov_pairs,1);
            spot = [ceil(nov_fix(:,nov_pairs(ip,1))) ceil(nov_fix(:,nov_pairs(ip,2)))];
            spot(2,:) = imageY-spot(2,:); %location
            spott = [(nov_fixtimes(1,nov_pairs(ip,1))+nov_fixtimes(2,nov_pairs(ip,1)))/2 ...
                (nov_fixtimes(1,nov_pairs(ip,2))+nov_fixtimes(2,nov_pairs(ip,2)))/2]; %time
            dist = sqrt((spot(1,1)-spot(1,2))^2+(spot(2,1)-spot(2,2))^2);
            spot(spot < 1) = 1;
            spot(1,spot(1,:) > imageX) = imageX;
            spot(2,spot(2,:) > imageY) = imageY;
            
            nov_returnfixsal{nov_pairs(ip,3)}(nov_count(nov_pairs(ip,3)),:) = [...
                spot(1,1) spot(2,1) spott(1) fullmap(spot(2,1),spot(1,1))...
                spot(1,2) spot(2,2) spott(2) fullmap(spot(2,2),spot(1,2))...
                dist diff(nov_fixtimes(:,nov_pairs(ip,1)))+1 ...
                diff(nov_fixtimes(:,nov_pairs(ip,2)))+1 nov_pairs(ip,1) nov_pairs(ip,2) novrep];
            %fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal fixdist
            %fix1dur fix2dur fixnum1 fixnum2 imgnum
            nov_count(nov_pairs(ip,3)) = nov_count(nov_pairs(ip,3))+1;
        end
    end
    
    %----Determine Repeat Fixation Pairs---%
    if size(rep_fix,2) > min_num_fix
        N=size(rep_fix,2);
        [x,y]=meshgrid(1:N);
        i=find(ones(N)-eye(N)); %forms rep_pairs except for self-pairing
        i=[x(i), y(i)];
        i(i(:,1) > i(:,2),:) = []; %repeat rep_pairs
        i(i(:,1)+1 == i(:,2),:) = []; %removes consecutive in time rep_pairs
        dist =sqrt((rep_fix(1,i(:,1))-rep_fix(1,i(:,2))).^2 +...
            (rep_fix(2,i(:,1))-rep_fix(2,i(:,2))).^2);
        wrep_count = 1;
        rep_pairs = NaN(ceil(size(rep_fix,2)/2),3);
        while ~isempty(i);
            [minn,mind] = min(dist);
            minn = minn(1); %if multiple just take the first
            mind = mind(1); %if multiple just take the first
            
            middlefixes = i(mind,1)+1:i(mind,2)-1;
            middist = sqrt((rep_fix(1,i(mind,1))-rep_fix(1,middlefixes)).^2+...
                (rep_fix(2,i(mind,2))-rep_fix(2,middlefixes)).^2);
            
            if any(middist >= min_middle_distance) %so at least 1 larger saccade away from area
                tind = find(minn > distancethreshold(1,:) & minn <= distancethreshold(2,:));
                if ~isempty(tind);
                    rep_pairs(wrep_count,:) = [i(mind,:) tind];
                end
                %remove all instances of fixations with prior or return number
                [rmvind1,~] = find(i(:,1) == i(mind,1));
                [rmvind2,~] = find(i(:,2) == i(mind,1));
                [rmvind3,~] = find(i(:,1) == i(mind,2));
                [rmvind4,~] = find(i(:,2) == i(mind,2));
                rmvind = [rmvind1; rmvind2; rmvind3; rmvind4];
                rmvind = unique(rmvind);
                i(rmvind,:) = [];
                dist(rmvind) = [];
                wrep_count = wrep_count+1;
            else %then remove pair
                dist(mind) = [];
                i(mind,:) = [];
            end
        end
        rep_pairs(isnan(rep_pairs(:,1)),:) = [];
        
        for ip = 1:size(rep_pairs,1);
            spot = [ceil(rep_fix(:,rep_pairs(ip,1))) ceil(rep_fix(:,rep_pairs(ip,2)))];
            spot(2,:) = imageY-spot(2,:); %location
            spott = [(rep_fixtimes(1,rep_pairs(ip,1))+rep_fixtimes(2,rep_pairs(ip,1)))/2 ...
                (rep_fixtimes(1,rep_pairs(ip,2))+rep_fixtimes(2,rep_pairs(ip,2)))/2]; %time
            dist = sqrt((spot(1,1)-spot(1,2))^2+(spot(2,1)-spot(2,2))^2);
            spot(spot < 1) = 1;
            spot(1,spot(1,:) > imageX) = imageX;
            spot(2,spot(2,:) > imageY) = imageY;
            
            rep_returnfixsal{rep_pairs(ip,3)}(rep_count(rep_pairs(ip,3)),:) = [...
                spot(1,1) spot(2,1) spott(1) fullmap(spot(2,1),spot(1,1))...
                spot(1,2) spot(2,2) spott(2) fullmap(spot(2,2),spot(1,2))...
                dist diff(rep_fixtimes(:,rep_pairs(ip,1)))+1 ...
                diff(rep_fixtimes(:,rep_pairs(ip,2)))+1 rep_pairs(ip,1) rep_pairs(ip,2) novrep];
            %fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal fixdist
            %fix1dur fix2dur fixnum1 fixnum2 imgnum
            rep_count(rep_pairs(ip,3)) = rep_count(rep_pairs(ip,3))+1;
        end
    end
end
%---remove excess NaNs---%
nov_returnfixsal = laundry(nov_returnfixsal);
rep_returnfixsal = laundry(rep_returnfixsal);


IORvariablenames = {
    'returnfixsal: [  %fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal...';
    'fixdist fix1dur fix2dur fixnum1 fixnum2 imagnum]';
    };

save([data_dir cortex_file(1:8) '-IOR'],'nov_returnfixsal','rep_returnfixsal',...
    'IORvariablenames','distancethreshold')
end