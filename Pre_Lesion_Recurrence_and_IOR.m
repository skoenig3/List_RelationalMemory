clear,clc
fixwin = 3.5; %size of the fixation window/2. Was a width of 7 dva
imageX = 800;
imageY = 600;
distance_threshold = 50;%in pixels, 24 pixels/dva
max_out_time = 100;%in ms, less than 100 ms could be blink

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Cortex Data\';
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Image Sets\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Figures\';


files = {'PW150416.4','PW150417.2','PW150420.2','PW150421.2',...
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

%define diagonals for 1-back,2-back,3-back
id0 =  eye(40);
id0 =  [id0(2:end,:); zeros(1,40)];
id1 = eye(40);
id1 = [id1(3:end,:); zeros(2,40)];
id2 = eye(40);
id2 = [id2(4:end,:); zeros(3,40)];
id3 = eye(40);
id3 = [id3(5:end,:); zeros(4,40)];

% sobel filters detect edges!
sobelx = [1     2   1;
    0     0   0;
    -1    -2  -1;];

sobely = [1     2   1;
    0     0   0;
    -1    -2  -1;];


nov_1backs = NaN(1,4000);
nov_2backs = NaN(1,4000);
nov_3backs = NaN(1,4000);

rep_1backs = NaN(1,4000);
rep_2backs = NaN(1,4000);
rep_3backs = NaN(1,4000);


all_nov_recurence_map = zeros(40,40);
all_rep_recurence_map = zeros(40,40);
all_nov_rep_recurence_map =  zeros(40,40);
all_nov_rep_count = zeros(1,40);
rep_recurence_measure = NaN(1,4000);
nov_recurence_measure =NaN(1,4000);
nov_rep_recurence_measure =NaN(1,4000);
nov_determinism = NaN(1,4000);
rep_determinism = NaN(1,4000);
nov_rep_determinism = NaN(1,4000);
nov_reverse_trace = NaN(1,4000);
rep_reverse_trace = NaN(1,4000);
nov_rep_reverse_trace = NaN(1,4000);
nov_laminar = NaN(1,4000);
rep_laminar = NaN(1,4000);
nov_rep_laminar = NaN(1,4000);
nov_corm = NaN(2,4000);
rep_corm = NaN(2,4000);
nov_rep_corm = NaN(2,4000);
nov_count = 1;
rep_count = 1;
image_complexity = NaN(1,4000);%percent of image containing edges

nov_back_durs1 = [];
rep_back_durs1 = [];
nov_back_ampangle1 = [];
rep_back_ampangle1 = [];

nov_back_sal1 = [];
rep_back_sal1 = [];
all_nov_sal = [];
all_rep_sal = [];


nov_all_fix_durs = NaN(4000,40);
rep_all_fix_durs = NaN(4000,40);
nov_trial_by_trial_sacamps = NaN(4000,40);
rep_trial_by_trial_sacamps = NaN(4000,40);
monk = NaN(1,4000);

for file = 1:length(files);
    load([data_dir files{file}(1:8) '_' files{file}(end) '-fixation.mat'])
    
    if strcmpi(files{file}(1:2),'PW')
        this_monkey = 1;
    elseif strcmpi(files{file}(1:2),'TO')
        this_monkey = 2;
    elseif strcmpi(files{file}(1:2),'RR')
        this_monkey = 3;
    elseif strcmpi(files{file}(1:2),'MF')
        this_monkey = 4;
    end
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
    if strcmpi('PW150416_2',files{file}(1:10))
        viewed(:,1) = 0;
        imgdur(:,1) = 0;
    end
    
    %Let's not even look at images in which the monkey looked away more than the
    %image was displayed for. Also should take care of EOG overflow trials
    viewed(viewed & img_dur > 14000) = 0; %1000 ms buffer for accidental lookaways and cortex timing
    viewed(:,viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    
    setnum = item_set(7:8);
    for novrep = 1:size(viewed,2);
        
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
        offscreen = find(isnan(nov_xy(1,:)));
        gaps = findgaps(offscreen);
        rmv = [];
        for g = 1:size(gaps,1)
            gp = gaps(g,:);
            gp(gp == 0) = [];
            if length(gp) < max_out_time%likely a blink
                rmv = [rmv g];
            end
        end
        gaps(rmv,:) = [];
        if ~isempty(gaps)
            first_offscreen = gaps(1,1);
            if first_offscreen-nov_img_on < 5000
                continue
            else
                nov_img_off = first_offscreen;
            end
        end
        
        nov_fix = round(fixationstats{viewed(1,novrep)}.fixations);
        nov_fix(nov_fix < 1) = 1;
        nov_fix(1,nov_fix(1,:) > imageX) = imageX;
        nov_fix(2,nov_fix(2,:) >= imageY) = imageY-1;
        nov_fixtimes = fixationstats{viewed(1,novrep)}.fixationtimes;
        nov_sac_times = fixationstats{viewed(1,novrep)}.saccadetimes;
        
        pre_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
        post_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
        nov_fix(:,post_img_fix) = [];
        nov_fix(:,pre_img_fix) = [];
        nov_fixtimes(:,post_img_fix) = [];
        nov_fixtimes(:,pre_img_fix ) = [];
        nov_fix_durs = diff(nov_fixtimes,1)+1;
        
        pre_img_sac = find(nov_sac_times(1,:) <= nov_img_on);
        post_img_sac = find(nov_sac_times(2,:) > nov_img_off);
        nov_sac_times(:,post_img_sac) = [];
        nov_sac_times(:,pre_img_sac) = [];
        
        
        %for repeat images
        rep_allval = per(viewed(2,novrep)).allval;
        rep_alltim = per(viewed(2,novrep)).alltim;
        rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        
        rep_xy = fixationstats{viewed(2,novrep)}.XY;
        rep_xy = rep_xy(:,rep_img_on:rep_img_off);
        if (rep_img_off-rep_img_on) > 10500
            continue
        end
        offscreen = find(isnan(rep_xy(1,:)));
        gaps = findgaps(offscreen);
        rmv = [];
        for g = 1:size(gaps,1)
            gp = gaps(g,:);
            gp(gp == 0) = [];
            if length(gp) < max_out_time%likely a blink
                rmv = [rmv g];
            end
        end
        gaps(rmv,:) = [];
        if ~isempty(gaps)
            first_offscreen = gaps(1,1);
            if first_offscreen-rep_img_on < 5000
                continue
            else
                rep_img_off = first_offscreen;
            end
        end
        
        
        rep_fix = round(fixationstats{viewed(2,novrep)}.fixations);
        rep_fix(rep_fix < 1) = 1;
        rep_fix(1,rep_fix(1,:) > imageX) = imageX;
        rep_fix(2,rep_fix(2,:) >= imageY) = imageY-1;
        rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
        rep_sac_times = fixationstats{viewed(2,novrep)}.saccadetimes;
        
        pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        rep_fix(:,post_img_fix) = [];
        rep_fix(:,pre_img_fix) = [];
        rep_fixtimes(:,post_img_fix) = [];
        rep_fixtimes(:,pre_img_fix ) = [];
        rep_fix_durs = diff(rep_fixtimes,1)+1;
        
        pre_img_sac = find(rep_sac_times(1,:) <= rep_img_on);
        post_img_sac = find(rep_sac_times(2,:) > rep_img_off);
        rep_sac_times(:,post_img_sac) = [];
        rep_sac_times(:,pre_img_sac) = [];
        
        %---Grab Important Vairables---%
        %load salience map
        if novrep < 10
            load([image_dir 'LRM' setnum '\S' setnum 'I0' num2str(novrep) '-saliencemap.mat'],'fullmap');
        else
            load([image_dir 'LRM' setnum '\S' setnum 'I' num2str(novrep) '-saliencemap.mat'],'fullmap');
        end
        
        %----All Fixation Durations---%
        nov_all_fix_durs(nov_count,1:size(nov_fixtimes,2)) = diff(nov_fixtimes,1)+1;
        rep_all_fix_durs(rep_count,1:size(rep_fixtimes,2)) = diff(rep_fixtimes,1)+1;
        
        %---Saccade Amplitudes---%
        nov_sac_times = nov_sac_times-nov_img_on;
        for sac = 1:size(nov_sac_times,2)
            sacx =  nov_xy(1,nov_sac_times(1,sac))- nov_xy(1,nov_sac_times(2,sac));
            sacy =  nov_xy(2,nov_sac_times(1,sac))- nov_xy(2,nov_sac_times(2,sac));
            nov_trial_by_trial_sacamps(nov_count,sac) = sqrt(sacx^2+sacy^2);
        end
        
        rep_sac_times = rep_sac_times-rep_img_on;
        for sac = 1:size(rep_sac_times,2)
            sacx =  rep_xy(1,rep_sac_times(1,sac))- rep_xy(1,rep_sac_times(2,sac));
            sacy =  rep_xy(2,rep_sac_times(1,sac))- rep_xy(2,rep_sac_times(2,sac));
            rep_trial_by_trial_sacamps(rep_count,sac) = sqrt(sacx^2+sacy^2);
        end
        
        %---all salience at fixation locations---%
        for f = 1:1:size(nov_fix,2)
            all_nov_sal(nov_count,f) = fullmap(imageY-nov_fix(2,f),nov_fix(1,f));
        end
        for f = 1:size(rep_fix,2)
            all_rep_sal(rep_count,f) = fullmap(imageY-rep_fix(2,f),rep_fix(1,f));
        end
        
        nov_rep_recurence_map = zeros(40,40);
        N1=size(nov_fix,2);
        N2=size(rep_fix,2);
        N = min([N1 N2]);
        N(N > 40) = 40;
        [x,y]=meshgrid(1:N);
        i=find(ones(N)); %forms pairs except for self-pairing
        i=[x(i), y(i)];
        dist =sqrt((nov_fix(1,i(:,1))-rep_fix(1,i(:,2))).^2 +...
            (nov_fix(2,i(:,1))-nov_fix(2,i(:,2))).^2);
        dind = find(dist <= distance_threshold);
        for d = 1:length(dind);
            nov_rep_recurence_map(i(dind(d),1),i(dind(d),2)) = nov_rep_recurence_map(i(dind(d),1),i(dind(d),2))+1;
        end
        R = sum(sum(nov_rep_recurence_map));
        nov_rep_recurence_measure(nov_count) = 100*R/N/N;
        
        determinism_count = 0;
        reverse_trace = 0;
        for i = 1:N-1
            for ii = 2:N-1
                if  nov_rep_recurence_map(ii,i) == 1 && nov_rep_recurence_map(ii+1,i+1) == 1%diagonal
                    determinism_count = determinism_count+1;
                end
                if  nov_rep_recurence_map(ii,i) == 1 && nov_rep_recurence_map(ii-1,i+1) == 1%diagonal
                    reverse_trace = reverse_trace+1;
                end
            end
        end
        nov_rep_determinism(nov_count) = 100*determinism_count/R;
        nov_rep_reverse_trace(nov_count) = 100*reverse_trace/R;
        laminar = nov_rep_recurence_map;
        laminarh = sum(laminar,1);
        laminarh(laminarh <= 1) = 0;
        laminarv = sum(laminar,2);
        laminarv(laminarv <= 1) = 0;
        nov_rep_laminar(nov_count) = 100*sum(laminarv+laminarh')/4/R;
        monk(nov_count) = this_monkey;
        
        if sum(laminarv+laminarh') < 3
            nov_rep_corm(:,nov_count)  = [NaN; NaN];
        else
            [xcm ycm] = centroid(laminar);
            nov_rep_corm(:,nov_count) = [xcm;ycm];
        end
        all_nov_rep_recurence_map = all_nov_rep_recurence_map+nov_rep_recurence_map;
        all_nov_rep_count(1:N) = all_nov_rep_count(1:N)+1;
        
        nov_recurence_map = zeros(40,40);
        N=size(nov_fix,2);
        N(N > 40) = 40;
        [x,y]=meshgrid(1:N);
        i=find(ones(N)); %forms pairs except for self-pairing
        i=[x(i), y(i)];
        dist =sqrt((nov_fix(1,i(:,1))-nov_fix(1,i(:,2))).^2 +...
            (nov_fix(2,i(:,1))-nov_fix(2,i(:,2))).^2);
        dind = find(dist <= distance_threshold);
        for d = 1:length(dind);
            nov_recurence_map(i(dind(d),1),i(dind(d),2)) = nov_recurence_map(i(dind(d),1),i(dind(d),2))+1;
        end
        R = sum(sum(triu(nov_recurence_map)))-N;
        nov_recurence_measure(nov_count) = 100*2*R/N/(N-1);
        
        if sum(nov_recurence_map(id1 == 1)) > 0
            b1 = find(nov_recurence_map(id1 == 1) == 1)+1;
            b0 =  find(nov_recurence_map(id1 == 1) == 1);
            for b = 1:length(b1)
                if nov_recurence_map(b0(b),b1(b)) == 1
                    continue %didn't leave still in same spot
                end
                if (b1(b) < 3) || (b1(b) > length(nov_fix_durs)-2)
                    continue
                end
                %check to make sure didn't look offscreen or blink
                %prior
                prior_sac = find(nov_sac_times(2,:)+1 == nov_fixtimes(1,b1(b)));
                if isempty(prior_sac)
                    continue
                else
                    prior_fix = find(nov_sac_times(1,prior_sac)-1 == nov_fixtimes(2,:));
                    if isempty(prior_fix)
                        continue
                    end
                end
                post_sac = find(nov_sac_times(1,:)-1 == nov_fixtimes(2,b1(b)));
                if isempty(post_sac)
                    continue
                else
                    post_fix = find(nov_sac_times(2,post_sac)+1 == nov_fixtimes(1,:));
                    if isempty(post_fix)
                        continue
                    end
                end
                
                nov_back_durs1 = [nov_back_durs1 ...
                    [b1(b)-2 nov_fix_durs(b1(b)-2)... %2 fixations back
                    b1(b)-1 nov_fix_durs(b1(b)-1)...%1 fixation back i.e. prior fixation
                    b1(b) nov_fix_durs(b1(b)) ...%out fixation
                    b1(b)+1 nov_fix_durs(b1(b)+1)...%return fixation
                    b1(b)+2 nov_fix_durs(b1(b)+2)...%2 fixations forward
                    this_monkey]'];
                
                nov_back_sal1 = [nov_back_sal1 ...
                    [b1(b)-2 fullmap(imageY-nov_fix(2,b1(b)-2),nov_fix(1,b1(b)-2)) ... %2 fixations back
                    b1(b)-1 fullmap(imageY-nov_fix(2,b1(b)-1),nov_fix(1,b1(b)-1))...%1 fixation back i.e. prior fixation
                    b1(b) fullmap(imageY-nov_fix(2,b1(b)),nov_fix(1,b1(b))) ...%out fixation
                    b1(b)+1 fullmap(imageY-nov_fix(2,b1(b)+1),nov_fix(1,b1(b)+1))...%return fixation
                    b1(b)+2 fullmap(imageY-nov_fix(2,b1(b)+2),nov_fix(1,b1(b)+2))...%2 fixations forward
                    this_monkey]'];
                
                prior_amp = sqrt(sum((nov_fix(:,b1(b))-nov_fix(:,prior_fix)).^2));
                post_amp =  sqrt(sum((nov_fix(:,b1(b))-nov_fix(:,post_fix)).^2));
                pre_angle =atan2d((nov_fix(2,b1(b))-nov_fix(2,prior_fix)),(nov_fix(1,b1(b))-nov_fix(1,prior_fix)));
                post_angle = atan2d((nov_fix(2,post_fix)-nov_fix(2,b1(b))),(nov_fix(1,post_fix)-nov_fix(1,b1(b))));
                nov_back_ampangle1 = [nov_back_ampangle1 [post_amp-prior_amp; post_angle-pre_angle]];
            end
        end
        
        
        nov_1backs(nov_count) = sum(nov_recurence_map(id1 == 1))/N;
        nov_2backs(nov_count) = sum(nov_recurence_map(id2 == 1))/N;
        nov_3backs(nov_count) = sum(nov_recurence_map(id3 == 1))/N;
        
        determinism_count = 0;
        reverse_trace = 0;
        for i = 1:N-1
            for ii = 2:N-1
                if ii > i
                    if  nov_recurence_map(ii,i) == 1 && nov_recurence_map(ii+1,i+1) == 1%diagonal
                        determinism_count = determinism_count+1;
                    end
                    if  nov_recurence_map(ii,i) == 1 && nov_recurence_map(ii-1,i+1) == 1%diagonal
                        reverse_trace = reverse_trace+1;
                    end
                    
                end
            end
        end
        nov_determinism(nov_count) = 100*determinism_count/R;
        nov_reverse_trace(nov_count) = 100*reverse_trace/R;
        laminar = nov_recurence_map;
        laminar(tril(laminar) == 1) = 0;
        nov_laminar(nov_count) = 100*(sum(sum(laminar,1) > 1)+  sum(sum(laminar,2) > 1))/2/R;
        monk(nov_count) = this_monkey;
        
        if sum(laminar(:))< 3
            nov_corm(:,nov_count)  = [NaN; NaN];
        else
            [xcm ycm] = centroid(laminar);
            nov_corm(:,nov_count) = [xcm;ycm];
        end
        all_nov_recurence_map = all_nov_recurence_map+nov_recurence_map;
        nov_count = nov_count+1;
        
        
        
        rep_recurence_map = zeros(40,40);
        N=size(rep_fix,2);
        N(N > 40) = 40;
        [x,y]=meshgrid(1:N);
        i=find(ones(N)); %forms pairs except for self-pairing
        i=[x(i), y(i)];
        dist =sqrt((rep_fix(1,i(:,1))-rep_fix(1,i(:,2))).^2 +...
            (rep_fix(2,i(:,1))-rep_fix(2,i(:,2))).^2);
        dind = find(dist <= distance_threshold);
        for d = 1:length(dind);
            rep_recurence_map(i(dind(d),1),i(dind(d),2)) = rep_recurence_map(i(dind(d),1),i(dind(d),2))+1;
        end
        R = sum(sum(triu(rep_recurence_map)))-N;
        rep_recurence_measure(rep_count) = 100*2*R/N/(N-1);
        
        if sum(rep_recurence_map(id1 == 1)) > 0
            b1 = find(rep_recurence_map(id1 == 1) == 1)+1;
            b0 =  find(rep_recurence_map(id1 == 1) == 1);
            for b = 1:length(b1)
                if rep_recurence_map(b0(b),b1(b)) == 1
                    continue %didn't leave still in same spot
                end
                if (b1(b) < 3) || (b1(b) > length(rep_fix_durs)-2)
                    continue
                end
                %check to make sure didn't look offscreen or blink
                %prior
                prior_sac = find(rep_sac_times(2,:)+1 == rep_fixtimes(1,b1(b)));
                if isempty(prior_sac)
                    continue
                else
                    prior_fix = find(rep_sac_times(1,prior_sac)-1 == rep_fixtimes(2,:));
                    if isempty(prior_fix)
                        continue
                    end
                end
                post_sac = find(rep_sac_times(1,:)-1 == rep_fixtimes(2,b1(b)));
                if isempty(post_sac)
                    continue
                else
                    post_fix = find(rep_sac_times(2,post_sac)+1 == rep_fixtimes(1,:));
                    if isempty(post_fix)
                        continue
                    end
                end
                
                rep_back_durs1 = [rep_back_durs1 ...
                    [b1(b)-2 rep_fix_durs(b1(b)-2)... %2 fixations back
                    b1(b)-1 rep_fix_durs(b1(b)-1)...%1 fixation back i.e. prior fixation
                    b1(b) rep_fix_durs(b1(b)) ...%out fixation
                    b1(b)+1 rep_fix_durs(b1(b)+1)...%return fixation
                    b1(b)+2 rep_fix_durs(b1(b)+2)...%2 fixations forward
                    this_monkey]'];
                
                rep_back_sal1 = [rep_back_sal1 ...
                    [b1(b)-2 fullmap(imageY-rep_fix(2,b1(b)-2),rep_fix(1,b1(b)-2)) ... %2 fixations back
                    b1(b)-1 fullmap(imageY-rep_fix(2,b1(b)-1),rep_fix(1,b1(b)-1))...%1 fixation back i.e. prior fixation
                    b1(b) fullmap(imageY-rep_fix(2,b1(b)),rep_fix(1,b1(b))) ...%out fixation
                    b1(b)+1 fullmap(imageY-rep_fix(2,b1(b)+1),rep_fix(1,b1(b)+1))...%return fixation
                    b1(b)+2 fullmap(imageY-rep_fix(2,b1(b)+2),rep_fix(1,b1(b)+2))...%2 fixations forward
                    this_monkey]'];
                
                
                prior_amp = sqrt(sum((rep_fix(:,b1(b))-rep_fix(:,prior_fix)).^2));
                post_amp =  sqrt(sum((rep_fix(:,b1(b))-rep_fix(:,post_fix)).^2));
                pre_angle =atan2d((rep_fix(2,b1(b))-rep_fix(2,prior_fix)),(rep_fix(1,b1(b))-rep_fix(1,prior_fix)));
                post_angle = atan2d((rep_fix(2,post_fix)-rep_fix(2,b1(b))),(rep_fix(1,post_fix)-rep_fix(1,b1(b))));
                rep_back_ampangle1 = [rep_back_ampangle1 [post_amp-prior_amp; post_angle-pre_angle]];
            end
        end
        
        rep_1backs(rep_count) = sum(rep_recurence_map(id1 == 1))/N;
        rep_2backs(rep_count) = sum(rep_recurence_map(id2 == 1))/N;
        rep_3backs(rep_count) = sum(rep_recurence_map(id3 == 1))/N;
        
        determinism_count = 0;
        reverse_trace = 0;
        for i = 1:N-1
            for ii = 2:N-1
                if ii > i
                    if  rep_recurence_map(ii,i) == 1 && rep_recurence_map(ii+1,i+1) == 1%diagonal
                        determinism_count = determinism_count+1;
                    end
                    if  rep_recurence_map(ii,i) == 1 && rep_recurence_map(ii-1,i+1) == 1%diagonal
                        reverse_trace = reverse_trace+1;
                    end
                    
                end
            end
        end
        rep_determinism(rep_count) = 100*determinism_count/R;
        rep_reverse_trace(rep_count) = 100*reverse_trace/R;
        
        laminar = rep_recurence_map;
        laminar(tril(laminar) == 1) = 0;
        rep_laminar(rep_count) = 100*(sum(sum(laminar,1) > 1)+  sum(sum(laminar,2) > 1))/2/R;
        
        if sum(laminar(:))< 3
            rep_corm(:,rep_count)  = [NaN; NaN];
        else
            [xcm ycm] = centroid(laminar);
            rep_corm(:,rep_count) = [xcm;ycm];
        end
        all_rep_recurence_map = all_rep_recurence_map+rep_recurence_map;
        rep_count = rep_count+1;
        
        
        nov_recurence_map = nov_recurence_map + nov_recurence_map;
        rep_recurence_map = rep_recurence_map + rep_recurence_map;
    end
end
%% 1-backs/2-backs/3-backs
figure
subplot(1,2,1)
errorbar([nanmean(nov_1backs) nanmean(nov_2backs)  nanmean(nov_3backs)],[nanstd(nov_1backs) nanstd(nov_2backs)  nanstd(nov_3backs)]./sqrt(69))
hold on
errorbar([nanmean(rep_1backs) nanmean(rep_2backs)  nanmean(rep_3backs)],[nanstd(rep_1backs) nanstd(rep_2backs)  nanstd(rep_3backs)]./sqrt(69))
hold off
set(gca,'Xtick',[1 2 3])
box off
xlabel('# of Backs')
ylabel('Proportion of Fixations')
box off

local_nov = nov_1backs+nov_2backs +nov_3backs;
local_rep = rep_1backs+rep_2backs+rep_3backs;

subplot(1,2,2)
bar([nanmean(local_nov) nanmean(local_rep)])
hold on
errorb([nanmean(local_nov) nanmean(local_rep)],...
    [nanstd(local_nov) nanstd(local_rep)]./sqrt(69))
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Novel','Repeat'})
box off


%%
figure
%---For Novel---%

rm_nov = all_nov_recurence_map;
for r = 1:length(rm_nov);
    fix_count = rm_nov(r,r);
    for i = 1:r-1
        rm_nov(r,i) = rm_nov(r,i)/fix_count;
        rm_nov(i,r) = rm_nov(i,r)/fix_count;
    end
    rm_nov(r,r) = 1;
end
rm_nov(find(eye(size(rm_nov)))) = NaN;
rm_nov2 = rm_nov;
rm_nov2(id0 == 1) = NaN;
rm_nov2(id0' == 1) = NaN;
rm_nov2(id1 == 1) = NaN;
rm_nov2(id1' == 1) = NaN;
% rm_nov2(id2 == 1) = NaN;
% rm_nov2(id2' == 1) = NaN;
rm_nov2 = 100*rm_nov2(1:25,1:25);

rm_rep = all_rep_recurence_map;
for r = 1:length(rm_rep);
    fix_count = rm_rep(r,r);
    for i = 1:r-1
        rm_rep(r,i) = rm_rep(r,i)/fix_count;
        rm_rep(i,r) = rm_rep(i,r)/fix_count;
    end
    rm_rep(r,r) = 1;
end
rm_rep(find(eye(size(rm_rep)))) = NaN;

rm_rep2 = rm_rep;
rm_rep2(id0 == 1) = NaN;
rm_rep2(id0' == 1) = NaN;
rm_rep2(id1 == 1) = NaN;
rm_rep2(id1' == 1) = NaN;
% rm_rep2(id2 == 1) = NaN;
% rm_rep2(id2' == 1) = NaN;
rm_rep2 = 100*rm_rep2(1:25,1:25);


rm_change = 100*(rm_rep-rm_nov);
rm_change(find(eye(size(rm_change)))) = NaN;
rm_change = rm_change(1:25,1:25);

rm_nov3 = rm_nov(1:25,1:25);
upperbound = prctile(rm_nov3(:),97.5);
lowerbound = prctile(rm_nov3(:),2.5);

subplot(2,2,1)
h = imagesc(100*rm_nov3);
axis xy, set(h,'alphadata',~isnan(rm_nov(1:25,1:25)));
caxis([100*lowerbound 100*upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
box off
colorbar
colormap('viridis')
title('Novel')
axis square
box off



upperbound = prctile(rm_nov2(:),97.5);
lowerbound = prctile(rm_nov2(:),2.5);

subplot(2,2,2)
h = imagesc(rm_nov2);
axis xy, set(h,'alphadata',~isnan(rm_nov2));
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
box off
colorbar
colormap('viridis')
title('Novel')
axis square
box off

%---For Repeat---%
upperbound = prctile(rm_rep2(:),97.5);
lowerbound = prctile(rm_rep2(:),2.5);

subplot(2,2,3)
h = imagesc(rm_rep2);
axis xy, set(h,'alphadata',~isnan(rm_rep2));
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
colormap('viridis')
colorbar
box off
title('Repeat')
axis square
box off

%---For Repeat minus Novel---%
upperbound = prctile(rm_change(:),97.5);
lowerbound = prctile(rm_change(:),2.5);

subplot(2,2,4)

f = fspecial('average',3);

rm2c = rm_change;
rm2c(1,1) = rm2c(1,2);
for i = 2:25
    rm2c(i,i) = rm2c(i,i-1);
end

h = imagesc(imfilter(rm2c,f));
axis xy, set(h,'alphadata',~isnan(rm_change));
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
colormap('viridis')
box off
colorbar
title('Repeat-Novel')
box off
axis square

%%
image_count = length(nov_recurence_measure);

figure
subplot(3,3,1)
hold on
[x,xi] = hist(nov_recurence_measure,0:0.25:100);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/image_count));
[x,xi] = hist(rep_recurence_measure,0:0.25:100);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/image_count));
[x,xi] = hist(nov_rep_recurence_measure,0:0.25:100);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/image_count));
hold off
xlabel('Recurrence (%)')
ylabel('Proportion of Images')
legend('Novel','Repeat','Novel/Repeat')
title('Recurrence Rate')
xlim([0 20])
ylim([0 .15])

subplot(3,3,2)
hold on
[x,xi] = hist(nov_determinism,0:0.25:100);
x(1) = max(x(2:end))/2;
plot(xi,filtfilt(ones(1,8)*1/8,1,x/image_count));
[x,xi] = hist(rep_determinism,0:0.25:100);
x(1) = max(x(2:end))/2;
plot(xi,filtfilt(ones(1,8)*1/8,1,x/image_count));
[x,xi] = hist(nov_rep_determinism,0:0.25:100);
x(1) = max(x(2:end))/2;
plot(xi,filtfilt(ones(1,8)*1/8,1,x/image_count));
hold off
xlabel('Determinism (%)')
ylabel('Proportion of Images')
title('Retrace Rate')
xlim([0 40])

subplot(3,3,3)
hold on
[x,xi] = hist(nov_reverse_trace,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/image_count));
[x,xi] = hist(rep_reverse_trace,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/image_count));
[x,xi] = hist(nov_rep_reverse_trace,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/image_count));
hold off
xlabel('Determinism (%)')
ylabel('Proportion of Images')
title('Reverse Retrace Rate')
xlim([0 60])

subplot(3,3,4)
hold on
[x,xi] = hist(nov_laminar,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/image_count));
[x,xi] = hist(rep_laminar,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/image_count));
[x,xi] = hist(nov_rep_laminar,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/image_count));
hold off
xlabel('Laminarity (%)')
ylabel('Proportion of Images')
title('Laminar Re-Repeated Fixations')
xlim([0 50])

subplot(3,3,5)
hold on
[x,xi] = hist(nov_corm(1,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/image_count));
[x,xi] = hist(rep_corm(1,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/image_count));
[x,xi] = hist(nov_rep_corm(1,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/image_count));
hold off
xlabel('Ordinal Fixation #')
ylabel('Proportion of Images')
title('Horizontal Center of Mass')
xlim([0 30])

center_nov_count = sum(~isnan(nov_corm(1,:)));
center_rep_count = sum(~isnan(rep_corm(1,:)));
center_nov_rep_count = sum(~isnan(nov_rep_corm(1,:)));


subplot(3,3,6)
hold on
[x,xi] = hist(nov_corm(2,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/center_nov_count));
[x,xi] = hist(rep_corm(2,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/center_rep_count));
[x,xi] = hist(nov_rep_corm(2,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/center_nov_rep_count));
hold off
xlabel('Ordinal Fixation #')
ylabel('Proportion of Images')
title('Vertical Center of Mass')
xlim([0 30])


subplot(3,3,7)
hold on
[x,xi] = hist(abs(nov_corm(1,:)-nov_corm(2,:)),0:0.25:20);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/center_nov_count));
[x,xi] = hist(abs(rep_corm(1,:)-rep_corm(2,:)),0:0.25:20);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/center_rep_count));
[x,xi] = hist(abs(nov_rep_corm(1,:)-nov_rep_corm(2,:)),0:0.25:20);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/center_nov_rep_count));
hold off
xlabel('Local --> Global')
ylabel('Proportion of Images')
title('Distance from Center: Local vs Global')


nov = abs(nov_corm(1,:)-nov_corm(2,:));
rep = abs(rep_corm(1,:)-rep_corm(2,:));

subplot(3,3,8)
hold on
[x,xi] = hist(rep-nov,-10:0.25:10);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/sum(~isnan(nov-rep))));
hold off
xlabel('Change in Local --> Global')
ylabel('Proportion of Images')
title('Change from Global to Local')


subtitle(['All Recurrence Measurse for ' files{1}(1:2)])
%%
nov_fix_dur_avg = NaN(max(monk),40);
rep_fix_dur_avg = NaN(max(monk),40);
nov_sal_avg = NaN(max(monk),40);
rep_sal_avg = NaN(max(monk),40);
image_count = NaN(1,max(monk));
median_num_fix = NaN(2,max(monk));
for m = 1:max(monk)
    these_images = find(monk == m);
    image_count(m) = length(these_images);
    
    nov_durs = nov_all_fix_durs(these_images,:);
    nov_fix_count = sum(~isnan(nov_durs'));
    median_num_fix(1,m) = median(nov_fix_count);
    nov_fix_dur_avg(m,1:median_num_fix(1,m)) = nanmean(nov_durs(:,1:median_num_fix(1,m)));
    rep_durs = rep_all_fix_durs(these_images,:);
    rep_fix_count = sum(~isnan(rep_durs'));
    median_num_fix(2,m) = median(rep_fix_count);
    rep_fix_dur_avg(m,1:median_num_fix(2,m)) = nanmean(rep_durs(:,1:median_num_fix(2,m)));
    
    nov_sal = all_nov_sal(these_images,:);
    nov_sal_avg(m,1:median_num_fix(1,m)) = nanmean(nov_sal(:,1:median_num_fix(1,m)));
    
    rep_sal = all_nov_sal(these_images,:);
    rep_sal_avg(m,1:median_num_fix(2,m)) = nanmean(rep_sal(:,1:median_num_fix(2,m)));
end

%%
nov_residual_duration = NaN(5,length(nov_back_durs1));
nov_residual_sal = NaN(5,length(nov_back_durs1));
for n = 1:length(nov_back_durs1)
    this_monkey = nov_back_durs1(end,n);
    if nov_back_durs1(9,n) <= median_num_fix(1,this_monkey)
        for f = 1:5
            nov_residual_duration(f,n) = nov_back_durs1(f*2,n)-nov_fix_dur_avg(this_monkey,f*2-1);
            nov_residual_sal(f,n) = nov_back_sal1(f*2,n)-nov_sal_avg(this_monkey,f*2-1);
        end
    else
        continue
    end
end

nov_residual_duration_best = NaN(5,length(nov_back_durs1));
nov_residual_sal_best =  NaN(5,length(nov_back_durs1));
for n = 1:length(nov_back_durs1)
    this_monkey = nov_back_durs1(end,n);
    if abs(nov_back_ampangle1(1,n)) > 1
        continue
    end
    if abs(nov_back_ampangle1(2,n)) < 160 || abs(nov_back_ampangle1(2,n)) > 200
        continue
    end
    if nov_back_durs1(9,n) <= median_num_fix(1,this_monkey)
        for f = 1:5
            nov_residual_duration_best(f,n) = nov_back_durs1(f*2,n)-nov_fix_dur_avg(this_monkey,f*2-1);
            nov_residual_sal_best(f,n) = nov_back_sal1(f*2,n)-nov_sal_avg(this_monkey,f*2-1);
        end
    else
        continue
    end
end

rep_residual_duration = NaN(5,length(rep_back_durs1));
rep_residual_sal = NaN(5,length(rep_back_durs1));
for n = 1:length(rep_back_durs1)
    this_monkey = rep_back_durs1(end,n);
    if rep_back_durs1(9,n) <= median_num_fix(1,this_monkey)
        for f = 1:5
            rep_residual_duration(f,n) = rep_back_durs1(f*2,n)-nanmean(rep_fix_dur_avg(this_monkey,f*2-1));
            rep_residual_sal(f,n) = rep_back_sal1(f*2,n)-rep_sal_avg(this_monkey,f*2-1);
        end
    else
        continue
    end
end

rep_residual_duration_best = NaN(5,length(rep_back_durs1));
rep_residual_sal_best = NaN(5,length(rep_back_durs1));
for n = 1:length(rep_back_durs1)
    this_monkey = rep_back_durs1(end,n);
    if abs(rep_back_ampangle1(1,n)) > 1
        continue
    end
    if abs(rep_back_ampangle1(2,n)) < 160 || abs(rep_back_ampangle1(2,n)) > 200
        continue
    end
    if rep_back_durs1(9,n) <= median_num_fix(1,this_monkey)
        for f = 1:5
            rep_residual_duration_best(f,n) = rep_back_durs1(f*2,n)-nanmean(rep_fix_dur_avg(this_monkey,f*2-1));
            rep_residual_sal_best(f,n) = rep_back_sal1(f*2,n)-rep_sal_avg(this_monkey,f*2-1);
        end
    else
        continue
    end
end


figure
subplot(2,3,1)
plot(nanmean(nov_residual_duration'))
hold on
plot(nanmean(rep_residual_duration'))
plot(nanmean(nov_residual_duration_best'))
plot(nanmean(rep_residual_duration_best'))
hold off
legend('Novel','Repeat','Best Novel','Best Repeat')
box off
set(gca,'Xtick',1:5)
set(gca,'XtickLabel',{'Prior-1','Prior','Out','Return','Return+1'})
ylabel(sprintf('Residual (Observed-Expected) \n Fixation Duration (ms)'))
title('Facilitation of Return(?)')

subplot(2,3,2)
hist([abs(nov_back_ampangle1(2,:)) abs(rep_back_ampangle1(2,:))] ,25)
hold on
plot([160 160],[0 4500],'k--')
plot([200 200],[0 4500],'k--')
hold off
xlabel('\Delta Saccade Angle')
ylabel('Fixation Count')
xlim([145 225])
box off
title('Distribution of \Delta Saccade Angles')

subplot(2,3,3)
hist([nov_back_ampangle1(1,:)/24 rep_back_ampangle1(1,:)/24] ,25)
hold on
plot([-1 -1],[0 1400],'k--')
plot([1 1],[0 1400],'k--')
hold off
xlabel('\Delta Saccade Amplitude (dva)')
ylabel('Fixation Count')
box off
title('Distribution of \Delta Saccade Amplitudes')
xlim([-2 2])

dirs = abs(nov_back_ampangle1(2,:));
subplot(2,3,4)
plot(abs(dirs-180),nov_residual_duration(3,:),'.k')
ylim([-150 150])
xlim([0 45])
xlabel('||\Delta|| Saccade Angle')
ylabel('Residual Fixation Duration (ms)')
box off
title('Relationship of \Delta Angle and Out Residual')


subplot(2,3,5)
plot(abs(nov_back_ampangle1(1,:)/24),nov_residual_duration(4,:),'.k')
ylim([-150 250])
xlim([0 2])
xlabel('||\Delta|| Saccade Amplitude')
ylabel('Residual Fixation Duration (ms)')
box off
title('Relationship of \Delta Amplitude and Return Residual')

subplot(2,3,6)
plot(abs(dirs-180),nov_residual_duration(4,:),'.k')
ylim([-150 250])
xlim([0 45])
xlabel('||\Delta|| Saccade Angle')
ylabel('Residual Fixation Duration (ms)')
box off
title('Relationship of \Delta Amplitude and Return Residual')

%%
figure
hold on
plot(nanmean(nov_residual_sal'))
plot(nanmean(rep_residual_sal'))
plot(nanmean(nov_residual_sal_best'))
plot(nanmean(rep_residual_sal_best'))
hold off
legend('Novel','Repeat','Best Novel','Best Repeat')
box off
set(gca,'Xtick',1:5)
set(gca,'XtickLabel',{'Prior-1','Prior','Out','Return','Return+1'})
ylabel(sprintf('Residual (Observed-Expected) \n Fixation Duration (ms)'))
title('Salience Around Return Fixation Locations')
%%
rm_nov_rep = all_nov_rep_recurence_map;
for r = 1:length(rm_nov_rep);
    for i = 1:r-1
        rm_nov_rep(r,i) = rm_nov_rep(r,i)/all_nov_rep_count(i);
        rm_nov_rep(i,r) = rm_nov_rep(i,r)/all_nov_rep_count(i);
    end
    rm_nov_rep(i,i) = rm_nov_rep(i,i)/all_nov_rep_count(i);
end
%%

rm2 = rm_nov_rep(1:25,1:25);
figure
subplot(2,2,1)
imagesc(rm2)
axis xy
colormap('viridis')
colorbar
axis square
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')

rm2 = rm_nov_rep;
rm2(eye(40) == 1) = NaN;
rm2 = rm2(1:25,1:25);
subplot(2,2,2)
axis xy
h = imagesc(rm2);
axis xy
colormap('viridis')
colorbar
set(h,'alphadata',~isnan(rm2));
axis square
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')

rm2 = rm_nov_rep;
rm2(eye(40) == 1) = NaN;
rm2(id0 == 1) = NaN;
rm2(id0' == 1) = NaN;
rm2 = rm2(1:25,1:25);
subplot(2,2,3)
axis xy
h = imagesc(rm2);
axis xy
colormap('viridis')
colorbar
set(h,'alphadata',~isnan(rm2));
axis square
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')

rm2 = rm_nov_rep;
rm2(eye(40) == 1) = NaN;
rm2(id0 == 1) = NaN;
rm2(id0' == 1) = NaN;
rm2(id1 == 1) = NaN;
rm2(id1' == 1) = NaN;
rm2 = rm2(1:25,1:25);
subplot(2,2,4)
axis xy
h = imagesc(rm2);
axis xy
colormap('viridis')
colorbar
set(h,'alphadata',~isnan(rm2));
axis square
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')

subtitle('Novel vs Repeat: Scan Path Similarity')
%%
%---Plot Fixation Durations By Fixation Number---%
nov_mean_fix_dur = nov_all_fix_durs;
rep_mean_fix_dur =rep_all_fix_durs;
nov_num_fix = NaN(1,4);
rep_num_fix = NaN(1,4);
for monks = 1:4
    these_sets = find(monk == monks);
    first_dur = nanmean(nov_mean_fix_dur(these_sets,1)); 
    nov_mean_fix_dur(these_sets,:) = nov_mean_fix_dur(these_sets,:)./first_dur;
    nov_num_fix(monks) = median(sum(~isnan(nov_mean_fix_dur(these_sets,:)),2));
    rep_mean_fix_dur(these_sets,:) = rep_mean_fix_dur(these_sets,:)./first_dur;
    rep_num_fix(monks) = median(sum(~isnan(rep_mean_fix_dur(these_sets,:)),2));
end
 
nov_num_fix = min(nov_num_fix);
rep_num_fix = min(rep_num_fix);

figure
hold on
errorbar(nanmean(nov_mean_fix_dur(:,1:nov_num_fix)),nanstd(nov_mean_fix_dur(:,1:nov_num_fix))...
    ./sqrt(sum(~isnan(nov_mean_fix_dur(:,1:nov_num_fix)))),'b')
errorbar(nanmean(rep_mean_fix_dur(:,1:rep_num_fix)),nanstd(rep_mean_fix_dur(:,1:rep_num_fix))...
    ./sqrt(sum(~isnan(rep_mean_fix_dur(:,1:rep_num_fix)))),'r')
hold off
xlabel('Ordinal Fixation #')
ylabel('Normalized Fixation Duration')
legend('Novel','Repeat')
ylim([0.85 1.25])
xlim([0 27.5])
%%
%---Plot Saccade Durations By Saccade Number---%
nov_mean_sac_amp = nov_trial_by_trial_sacamps;
rep_mean_sac_amp =rep_trial_by_trial_sacamps;
nov_num_fix = NaN(1,4);
rep_num_fix = NaN(1,4);
for monks = 1:4
    these_sets = find(monk == monks);
    first_dur = nanmean(nov_mean_sac_amp(these_sets,1)); 
    nov_mean_sac_amp(these_sets,:) = nov_mean_sac_amp(these_sets,:)./first_dur;
    nov_num_fix(monks) = median(sum(~isnan(nov_mean_sac_amp(these_sets,:)),2));
    rep_mean_sac_amp(these_sets,:) = rep_mean_sac_amp(these_sets,:)./first_dur;
    rep_num_fix(monks) = median(sum(~isnan(rep_mean_sac_amp(these_sets,:)),2));
end
 
nov_num_fix = min(nov_num_fix);
rep_num_fix = min(rep_num_fix);

figure
hold on
errorbar(nanmean(nov_mean_sac_amp(:,1:nov_num_fix)),nanstd(nov_mean_sac_amp(:,1:nov_num_fix))...
    ./sqrt(sum(~isnan(nov_mean_sac_amp(:,1:nov_num_fix)))),'b')
errorbar(nanmean(rep_mean_sac_amp(:,1:rep_num_fix)),nanstd(rep_mean_sac_amp(:,1:rep_num_fix))...
    ./sqrt(sum(~isnan(rep_mean_sac_amp(:,1:rep_num_fix)))),'r')
hold off
xlabel('Ordinal Saccade #')
ylabel('Normalized Saccade Amplitude')
legend('Novel','Repeat')
ylim([0.95 1.6])
xlim([0 28.5])