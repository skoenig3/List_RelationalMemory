%written 10/26/2017 by Seth Konig
%I'm worried that the saccade amplitudes from Tobii are odd becasue of some
%import function issues

clar

pre_files = {'TO150513.2','TO150514.2','TO150515.2','TO150518.2','TO150519.2',...
    'TO150520.2','TO150521.2','TO150522.2','TO150526.2','TO150527.2',...
    'TO150528.2','TO150529.2','TO150601.2','TO150602.2','TO150603.2',...
    'TO150604.2','TO150605.2'};

pre_item_files = {'ListRM01.itm','ListRM02.itm','ListRM03.itm','ListRM04.itm','ListRM05.itm',...
    'ListRM06.itm','ListRM07.itm','ListRM08.itm','ListRM09.itm','ListRM10.itm',...
    'ListRM11.itm','ListRM12.itm','ListRM13.itm','ListRM14.itm','ListRM15.itm',...
    'ListRM16.itm','ListRM17.itm'};

post_files = {'TO170404.2','TO170405.2','TO170406.2','TO170407.2',...
    'TO170410.2','TO170411.2','TO170412.2','TO170414.2',...
    'TO170419.2','TO170421.2','TO170425.2','TO170426.2',...
    'TO170427.2','TO170428.2','TO170501.2'};

post_item_files = {'ListRM22.itm','ListRM23.itm','ListRM24.itm','ListRM25.itm',...
    'ListRM26.itm','ListRM27.itm','ListRM28.itm','ListRM30.itm',...
    'ListRM33.itm','ListRM35.itm','ListRM36.itm','ListRM37.itm',...
    'ListRM38.itm','ListRM39.itm','ListRM40.itm'};


imageX = 800;
imageY = 600;
samprate = 5;

pre_nov_sac_amps = NaN(length(pre_files),50);
all_pre_nov_sac_amps = [];
pre_sac_nov_amps_volt = NaN(length(pre_files),50);
all_pre_nov_sac_amps_volt = [];
pre_rep_sac_amps = NaN(length(pre_files),50);
all_pre_rep_sac_amps = [];
pre_sac_rep_amps_volt = NaN(length(pre_files),50);
all_pre_rep_sac_amps_volt = [];

for file = 1:length(pre_files)
    disp(['Analyzing ' pre_files{file}])
    item_set = pre_item_files{file};
    
    %---Import Cortex data file---%
    cortexfile = ['R:\Cortex Data\Tobii\' pre_files{file}];
    
    cnd_file = [item_set(1:end-4) '.cnd'];
    
    [time_arr,event_arr,eog_arr,epp_arr,~,~]  = get_ALLdata(cortexfile);
    [itmlist,clrchng_locations,first_img_item,imgs] = read_ListRM_itm_and_cnd_files(item_set,cnd_file);
    
    %---Get the calibration for the eye data---%
    %essentially the same as all other task's we use
    numrpt = size(event_arr,2);
    valrptcnt = 0;
    clear per clrchgind
    for rptlop = 1:numrpt
        if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) < first_img_item
            if size(find(event_arr(:,rptlop) == 3)) ~=0 %if clrchng trial was rewarded
                perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
                perendind = find(event_arr(:,rptlop) == 24);
                if length( perbegind) > 1
                    perbegind = perbegind(2);
                    perendind = perendind(2);
                end
                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                begtimdum = time_arr(perbegind,rptlop)-100;
                endtimdum = time_arr(perendind,rptlop);
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    clrchgind(valrptcnt)=rptlop;
                    per(valrptcnt).begsmpind = begtimdum;
                    per(valrptcnt).endsmpind = endtimdum;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                    per(valrptcnt).event = rptlop;
                end
            end
        end
    end
    
    %Don't keep first 6 successful trials. These trials are all central
    %cailbration trials for calibration offset correct
    per(1:5) = [];
    
    clear cnd
    numrpt = size(per,2);
    cnd = zeros(1,numrpt);
    for rptlop = 1:numrpt
        cnd(rptlop)=per(rptlop).cnd-1000;
    end
    
    % Create structures x and y of the corresponding average eye data for each trial
    % instance (l) of each condition (k)
    spacex = [-12,-6,0,6,12];%what actually gets displayed
    spacey = [-8,-4,0,4,8];%what actually gets displayed
    x = cell(length(spacex),length(spacey));%---For Calibration with Eye tracking data with cp2tform---%
    y = cell(length(spacex),length(spacey));
    control = NaN(length(cnd),2);
    for k = 1:length(cnd)
        control(k,:) = clrchng_locations{cnd(k)}';
        
        xi = find(clrchng_locations{cnd(k)}(1) == spacex);
        yi = find(clrchng_locations{cnd(k)}(2) == spacey);
        eyeind = floor(((per(k).begsmpind-1000)/samprate)*2):(floor((per(k).endsmpind-1000)/samprate))*2;
        evenind = eyeind(logical(~rem(eyeind,2)));
        oddind =  eyeind(logical(rem(eyeind,2)));
        x{xi,yi} = [x{xi,yi} mean(eog_arr(oddind,per(k).event))];
        y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))];
    end
    
    %Test for errors%
    count = zeros(length(spacey),length(spacex));
    for xi = 1:length(spacex);
        for yi = 1:length(spacey);
            count(yi,xi) = sum(control(:,1) == spacex(xi) & control(:,2) == spacey(yi));
        end
    end
    if any(count < 10);
        disp('Calibration trial analysis incomplete or error')
        disp('Check number of calibration pionts or task not finished')
    end
    
    clear meanx meany
    for k=1:numel(x)
        xss = x{k};
        low = mean(xss)-std(xss);
        high = mean(xss)+std(xss);
        xss(xss < low) = [];
        xss(xss > high) = [];
        meanx(k)=median(xss);
    end
    for k=1:numel(y)
        yss = y{k};
        low = mean(yss)-std(yss);
        high = mean(yss)+std(yss);
        yss(yss < low) = [];
        yss(yss > high) = [];
        meany(k)=median(y{k});
    end
    
    controlx = [];
    controly = [];
    for i = 1:length(spacex);
        for ii = 1:length(spacey);
            controly = [controly spacey(i)];
            controlx = [controlx spacex(ii)];
        end
    end
    
    %want to change this to MSE estimate with different calibration functions,
    %probably affine and polynomial 3 or 4
    tform = cp2tform([controlx' controly'], [meanx' meany'],'affine');
    tform.forward_fcn = tform.inverse_fcn;
    
    newx = [];
    newy = [];
    figure(101)
    hold on
    for i = 1:length(controlx);
        [x,y] = tformfwd(tform,meanx(i),meany(i));
        plot(x,y,'*b')
        newx(i) = x;
        newy(i) = y;
    end
    hold off
    xlim([-17.5 17.5])
    ylim([-12.5 12.5])
    
    
    figure(102)
    hold on
    plot(meanx,meany,'*b')
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Import Image Viewing Data---%
    
    new_eog_arr = [];
    numrpt = size(event_arr,2);
    new_eog_arr=[];
    valrptcnt = 0;
    clear per clrchgind
    for rptlop = 1:numrpt
        if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) >= first_img_item
            if size(find(event_arr(:,rptlop) == 23)) ~=0 %if image actually turned on
                perbegind = find(event_arr(:,rptlop) == 100,1,'first'); %eye data on
                perendind = find(event_arr(:,rptlop) == 101,1,'first'); %eye data off
                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                begtimdum = time_arr(perbegind,rptlop);
                endtimdum = time_arr(perendind,rptlop);
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    vpcind(valrptcnt)=rptlop;
                    per(valrptcnt).begsmpind = begtimdum;
                    per(valrptcnt).endsmpind = endtimdum;
                    per(valrptcnt).begpos = 1;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                    new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                end
            end
        end
    end
    
    %---get eye data for only when fixation cross or picture is displayed---%
    eyedat = cell(1,length(per));
    for trlop=1:size(per,2)
        trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
        horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
        vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
        
        picstart=1*samprate;
        picend=per(trlop).endsmpind-per(trlop).begsmpind;%added in case sometimes get weird
        %indexing artifacts that are off by 1 index due to cortex having a
        %clock speed with a 1 ms resoultion and the eye data collected at a 5 ms
        %resoultion
        
        eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate));
        eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate));
    end
    
    %---Recalibrate and automatically scale eye data---%
    fixationstats = cell(1,length(eyedat));
    raw_eyedat = cell(1,length(eyedat));
    for eye = 1:length(eyedat)
        x = eyedat{eye}(1,:);
        y = eyedat{eye}(2,:);
        x_raw = x;
        y_raw = y;
        [x,y] = tformfwd(tform,x,y);
        
        x = 24*x; %convert from cortex dva to pixels
        y = 24*y; %convert from cortex dva to pixels
        x = x+imageX/2;
        y = y+imageY/2;
        
        %---Remove Data outside of image---%
        y_raw(x < -24) = NaN;
        x_raw(x < -24) = NaN;
        x_raw(y < -24) = NaN;
        y_raw(y < -24) = NaN;
        y_raw(x > imageX+24) = NaN;
        x_raw(x > imageX+24)= NaN;
        x_raw(y > imageY+24) = NaN;
        y_raw(y > imageY+24) = NaN;
        
        y(x < -24) = NaN;
        x(x < -24) = NaN;
        y(x > imageX+24) = NaN;
        x(x > imageX+24)= NaN;
        x(y < -24) = NaN;
        y(y < -24) = NaN;
        x(y > imageY+24) = NaN;
        y(y > imageY+24) = NaN;
        
        raw_eyedat{eye} = [x_raw; y_raw];
        eyedat{eye} = [x;y];
        
        [fixationtimes,saccadetimes] = adaptive_VT_Recording_Data2(eyedat{eye});
        
        fixationstats{eye}.fixationtimes = fixationtimes;
        fixationstats{eye}.saccadetimes = saccadetimes;
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
    imgnum = 1:96;
    
    %Let's not even look at images in which the monkey looked away more than the
    %image was displayed for. Also should take care of EOG overflow trials
    viewed(viewed & img_dur > 12000) = 0; %1000 ms buffer for accidental lookaways and cortex timing
    imgnum(viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    viewed(:,viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    
    nov_trial_by_trial_sacamps = NaN(96,50); %each trials saccade amplitudes
    rep_trial_by_trial_sacamps = NaN(96,50); %each trials fixation amplitudes
    nov_trial_by_trial_sacamps_volts = NaN(96,50); %each trials saccade amplitudes
    rep_trial_by_trial_sacamps_volts = NaN(96,50); %each trials fixation amplitudes
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
        
        nov_sactimes = fixationstats{viewed(1,novrep)}.saccadetimes;
        prep_img_saccades = find(nov_sactimes(1,:) <= nov_img_on);
        postp_img_saccades = find(nov_sactimes(2,:) > nov_img_off);
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
        
        rep_sactimes = fixationstats{viewed(2,novrep)}.saccadetimes;
        prep_img_saccades = find(rep_sactimes(1,:) <= rep_img_on);
        postp_img_saccades = find(rep_sactimes(2,:) > rep_img_off);
        rep_sactimes(:,postp_img_saccades) = [];
        rep_sactimes(:,prep_img_saccades) = [];
        
        
        nov_x = eyedat{viewed(1,novrep)}(1,:);
        nov_y = eyedat{viewed(1,novrep)}(2,:);
        nov_x_raw = raw_eyedat{viewed(1,novrep)}(1,:);
        nov_y_raw = raw_eyedat{viewed(1,novrep)}(2,:);
        for sac = 1:size(nov_sactimes,2)
            sacx = nov_x(nov_sactimes(2,sac)) - nov_x(nov_sactimes(1,sac));
            sacy = nov_y(nov_sactimes(2,sac)) - nov_y(nov_sactimes(1,sac));
            nov_trial_by_trial_sacamps(novrep,sac) = sqrt(sacx^2+sacy^2);
            
            sacx = nov_x_raw(nov_sactimes(2,sac)) - nov_x_raw(nov_sactimes(1,sac));
            sacy = nov_y_raw(nov_sactimes(2,sac)) - nov_y_raw(nov_sactimes(1,sac));
            nov_trial_by_trial_sacamps_volts(novrep,sac) = sqrt(sacx^2+sacy^2);
        end
        
        rep_x = eyedat{viewed(2,novrep)}(1,:);
        rep_y = eyedat{viewed(2,novrep)}(2,:);
        rep_x_raw = raw_eyedat{viewed(2,novrep)}(1,:);
        rep_y_raw = raw_eyedat{viewed(2,novrep)}(2,:);
        for sac = 1:size(rep_sactimes,2)
            sacx = rep_x(rep_sactimes(2,sac)) - rep_x(rep_sactimes(1,sac));
            sacy = rep_y(rep_sactimes(2,sac)) - rep_y(rep_sactimes(1,sac));
            rep_trial_by_trial_sacamps(novrep,sac) = sqrt(sacx^2+sacy^2);
            
            sacx = rep_x_raw(rep_sactimes(2,sac)) - rep_x_raw(rep_sactimes(1,sac));
            sacy = rep_y_raw(rep_sactimes(2,sac)) - rep_y_raw(rep_sactimes(1,sac));
            rep_trial_by_trial_sacamps_volts(novrep,sac) = sqrt(sacx^2+sacy^2);
        end
    end
    
    nov_median = round(median(sum(~isnan(nov_trial_by_trial_sacamps)')));
    pre_nov_sac_amps(file,1:nov_median) = nanmean(nov_trial_by_trial_sacamps(:,1:nov_median));
    all_pre_nov_sac_amps = [all_pre_nov_sac_amps; nov_trial_by_trial_sacamps];
    pre_nov_sac_amps_volt(file,1:nov_median) = nanmean(nov_trial_by_trial_sacamps_volts(:,1:nov_median));
    all_pre_nov_sac_amps_volt = [all_pre_nov_sac_amps_volt; nov_trial_by_trial_sacamps_volts];
    
    rep_median = round(median(sum(~isnan(rep_trial_by_trial_sacamps)')));
    pre_rep_sac_amps(file,1:rep_median) = nanmean(rep_trial_by_trial_sacamps(:,1:rep_median));
    all_pre_rep_sac_amps = [all_pre_rep_sac_amps; rep_trial_by_trial_sacamps];
    pre_rep_sac_amps_volt(file,1:rep_median) = nanmean(rep_trial_by_trial_sacamps_volts(:,1:rep_median));
    all_pre_rep_sac_amps_volt = [all_pre_rep_sac_amps_volt; rep_trial_by_trial_sacamps_volts];
end

post_nov_sac_amps = NaN(length(post_files),50);
all_post_nov_sac_amps = [];
post_sac_nov_amps_volt = NaN(length(post_files),50);
all_post_nov_sac_amps_volt = [];
post_rep_sac_amps = NaN(length(post_files),50);
all_post_rep_sac_amps = [];
post_sac_rep_amps_volt = NaN(length(post_files),50);
all_post_rep_sac_amps_volt = [];

for file = 1:length(post_files)
    disp(['Analyzing ' post_files{file}])
    item_set = post_item_files{file};
    
    %---Import Cortex data file---%
    cortexfile = ['R:\Cortex Data\Tobii\' post_files{file}];
    
    cnd_file = [item_set(1:end-4) '.cnd'];
    
    [time_arr,event_arr,eog_arr,epp_arr,~,~]  = get_ALLdata(cortexfile);
    [itmlist,clrchng_locations,first_img_item,imgs] = read_ListRM_itm_and_cnd_files(item_set,cnd_file);
    
    %---Get the calibration for the eye data---%
    %essentially the same as all other task's we use
    numrpt = size(event_arr,2);
    valrptcnt = 0;
    clear per clrchgind
    for rptlop = 1:numrpt
        if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) < first_img_item
            if size(find(event_arr(:,rptlop) == 3)) ~=0 %if clrchng trial was rewarded
                perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
                perendind = find(event_arr(:,rptlop) == 24);
                if length( perbegind) > 1
                    perbegind = perbegind(2);
                    perendind = perendind(2);
                end
                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                begtimdum = time_arr(perbegind,rptlop)-100;
                endtimdum = time_arr(perendind,rptlop);
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    clrchgind(valrptcnt)=rptlop;
                    per(valrptcnt).begsmpind = begtimdum;
                    per(valrptcnt).endsmpind = endtimdum;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                    per(valrptcnt).event = rptlop;
                end
            end
        end
    end
    
    %Don't keep first 6 successful trials. These trials are all central
    %cailbration trials for calibration offset correct
    per(1:5) = [];
    
    clear cnd
    numrpt = size(per,2);
    cnd = zeros(1,numrpt);
    for rptlop = 1:numrpt
        cnd(rptlop)=per(rptlop).cnd-1000;
    end
    
    % Create structures x and y of the corresponding average eye data for each trial
    % instance (l) of each condition (k)
    spacex = [-12,-6,0,6,12];%what actually gets displayed
    spacey = [-8,-4,0,4,8];%what actually gets displayed
    x = cell(length(spacex),length(spacey));%---For Calibration with Eye tracking data with cp2tform---%
    y = cell(length(spacex),length(spacey));
    control = NaN(length(cnd),2);
    for k = 1:length(cnd)
        control(k,:) = clrchng_locations{cnd(k)}';
        
        xi = find(clrchng_locations{cnd(k)}(1) == spacex);
        yi = find(clrchng_locations{cnd(k)}(2) == spacey);
        eyeind = floor(((per(k).begsmpind-1000)/samprate)*2):(floor((per(k).endsmpind-1000)/samprate))*2;
        evenind = eyeind(logical(~rem(eyeind,2)));
        oddind =  eyeind(logical(rem(eyeind,2)));
        x{xi,yi} = [x{xi,yi} mean(eog_arr(oddind,per(k).event))];
        y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))];
    end
    
    %Test for errors%
    count = zeros(length(spacey),length(spacex));
    for xi = 1:length(spacex);
        for yi = 1:length(spacey);
            count(yi,xi) = sum(control(:,1) == spacex(xi) & control(:,2) == spacey(yi));
        end
    end
    if any(count < 10);
        disp('Calibration trial analysis incomplete or error')
        disp('Check number of calibration pionts or task not finished')
    end
    
    clear meanx meany
    for k=1:numel(x)
        xss = x{k};
        low = mean(xss)-std(xss);
        high = mean(xss)+std(xss);
        xss(xss < low) = [];
        xss(xss > high) = [];
        meanx(k)=median(xss);
    end
    for k=1:numel(y)
        yss = y{k};
        low = mean(yss)-std(yss);
        high = mean(yss)+std(yss);
        yss(yss < low) = [];
        yss(yss > high) = [];
        meany(k)=median(y{k});
    end
    
    controlx = [];
    controly = [];
    for i = 1:length(spacex);
        for ii = 1:length(spacey);
            controly = [controly spacey(i)];
            controlx = [controlx spacex(ii)];
        end
    end
    
    %want to change this to MSE estimate with different calibration functions,
    %probably affine and polynomial 3 or 4
    tform = cp2tform([controlx' controly'], [meanx' meany'],'affine');
    tform.forward_fcn = tform.inverse_fcn;
    
    newx = [];
    newy = [];
    figure(101)
    hold on
    for i = 1:length(controlx);
        [x,y] = tformfwd(tform,meanx(i),meany(i));
        plot(x,y,'*r')
        newx(i) = x;
        newy(i) = y;
    end
    hold off
    xlim([-17.5 17.5])
    ylim([-12.5 12.5])
    
    
    figure(102)
    hold on
    plot(meanx,meany,'*r')
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Import Image Viewing Data---%
    
    new_eog_arr = [];
    numrpt = size(event_arr,2);
    new_eog_arr=[];
    valrptcnt = 0;
    clear per clrchgind
    for rptlop = 1:numrpt
        if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) >= first_img_item
            if size(find(event_arr(:,rptlop) == 23)) ~=0 %if image actually turned on
                perbegind = find(event_arr(:,rptlop) == 100,1,'first'); %eye data on
                perendind = find(event_arr(:,rptlop) == 101,1,'first'); %eye data off
                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                begtimdum = time_arr(perbegind,rptlop);
                endtimdum = time_arr(perendind,rptlop);
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    vpcind(valrptcnt)=rptlop;
                    per(valrptcnt).begsmpind = begtimdum;
                    per(valrptcnt).endsmpind = endtimdum;
                    per(valrptcnt).begpos = 1;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                    new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                end
            end
        end
    end
    
    %---get eye data for only when fixation cross or picture is displayed---%
    eyedat = cell(1,length(per));
    for trlop=1:size(per,2)
        trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
        horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
        vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
        
        picstart=1*samprate;
        picend=per(trlop).endsmpind-per(trlop).begsmpind;%added in case sometimes get weird
        %indexing artifacts that are off by 1 index due to cortex having a
        %clock speed with a 1 ms resoultion and the eye data collected at a 5 ms
        %resoultion
        
        try
            eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate));
            eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate));
        catch
            eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate)-1);
            eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate)-1);
        end
    end
    
    %---Recalibrate and automatically scale eye data---%
    fixationstats = cell(1,length(eyedat));
    raw_eyedat = cell(1,length(eyedat));
    for eye = 1:length(eyedat)
        x = eyedat{eye}(1,:);
        y = eyedat{eye}(2,:);
        x_raw = x;
        y_raw = y;
        [x,y] = tformfwd(tform,x,y);
        
        x = 24*x; %convert from cortex dva to pixels
        y = 24*y; %convert from cortex dva to pixels
        x = x+imageX/2;
        y = y+imageY/2;
        
        %---Remove Data outside of image---%
        y_raw(x < -24) = NaN;
        x_raw(x < -24) = NaN;
        x_raw(y < -24) = NaN;
        y_raw(y < -24) = NaN;
        y_raw(x > imageX+24) = NaN;
        x_raw(x > imageX+24)= NaN;
        x_raw(y > imageY+24) = NaN;
        y_raw(y > imageY+24) = NaN;
        
        y(x < -24) = NaN;
        x(x < -24) = NaN;
        y(x > imageX+24) = NaN;
        x(x > imageX+24)= NaN;
        x(y < -24) = NaN;
        y(y < -24) = NaN;
        x(y > imageY+24) = NaN;
        y(y > imageY+24) = NaN;
        
        raw_eyedat{eye} = [x_raw; y_raw];
        eyedat{eye} = [x;y];
        
        [fixationtimes,saccadetimes] = adaptive_VT_Recording_Data2(eyedat{eye});
        
        fixationstats{eye}.fixationtimes = fixationtimes;
        fixationstats{eye}.saccadetimes = saccadetimes;
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
    imgnum = 1:96;
    
    %Let's not even look at images in which the monkey looked away more than the
    %image was displayed for. Also should take care of EOG overflow trials
    viewed(viewed & img_dur > 12000) = 0; %1000 ms buffer for accidental lookaways and cortex timing
    imgnum(viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    viewed(:,viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    
    nov_trial_by_trial_sacamps = NaN(96,50); %each trials saccade amplitudes
    rep_trial_by_trial_sacamps = NaN(96,50); %each trials fixation amplitudes
    nov_trial_by_trial_sacamps_volts = NaN(96,50); %each trials saccade amplitudes
    rep_trial_by_trial_sacamps_volts = NaN(96,50); %each trials fixation amplitudes
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
        
        nov_sactimes = fixationstats{viewed(1,novrep)}.saccadetimes;
        prep_img_saccades = find(nov_sactimes(1,:) <= nov_img_on);
        postp_img_saccades = find(nov_sactimes(2,:) > nov_img_off);
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
        
        rep_sactimes = fixationstats{viewed(2,novrep)}.saccadetimes;
        prep_img_saccades = find(rep_sactimes(1,:) <= rep_img_on);
        postp_img_saccades = find(rep_sactimes(2,:) > rep_img_off);
        rep_sactimes(:,postp_img_saccades) = [];
        rep_sactimes(:,prep_img_saccades) = [];
        
        
        nov_x = eyedat{viewed(1,novrep)}(1,:);
        nov_y = eyedat{viewed(1,novrep)}(2,:);
        nov_x_raw = raw_eyedat{viewed(1,novrep)}(1,:);
        nov_y_raw = raw_eyedat{viewed(1,novrep)}(2,:);
        for sac = 1:size(nov_sactimes,2)
            sacx = nov_x(nov_sactimes(2,sac)) - nov_x(nov_sactimes(1,sac));
            sacy = nov_y(nov_sactimes(2,sac)) - nov_y(nov_sactimes(1,sac));
            nov_trial_by_trial_sacamps(novrep,sac) = sqrt(sacx^2+sacy^2);
            
            sacx = nov_x_raw(nov_sactimes(2,sac)) - nov_x_raw(nov_sactimes(1,sac));
            sacy = nov_y_raw(nov_sactimes(2,sac)) - nov_y_raw(nov_sactimes(1,sac));
            nov_trial_by_trial_sacamps_volts(novrep,sac) = sqrt(sacx^2+sacy^2);
        end
        
        rep_x = eyedat{viewed(2,novrep)}(1,:);
        rep_y = eyedat{viewed(2,novrep)}(2,:);
        rep_x_raw = raw_eyedat{viewed(2,novrep)}(1,:);
        rep_y_raw = raw_eyedat{viewed(2,novrep)}(2,:);
        for sac = 1:size(rep_sactimes,2)
            sacx = rep_x(rep_sactimes(2,sac)) - rep_x(rep_sactimes(1,sac));
            sacy = rep_y(rep_sactimes(2,sac)) - rep_y(rep_sactimes(1,sac));
            rep_trial_by_trial_sacamps(novrep,sac) = sqrt(sacx^2+sacy^2);
            
            sacx = rep_x_raw(rep_sactimes(2,sac)) - rep_x_raw(rep_sactimes(1,sac));
            sacy = rep_y_raw(rep_sactimes(2,sac)) - rep_y_raw(rep_sactimes(1,sac));
            rep_trial_by_trial_sacamps_volts(novrep,sac) = sqrt(sacx^2+sacy^2);
        end
    end
    
    nov_median = round(median(sum(~isnan(nov_trial_by_trial_sacamps)')));
    post_nov_sac_amps(file,1:nov_median) = nanmean(nov_trial_by_trial_sacamps(:,1:nov_median));
    all_post_nov_sac_amps = [all_post_nov_sac_amps; nov_trial_by_trial_sacamps];
    post_nov_sac_amps_volt(file,1:nov_median) = nanmean(nov_trial_by_trial_sacamps_volts(:,1:nov_median));
    all_post_nov_sac_amps_volt = [all_post_nov_sac_amps_volt; nov_trial_by_trial_sacamps_volts];
    
    rep_median = round(median(sum(~isnan(rep_trial_by_trial_sacamps)')));
    post_rep_sac_amps(file,1:rep_median) = nanmean(rep_trial_by_trial_sacamps(:,1:rep_median));
    all_post_rep_sac_amps = [all_post_rep_sac_amps; rep_trial_by_trial_sacamps];
    post_rep_sac_amps_volt(file,1:rep_median) = nanmean(rep_trial_by_trial_sacamps_volts(:,1:rep_median));
    all_post_rep_sac_amps_volt = [all_post_rep_sac_amps_volt; rep_trial_by_trial_sacamps_volts];
end
%%
pre_nov_median = round(median(sum(~isnan(pre_nov_sac_amps'))));
pre_rep_median = round(median(sum(~isnan(pre_rep_sac_amps'))));
post_nov_median = round(median(sum(~isnan(post_nov_sac_amps'))));
post_rep_median = round(median(sum(~isnan(post_rep_sac_amps'))));


figure
subplot(2,2,1)
hold on
plot(nanmean(pre_nov_sac_amps(:,1:pre_nov_median))/24)
plot(nanmean(pre_rep_sac_amps(:,1:pre_rep_median))/24)
plot(nanmean(post_nov_sac_amps(:,1:post_nov_median))/24)
plot(nanmean(post_rep_sac_amps(:,1:post_rep_median))/24)
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
legend('Pre-Nov','Pre-Rep','Post-Nov','Post-Rep')

subplot(2,2,3)
hold on
histogram(all_pre_nov_sac_amps(:)/24,0:0.5:25)
histogram(all_pre_rep_sac_amps(:)/24,0:0.5:25)
histogram(all_post_nov_sac_amps(:)/24,0:0.5:25)
histogram(all_post_rep_sac_amps(:)/24,0:0.5:25)
hold off
xlim([0 20])
xlabel('Saccade Amplitude (dva)')
ylabel('Saccade Count')

subplot(2,2,2)
hold on
plot(nanmean(pre_nov_sac_amps_volt(:,1:pre_nov_median))/24)
plot(nanmean(pre_rep_sac_amps_volt(:,1:pre_rep_median))/24)
plot(nanmean(post_nov_sac_amps_volt(:,1:post_nov_median))/24)
plot(nanmean(post_rep_sac_amps_volt(:,1:post_rep_median))/24)
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (mV)')

subplot(2,2,4)
hold on
histogram(all_pre_nov_sac_amps_volt(:)/24,0:5:200)
histogram(all_pre_rep_sac_amps_volt(:)/24,0:5:200)
histogram(all_post_nov_sac_amps_volt(:)/24,0:5:200)
histogram(all_post_rep_sac_amps_volt(:)/24,0:5:200)
hold off
xlim([0 200])
ylabel('Saccade Amplitude (mV)')
ylabel('Saccade Count')