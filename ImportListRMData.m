function ImportListRMData(cortexfile,item_set,figure_dir,data_dir)
% written May, 2015 by Seth Konig
% Imports cortex file for the ListRM (the List task for the relational
% memory projcect) task into Matlab. Then the function parses the task data (i.e. when cross
% hairs come on and off and detects fixations and saccades using Cluster
% Fix.
%
% Inputs:
%   1) cortexfile: cortex file for the ListRM task
%   2) Item set for the listsq task
%   3) figure_dir: where to save a figure of the quality of calibration
% Outputs:
%   1) Matlab file with the save fixations and saccade times as well when
%   the image turned on and off.
%   2) 2 figures of calibration quality saved in the figure_dir

imageX = 800;
imageY = 600;

samprate = 5;
%---Import Cortex data file---%
if strcmpi(cortexfile(1:2),'PW')
    cortexfile = ['\\towerexablox.wanprc.org\Buffalo\Cortex Data\Vivian\' cortexfile];
elseif strcmpi(cortexfile(1:2),'TT')
    cortexfile = ['\\towerexablox.wanprc.org\Buffalo\Cortex Data\Timmy\' cortexfile];
elseif strcmpi(cortexfile(1:2),'RR')
    cortexfile = ['\\towerexablox.wanprc.org\Buffalo\Cortex Data\Red\' cortexfile];
elseif strcmpi(cortexfile(1:2),'TO')
    cortexfile = ['\\towerexablox.wanprc.org\Buffalo\Cortex Data\Tobii\' cortexfile];
elseif  strcmpi(cortexfile(1:2),'MF')
     cortexfile = ['\\towerexablox.wanprc.org\Buffalo\Cortex Data\Manfred\' cortexfile];
end

[time_arr,event_arr,eog_arr,epp_arr,~,~]  = get_ALLdata(cortexfile);

[~,block_1] = find(event_arr == 5001);
if length(block_1) > 1 %file accidentally appended
    disp('Multiple Files Appended')
   event_arr(:,block_1(2):end) = [];
   time_arr(:,block_1(2):end) = [];
   eog_arr(:,block_1(2):end) = [];
   epp_arr(:,block_1(2):end) = [];
end

%---Import Item/Condition file informoation---%
%wrong CND file loaded for these sets. CND file determines calibration
%trial order not which images are displayed that is the item file
if strcmpi(cortexfile(end-9:end),'PW160325.2') || strcmpi(cortexfile(end-9:end),'PW160329.2')
    cnd_file = 'ListRM01.cnd';
else
    cnd_file = [item_set(1:end-4) '.cnd'];
end

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
clr = 'rgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmk';
figure
hold on
for k = 1:length(cnd)
    control(k,:) = clrchng_locations{cnd(k)}';
    
    xi = find(clrchng_locations{cnd(k)}(1) == spacex);
    yi = find(clrchng_locations{cnd(k)}(2) == spacey);
    eyeind = floor(((per(k).begsmpind-1000)/samprate)*2):(floor((per(k).endsmpind-1000)/samprate))*2;
    evenind = eyeind(logical(~rem(eyeind,2)));
    oddind =  eyeind(logical(rem(eyeind,2)));
    x{xi,yi} = [x{xi,yi} mean(eog_arr(oddind,per(k).event))];
    y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))];
    plot(mean(eog_arr(oddind,per(k).event)),mean(eog_arr(evenind,per(k).event)),[clr(xi*yi) '+'])
end
title(['Calibration transformation for ' cortexfile(end-9:end)])
save_and_close_fig(figure_dir,['Calibration_'  cortexfile(end-9:end-2)])

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
figure
hold on
for i = 1:length(controlx);
    plot(controlx(i),controly(i),'r+')
    [x,y] = tformfwd(tform,meanx(i),meany(i));
    plot(x,y,'*b')
    newx(i) = x;
    newy(i) = y;
end
title(['Calibration transformation for ' cortexfile(end-9:end)])
xlim([-17.5 17.5])
ylim([-12.5 12.5])
save_and_close_fig(figure_dir,['Estimated_Calibration_'  cortexfile(end-9:end-2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Import Image Viewing Data---%

new_eog_arr = [];
if ~isempty(epp_arr);%if there is pupil data
    new_epp_arr = [];
end

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
                new_epp_arr = cat(2,new_epp_arr,epp_arr(:,rptlop));
            end
        end
    end
end

%---get eye data for only when fixation cross or picture is displayed---%
eyedat = cell(1,length(per));
if ~isempty(epp_arr);
    pupildata = cell(1,length(per)); %horizontal pupil diameter
end
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
        if picend > 20000
            disp(['EOG overflow. Image Trial ' num2str(trlop) '/' num2str(size(per,2))])
            eyedat{trlop}(1,:) = horeog;
            eyedat{trlop}(2,:) = vrteog;
        else
            try
                eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate)-1);
                eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate)-1);
            catch
                error('Unknown reason EOG and pic duration do not line up')
            end
        end
    end
    try
        trlepp=new_epp_arr(~isnan(new_epp_arr(:,trlop)),trlop); % epp for this trial
        pupildata{trlop} = trlepp(2:2:2*floor(picend/samprate));%odd indexes contains nothing but noise
    catch
        if picend > 20000
            pupildata{trlop} = trlepp(2:2:end); %epp for this trial
        else
            try
                pupildata{trlop} = trlepp(2:2:2*floor(picend/samprate)-1);%odd indexes contains nothing but noise
            catch
                error('Unknown reason EOG and pic duration do not line up')
            end
        end
    end
end

%---Recalibrate and automatically scale eye data---%
fixationstats = {};
raw_eyedat = eyedat;
for eye = 1:length(eyedat)
    x = eyedat{eye}(1,:);
    y = eyedat{eye}(2,:);
    [x,y] = tformfwd(tform,x,y);
    
    x = 24*x; %convert from cortex dva to pixels
    y = 24*y; %convert from cortex dva to pixels
    x = x+imageX/2;
    y = y+imageY/2;
    
    raw_eyedat{eye} = [x;y]; 
    y(x < -24) = NaN;
    x(x < -24) = NaN;
    y(x > imageX+24) = NaN;
    x(x > imageX+24)= NaN;
    x(y < -24) = NaN;
    y(y < -24) = NaN;
    x(y > imageY+24) = NaN;
    y(y > imageY+24) = NaN;
    
    eyedat{eye} = [x;y];
    [fixationstats{eye}] = ClusterFix_Outside(eyedat{eye});
end
save([data_dir cortexfile(end-9:end-2) '_' cortexfile(end) '-fixation'],...
    'fixationstats','per','item_set','pupildata','imgs','raw_eyedat')
end