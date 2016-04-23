function ImportListSqData(clrchng_cortexfile,ListSQ_cortexfile,item_set)
% written December, 2014 by Seth Konig
% Imports a color change calibration file and cortex file for the ListSQ
% task into Matlab. Then the function parses the task data (i.e. when cross
% hairs come on and off and detects fixations and saccades using Cluster
% Fix.
%
% Inputs:
%   1) listsq_cortexfile: cortex file for the ListSQ task
%   2) clrchng_cortexfile: cortex file for the clrchng calibration task
%   3) Item set for the listsq task
%
% Outputs:
%   1) Matlab file with the save fixations and saccade times as well the
%   sequence task data containing the times at which the crosshairs turned
%   on and off

samprate = 5;
imageX = 800;
imageY = 600;

if strcmpi(clrchng_cortexfile(1:2),'PW')
    clrchng_cortexfile = ['R:\Buffalo Lab\Cortex Data\Vivian\' clrchng_cortexfile];
elseif strcmpi(clrchng_cortexfile(1:2),'TT')
    clrchng_cortexfile = ['R:\Buffalo Lab\Cortex Data\Timmy\' clrchng_cortexfile];
elseif strcmpi(clrchng_cortexfile(1:2),'RR')
    clrchng_cortexfile = ['R:\Buffalo Lab\Cortex Data\Red\' clrchng_cortexfile];
end

%%----Color Change Calibration----%%
ITMFile = 'R:\Buffalo Lab\eblab\Cortex Programs\ClrChng\cch25.itm';
CNDFile = 'R:\Buffalo Lab\eblab\Cortex Programs\ClrChng\cch25.cnd';
% this is different becasue the spacing is different and I don't have
% a new item file on the network for the new spacing
ind_spacex = [-6,-3,0,3,6]; %whats on the network
ind_spacey = [-6,-3,0,3,6];%whats on the network
spacex = [-12,-6,0,6,12];%what actually gets displayed
spacey = [-8,-4,0,4,8];%what actually gets displayed
[time_arr,event_arr,eog_arr,~,~,~] = get_ALLdata(clrchng_cortexfile);

itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) <= 189
        if size(find(event_arr(:,rptlop) == 200)) ~=0
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

%Don't keep first 2 calibration pionts these are for offset correction at
%start of task
clear cnd
numrpt = size(per,2);
cnd = zeros(1,numrpt);
for rptlop = 1:numrpt
    cnd(rptlop)=per(rptlop).cnd;
end

% Create structures x and y of the corresponding average eye data for each trial
% instance (l) of each condition (k)

x = cell(length(spacex),length(spacey));%---For Calibration with Eye tracking data with cp2tform---%
y = cell(length(spacex),length(spacey));
control = NaN(length(cnd),2);
clr = ['rgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmk'];
figure
hold on
for k = 1:length(cnd)
    C = textscan(itmfil(itmlist(cnd(k)-1000)+5,:),'%d');
    control(k,:) = C{1}(9:10)';
    
    xi = find(C{1}(9) == ind_spacex);
    yi = find(C{1}(10) == ind_spacey);
    eyeind = floor(((per(k).begsmpind-1000)/samprate)*2):(floor((per(k).endsmpind-1000)/samprate))*2;
    evenind = eyeind(logical(~rem(eyeind,2)));
    oddind =  eyeind(logical(rem(eyeind,2)));
    x{xi,yi} = [x{xi,yi} mean(eog_arr(oddind,per(k).event))];
    y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))];
    plot(mean(eog_arr(oddind,per(k).event)),mean(eog_arr(evenind,per(k).event)),[clr(xi*yi) '+'])
end
if iscell(clrchng_cortexfile)
    title(['Calibration transformation for ' clrchng_cortexfile{1}(end-9:end)])
else
    title(['Calibration transformation for ' clrchng_cortexfile(end-9:end)])
end

%Test for errors%
count = zeros(length(spacey),length(spacex));
for xi = 1:length(spacex);
    for yi = 1:length(spacey);
        count(yi,xi) = sum(control(:,1) == spacex(xi) & control(:,2) == spacey(yi));
    end
end
if any(count < 5);
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
if iscell(clrchng_cortexfile)
    title(['Calibration transformation for ' clrchng_cortexfile{1}(end-9:end)])
else
    title(['Calibration transformation for ' clrchng_cortexfile(end-9:end)])
end
xlim([-17.5 17.5])
ylim([-12.5 12.5])


%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Import LastSQ Task---%
if strcmpi(ListSQ_cortexfile(1:2),'PW')
    cortexfile = ['R:\Buffalo Lab\Cortex Data\Vivian\' ListSQ_cortexfile];
elseif strcmpi(ListSQ_cortexfile(1:2),'TT')
    cortexfile = ['R:\Buffalo Lab\Cortex Data\Timmy\' ListSQ_cortexfile];
elseif strcmpi(ListSQ_cortexfile(1:2),'RR')
    cortexfile = ['R:\Buffalo Lab\Cortex Data\Red\' ListSQ_cortexfile];
end

[time_arr,event_arr,eog_arr,~,~,~]  = get_ALLdata(cortexfile);

[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_set);

numrpt = size(event_arr,2);
new_eog_arr=[];
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if event_arr((find(event_arr(:,rptlop)>500,1,'last')),rptlop)-500 > 1 %1st block is color change always
        cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
        cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
        cnd = event_arr(cndnumind,rptlop)-1000;
        priorcndnumind = find(event_arr(:,rptlop-1) >= 1000 & event_arr(:,rptlop-1) <=2000);
        priorcnd =  event_arr(priorcndnumind,rptlop-1)-1000;
        if any(itmlist(cnd) == sequence_items)  && numel(find(event_arr(:,rptlop) == 3)) == 0
            continue;%sequence trial but not reward
        elseif any(itmlist(cnd) == sequence_items) && cnd == priorcnd
            continue;% same sequence trial was shown  at least 2x in a row so dont want
        elseif ~any(itmlist(cnd) == sequence_items) && numel(find(event_arr(:,rptlop) == 23)) == 0
            continue; %image trial but broke on fixation cross
        end
        perbegind = find(event_arr(:,rptlop) == 100); % eye data starts at 100; might have predictive looking
        perendind = find(event_arr(:,rptlop) == 101); % eye data stops collecting after rewards so can stop here
        cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
        blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
        begtimdum = time_arr(perbegind,rptlop);
        endtimdum = time_arr(perendind(end),rptlop);
        
         test0_on = find(event_arr(:,rptlop) == 23);
            test0_off = find(event_arr(:,rptlop) == 24);
            test1_on = find(event_arr(:,rptlop) == 25);
            test1_off = find(event_arr(:,rptlop) == 26);
            test2_on = find(event_arr(:,rptlop) == 27);
            test2_off = find(event_arr(:,rptlop) == 28);
            test3_on = find(event_arr(:,rptlop) == 29);
            test3_off = find(event_arr(:,rptlop) == 30);
            
            %for looking at trials with break fixations 
            test0_on = time_arr(test0_on,rptlop)-begtimdum;
            if isempty(test0_off);
                test0_off = NaN;
            else
                test0_off = time_arr(test0_off,rptlop)-begtimdum;
            end
            if isempty(test1_off);
                test1_off = NaN;
            else
                test1_off = time_arr(test1_off,rptlop)-begtimdum;
            end
            if isempty(test2_off);
                test2_off = NaN;
            else
                test2_off = time_arr(test2_off,rptlop)-begtimdum;
            end
            if isempty(test3_off);
                test3_off = NaN;
            else
                test3_off = time_arr(test3_off,rptlop)-begtimdum;
            end
            if isempty(test1_on);
                test1_on = NaN;
            else
                test1_on = time_arr(test1_on,rptlop)-begtimdum;
            end
            if isempty(test2_on);
                test2_on = NaN; 
            else
                test2_on = time_arr(test2_on,rptlop)-begtimdum;
            end
            if isempty(test3_on);
                test3_on = NaN;
             else
                test3_on = time_arr(test3_on,rptlop)-begtimdum;
            end
        
        if endtimdum > begtimdum
            valrptcnt = valrptcnt + 1;
            per(valrptcnt).begsmpind = begtimdum;
            per(valrptcnt).endsmpind = endtimdum(1);
            per(valrptcnt).begpos = 1;
            per(valrptcnt).cnd = event_arr(cndnumind,rptlop)-1000;
            per(valrptcnt).blk = event_arr(blknumind,rptlop)-500;
            per(valrptcnt).allval = event_arr(:,rptlop);
            per(valrptcnt).alltim = time_arr(:,rptlop);
            per(valrptcnt).event = rptlop;
            per(valrptcnt).test= ...
                [test0_on test0_off;...
                test1_on test1_off;...
                test2_on test2_off;...
                test3_on test3_off];
            new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
        end
    end
end

%---get eye data for only when fixation cross or picture is displayed---%
eyedat = cell(1,length(per));
cnd=[];
teststart = [];
for trlop=1:size(per,2)
    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
    picstart=1*samprate;
    picend=per(trlop).endsmpind-per(trlop).begsmpind;
    picend(picend > samprate*length(horeog)) = length(horeog)*samprate; 
    
    eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate));
    eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate));
    cnd(trlop)=per(trlop).cnd;
end

%---Recalibrate and automatically scale eye data---%
for eye = 1:length(eyedat)
    x = eyedat{eye}(1,:);
    y = eyedat{eye}(2,:);
    [x,y] = tformfwd(tform,x,y);
    eyedat{eye} = [x;y];
end

% Note needed to redo cluster Fix for small amounts of data: over lap in clusters in velocity/acceleration state space
fixationstats = ClusterFixation_Short(eyedat);
save(['C:\Users\seth.koenig\Documents\MATLAB\ListSQ_RelationalMemory\Cortex Data\' ...
    ListSQ_cortexfile(1:8) '_' ListSQ_cortexfile(end) '-fixation'],...
    'fixationstats','per','item_set')

end