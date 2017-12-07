% test to see if ListRM was running correctly without a monkey present
cd('C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\');
data_dir = '\\towerexablox.wanprc.org\Buffalo\Cortex Data\Vivian\';
[time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata([data_dir 'PW150417.2']);

ITMFile = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Item and CND files\ListRM10.itm';
CNDFile = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Item and CND files\ListRM10.cnd';


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
    if event_arr((find(event_arr(:,rptlop)>500,1,'last')),rptlop)-500 == 2 %1st block is color change always
        if itmlist(event_arr(find(event_arr(:,rptlop)>1000,1,'last'),rptlop)-1000) <= 75 %clrchng trials
            if size(find(event_arr(:,rptlop) == 200)) ~=0
                perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
                perendind = find(event_arr(:,rptlop) == 24);
                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                begtimdum = time_arr(perbegind,rptlop)-100;
                endtimdum = time_arr(perendind,rptlop);
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    clrchgind(valrptcnt)=rptlop;
                    per(valrptcnt).begsmpind = begtimdum;
                    per(valrptcnt).endsmpind = endtimdum;
                    per(valrptcnt).begpos = 1;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                end
            end
        end
    end
end

numrpt = size(per,2);
for rptlop = 1:numrpt
    clrcnd(rptlop)=itmlist(per(rptlop).cnd-1000);
end
count = NaN(1,75);
for item = 3:3:75;
   count(item) = sum(clrcnd == item);  
end
count(isnan(count)) = [];
%%
numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if event_arr((find(event_arr(:,rptlop)>500,1,'last')),rptlop)-500 == 2 %1st block is color change always
        if itmlist(event_arr(find(event_arr(:,rptlop)>1000,1,'last'),rptlop)-1000) > 75 %clrchng trials
            if size(find(event_arr(:,rptlop) == 200)) ~=0
                perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
                perendind = find(event_arr(:,rptlop) == 24);
                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                begtimdum = time_arr(perbegind,rptlop)-100;
                endtimdum = time_arr(perendind,rptlop);
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    clrchgind(valrptcnt)=rptlop;
                    per(valrptcnt).begsmpind = begtimdum;
                    per(valrptcnt).endsmpind = endtimdum;
                    per(valrptcnt).begpos = 1;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                    
                    cross_on = find(event_arr(:,rptlop) == 9);
                    cross_off = find(event_arr(:,rptlop) == 10);
                    img_on = find(event_arr(:,rptlop) == 23);
                    img_off = find(event_arr(:,rptlop) == 24);
                    
                    per(valrptcnt).crossfixdur = time_arr(cross_off,rptlop)-time_arr(cross_on,rptlop);
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop)-1000;
                    per(valrptcnt).imgdur = time_arr(img_off,rptlop)-time_arr(img_on,rptlop);
                end
            end
        end
    end
end


imgdur = [];
crossdur = [];
for trial = 1:length(per);
    imgdur = [imgdur per(trial).imgdur];
    crossdur = [crossdur per(trial).crossfixdur];
end
%%
imgs = zeros(2,96);
total_images_trails = zeros(2,96); 
img_trial = 1;
for trial = 1:length(per);
    cndline = textscan(cndfil(per(trial).cnd+1,:),'%d');
    imgnum = cndline{1}(end)-75;
    if imgnum >= 1;
        if imgs(1,imgnum) == 0;
            imgs(1,imgnum) = img_trial;
            total_images_trails(1,imgnum) = total_images_trails(1,imgnum)+1; 
        else
            imgs(2,imgnum) = img_trial;
              total_images_trails(2,imgnum) = total_images_trails(2,imgnum)+1; 
        end
        img_trial = img_trial+1;
    end
end