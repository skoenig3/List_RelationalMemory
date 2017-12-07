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


nov_1backs = [];
nov_2backs = [];
nov_3backs = [];

rep_1backs = [];
rep_2backs = [];
rep_3backs = [];


all_nov_recurence_map = zeros(40,40);
all_rep_recurence_map = zeros(40,40);
rep_recurence_measure = [];
nov_recurence_measure = [];
nov_determinism = [];
rep_determinism = [];
nov_reverse_trace = [];
rep_reverse_trace = [];
nov_laminar = [];
rep_laminar = [];
nov_corm = [];
rep_corm = [];
nov_count = 1;
rep_count = 1;

nov_back_durs1 = [];
rep_back_durs1 = [];

nov_all_fix_durs = NaN(4000,40);
rep_all_fix_durs = NaN(4000,40);
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
    
    img_count = 0;
    n1backs = 0;
    n2backs = 0;
    n3backs = 0;
    r1backs = 0;
    r2backs = 0;
    r3backs = 0;
    
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
    
    for novrep = 1:size(viewed,2);
        %---Grab Important Vairables---%
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
        
        nov_fix = fixationstats{viewed(1,novrep)}.fixations;
        nov_fixtimes = fixationstats{viewed(1,novrep)}.fixationtimes;
        
        img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
        post_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
        nov_fix(:,post_img_fix) = [];
        nov_fix(:,img_fix) = [];
        nov_fixtimes(:,post_img_fix) = [];
        nov_fixtimes(:,img_fix ) = [];
        nov_fix_durs = diff(nov_fixtimes,1)+1;
        
        randind = randperm(size(nov_fix,2));
        nov_fix = nov_fix(:,randind);
        
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
        
        
        rep_fix = fixationstats{viewed(2,novrep)}.fixations;
        rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
        
        img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        rep_fix(:,post_img_fix) = [];
        rep_fix(:,img_fix) = [];
        rep_fixtimes(:,post_img_fix) = [];
        rep_fixtimes(:,img_fix ) = [];
        rep_fix_durs = diff(rep_fixtimes,1)+1;
        
        randind = randperm(size(rep_fix,2));
        rep_fix = rep_fix(:,randind);
        
        %----All Fixation Durations---%
        nov_all_fix_durs(nov_count,1:size(nov_fixtimes,2)) = diff(nov_fixtimes,1)+1;
        rep_all_fix_durs(rep_count,1:size(rep_fixtimes,2)) = diff(rep_fixtimes,1)+1;
        

        img_count = img_count+1;
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
            for b = 1:length(b1)
                if b1(b) > 2 && b1(b) < length(nov_fix_durs)-1
                    nov_back_durs1 = [nov_back_durs1 ...
                        [b1(b)-2 nov_fix_durs(b1(b)-2)... %2 fixations back
                        b1(b)-1 nov_fix_durs(b1(b)-1)...%1 fixation back i.e. prior fixation
                        b1(b) nov_fix_durs(b1(b)) ...%out fixation
                        b1(b)+1 nov_fix_durs(b1(b)+1)...%return fixation
                        b1(b)+2 nov_fix_durs(b1(b)+2)...%2 fixations forward
                        this_monkey]'];
                end
            end
        end
        
        n1backs = n1backs+sum(nov_recurence_map(id1 == 1))/N;
        n2backs = n2backs+sum(nov_recurence_map(id2 == 1))/N;
        n3backs = n3backs+sum(nov_recurence_map(id3 == 1))/N;
                
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
 
        if sum(laminar) < 3
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
            for b = 1:length(b1)
                if b1(b) > 2 && b1(b) < length(rep_fix_durs)-1
                    rep_back_durs1 = [rep_back_durs1 ...
                        [b1(b)-2 rep_fix_durs(b1(b)-2)... %2 fixations back
                        b1(b)-1 rep_fix_durs(b1(b)-1)...%1 fixation back i.e. prior fixation
                        b1(b) rep_fix_durs(b1(b)) ...%out fixation
                        b1(b)+1 rep_fix_durs(b1(b)+1)...%return fixation
                        b1(b)+2 rep_fix_durs(b1(b)+2)...%2 fixations forward
                        this_monkey]'];
                end
            end
        end
        
        r1backs = r1backs+sum(rep_recurence_map(id1 == 1))/N;
        r2backs = r2backs+sum(rep_recurence_map(id2 == 1))/N;
        r3backs = r3backs+sum(rep_recurence_map(id3 == 1))/N;
        
        
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
        
        if sum(laminar) < 3
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
    
    
    nov_1backs = [nov_1backs n1backs/img_count];
    nov_2backs = [nov_2backs n2backs/img_count];
    nov_3backs = [nov_3backs n3backs/img_count];
    
    rep_1backs = [rep_1backs r1backs/img_count];
    rep_2backs = [rep_2backs r2backs/img_count];
    rep_3backs = [rep_3backs r3backs/img_count];
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
ylabel('Percentage of Fixations')
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
hold off
xlabel('Recurrence (%)')
ylabel('Proportion of Images') 
legend('Novel','Repeat')
title('Recurrence Rate')
xlim([0 20])

subplot(3,3,2)
hold on
[x,xi] = hist(nov_determinism,0:0.25:100);
x(1) = max(x(2:end))/2;
plot(xi,filtfilt(ones(1,8)*1/8,1,x/image_count));
[x,xi] = hist(rep_determinism,0:0.25:100);
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
hold off
xlabel('Laminarity (%)')
ylabel('Proportion of Images') 
title('Laminar Re-Re-Exploration')
xlim([0 50])

subplot(3,3,5)
hold on
[x,xi] = hist(nov_corm(1,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/image_count));
[x,xi] = hist(rep_corm(1,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/image_count));
hold off
xlabel('Ordinal Fixation #')
ylabel('Proportion of Images') 
title('Horizontal Center of Mass')
xlim([0 30])

center_nov_count = sum(~isnan(nov_corm(1,:)));
center_rep_count = sum(~isnan(rep_corm(1,:)));

subplot(3,3,6)
hold on
[x,xi] = hist(nov_corm(2,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/center_nov_count));
[x,xi] = hist(rep_corm(2,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/center_rep_count));
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
end

%%
nov_redisual_duration = NaN(5,length(nov_back_durs1));
for n = 1:length(nov_back_durs1)
    this_monkey = nov_back_durs1(end,n);
    if nov_back_durs1(9,n) <= median_num_fix(1,this_monkey)
        for f = 1:5
            nov_redisual_duration(f,n) = nov_back_durs1(f*2,n)-nanmean(nov_fix_dur_avg(this_monkey,f*2-1));
        end
    else
        continue
    end
end
%%
rep_redisual_duration = NaN(5,length(rep_back_durs1));
for n = 1:length(rep_back_durs1)
    this_monkey = rep_back_durs1(end,n);
    if rep_back_durs1(9,n) <= median_num_fix(1,this_monkey)
        for f = 1:5
            rep_redisual_duration(f,n) = rep_back_durs1(f*2,n)-nanmean(rep_fix_dur_avg(this_monkey,f*2-1));
        end
    else
        continue
    end
end

figure
plot(nanmean(nov_redisual_duration'))
hold on
plot(nanmean(rep_redisual_duration'))