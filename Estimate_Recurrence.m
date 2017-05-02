% Calculate Recurence
%written by Seth Konig 4.16.2017
clear,clc
fixwin = 3.5; %size of the fixation window/2. Was a width of 7 dva
imageX = 800;
imageY = 600;

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Cortex Data\';
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Image Sets\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Figures\';

%---Vivian---%
% pre_files = {'PW150416.4','PW150417.2','PW150420.2','PW150421.2',...
%     'PW150422.2','PW150423.2','PW150424.2','PW150427.2',...
%     'PW150428.2','PW150429.2','PW150430.2','PW150501.2',...
%     'PW150504.2','PW150505.2','PW150506.2','PW150507.2',...
%     'PW150511.2'};
% 
% post_files = {'PW160310.2','PW160311.2','PW160314.2','PW160315.2','PW160316.2',...
%     'PW160317.2','PW160318.2','PW160321.2','PW160322.2','PW160323.2',...
%     'PW160324.2','PW160325.2','PW160328.2','PW160329.2','PW160330.2'};


%---Red---%
pre_files = {'RR150423.2','RR150424.2','RR150427.2','RR150428.2','RR150429.2',...
    'RR150430.2','RR150501.2','RR150504.2','RR150505.2','RR150506.2',...
    'RR150507.2','RR150508.2','RR150511.2','RR150512.2','RR150513.2',...
    'RR150515.2','RR150518.2','RR150519.2','RR150520.2','RR150521.2'};
post_files = {'RR160324.2','RR160325.1','RR160328.2','RR160330.2',...
                'RR160331.2','RR160401.2','RR160405.2','RR160406.2',...
                'RR160407.2','RR160408.2','RR160411.2','RR160412.2',...
                'RR160413.2','RR160414.2','RR160415.2','RR160418.2',...
                'RR160419.2','RR160420.2','RR160421.2','RR160422.2'};


imageX = 800;
imageY = 600;


nov_pre_recurence_map = zeros(50,50);
rep_pre_recurence_map = zeros(50,50);
rep_pre_recurence_measure = [];
nov_pre_recurence_measure = [];
nov_pre_determinism = [];
rep_pre_determinism = [];
nov_pre_reverse_trace = [];
rep_pre_reverse_trace = [];
nov_pre_laminar = [];
rep_pre_laminar = [];
nov_pre_corm = [];
rep_pre_corm = [];
nov_pre_count = 1;
rep_pre_count = 1;

nov_post_recurence_map = zeros(50,50);
rep_post_recurence_map = zeros(50,50);
rep_post_recurence_measure = [];
nov_post_recurence_measure = [];
nov_post_determinism = [];
rep_post_determinism = [];
nov_post_reverse_trace = [];
rep_post_reverse_trace = [];
nov_post_laminar = [];
rep_post_laminar = [];
nov_post_corm = [];
rep_post_corm = [];
nov_post_count = 1;
rep_post_count = 1;


for file = 1:length(pre_files);
    load([data_dir pre_files{file}(1:8) '_' pre_files{file}(end) '-fixation.mat'])
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
    if strcmpi('PW150416_2',pre_files{file}(1:10))
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
        
        rep_fix = fixationstats{viewed(2,novrep)}.fixations;
        rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
        
        pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        rep_fix(:,post_img_fix) = [];
        rep_fix(:,pre_img_fix) = [];
        
        nov_recurence_map = zeros(50,50);
        N=size(nov_fix,2);
        N(N > 50) = 50;
        [x,y]=meshgrid(1:N);
        i=find(ones(N)); %forms pairs except for self-pairing
        i=[x(i), y(i)];
        dist =sqrt((nov_fix(1,i(:,1))-nov_fix(1,i(:,2))).^2 +...
            (nov_fix(2,i(:,1))-nov_fix(2,i(:,2))).^2);
        dind = find(dist <= 48);
        for d = 1:length(dind);
            nov_recurence_map(i(dind(d),1),i(dind(d),2)) = nov_recurence_map(i(dind(d),1),i(dind(d),2))+1;
        end
        R = sum(sum(triu(nov_recurence_map)))-N;
        nov_pre_recurence_measure(nov_pre_count) = 100*2*R/N/(N-1);
        
                
        determinism_count = 0;
        reverse_trace = 0;
        for i = 1:N-1
            for ii = 2:N
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
        nov_pre_determinism(nov_pre_count) = 100*determinism_count/R;
        nov_pre_reverse_trace(nov_pre_count) = 100*reverse_trace/R;
        laminar = nov_recurence_map;
        laminar(tril(laminar) == 1) = 0;
        nov_pre_laminar(nov_pre_count) = 100*(sum(sum(laminar,1) > 1)+  sum(sum(laminar,2) > 1))/2/R;
        
        if sum(laminar) < 3
            nov_pre_corm(:,nov_pre_count)  = [NaN; NaN];
        else
            [xcm ycm] = centroid(laminar);
            nov_pre_corm(:,nov_pre_count) = [xcm;ycm];
        end
        
        nov_pre_count = nov_pre_count+1;
        
        rep_recurence_map = zeros(50,50);
        N=size(rep_fix,2);
        N(N > 50) = 50;
        [x,y]=meshgrid(1:N);
        i=find(ones(N)); %forms pairs except for self-pairing
        i=[x(i), y(i)];
        dist =sqrt((rep_fix(1,i(:,1))-rep_fix(1,i(:,2))).^2 +...
            (rep_fix(2,i(:,1))-rep_fix(2,i(:,2))).^2);
        dind = find(dist <= 48);
        for d = 1:length(dind);
            rep_recurence_map(i(dind(d),1),i(dind(d),2)) = rep_recurence_map(i(dind(d),1),i(dind(d),2))+1;
        end
        R = sum(sum(triu(rep_recurence_map)))-N;
        rep_pre_recurence_measure(rep_pre_count) = 100*2*R/N/(N-1);
        
        determinism_count = 0;
        reverse_trace = 0;
        for i = 1:N-1
            for ii = 2:N
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
        rep_pre_determinism(rep_pre_count) = 100*determinism_count/R;
        rep_pre_reverse_trace(rep_pre_count) = 100*reverse_trace/R;
        
        laminar = rep_recurence_map;
        laminar(tril(laminar) == 1) = 0;
        rep_pre_laminar(rep_pre_count) = 100*(sum(sum(laminar,1) > 1)+  sum(sum(laminar,2) > 1))/2/R;
        
        if sum(laminar) < 3
            rep_pre_corm(:,rep_pre_count)  = [NaN; NaN];
        else
            [xcm ycm] = centroid(laminar);
            rep_pre_corm(:,rep_pre_count) = [xcm;ycm];
        end
        
        rep_pre_count = rep_pre_count+1;

        
        nov_pre_recurence_map = nov_pre_recurence_map + nov_recurence_map;
        rep_pre_recurence_map = rep_pre_recurence_map + rep_recurence_map;
    end
end
%%
for file = 1:length(post_files);
    load([data_dir post_files{file}(1:8) '_' post_files{file}(end) '-fixation.mat'])
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
    if strcmpi('PW150416_2',post_files{file}(1:10))
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
        
        rep_fix = fixationstats{viewed(2,novrep)}.fixations;
        rep_fixtimes = fixationstats{viewed(2,novrep)}.fixationtimes;
        
        pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        rep_fix(:,post_img_fix) = [];
        rep_fix(:,pre_img_fix) = [];
        
        nov_recurence_map = zeros(50,50);
        N=size(nov_fix,2);
        N(N > 50) = 50;
        [x,y]=meshgrid(1:N);
        i=find(ones(N)); %forms pairs except for self-pairing
        i=[x(i), y(i)];
        dist =sqrt((nov_fix(1,i(:,1))-nov_fix(1,i(:,2))).^2 +...
            (nov_fix(2,i(:,1))-nov_fix(2,i(:,2))).^2);
        dind = find(dist <= 48);
        for d = 1:length(dind);
            nov_recurence_map(i(dind(d),1),i(dind(d),2)) = nov_recurence_map(i(dind(d),1),i(dind(d),2))+1;
        end
        R = sum(sum(triu(nov_recurence_map)))-N;
        nov_post_recurence_measure(nov_post_count) = 100*2*R/N/(N-1);
        
                
        determinism_count = 0;
        reverse_trace = 0;
        for i = 1:N-1
            for ii = 2:N
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
        nov_post_determinism(nov_post_count) = 100*determinism_count/R;
        nov_post_reverse_trace(nov_post_count) = 100*reverse_trace/R;
        laminar = nov_recurence_map;
        laminar(tril(laminar) == 1) = 0;
        nov_post_laminar(nov_post_count) = 100*(sum(sum(laminar,1) > 1)+  sum(sum(laminar,2) > 1))/2/R;
        
        if sum(laminar) < 3
            nov_post_corm(:,nov_post_count)  = [NaN; NaN];
        else
            [xcm ycm] = centroid(laminar);
            nov_post_corm(:,nov_post_count) = [xcm;ycm];
        end
        
        nov_post_count = nov_post_count+1;
        
        rep_recurence_map = zeros(50,50);
        N=size(rep_fix,2);
        N(N > 50) = 50;
        [x,y]=meshgrid(1:N);
        i=find(ones(N)); %forms pairs except for self-pairing
        i=[x(i), y(i)];
        dist =sqrt((rep_fix(1,i(:,1))-rep_fix(1,i(:,2))).^2 +...
            (rep_fix(2,i(:,1))-rep_fix(2,i(:,2))).^2);
        dind = find(dist <= 48);
        for d = 1:length(dind);
            rep_recurence_map(i(dind(d),1),i(dind(d),2)) = rep_recurence_map(i(dind(d),1),i(dind(d),2))+1;
        end
        R = sum(sum(triu(rep_recurence_map)))-N;
        rep_post_recurence_measure(rep_post_count) = 100*2*R/N/(N-1);
        
        determinism_count = 0;
        reverse_trace = 0;
        for i = 1:N-1
            for ii = 2:N
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
        rep_post_determinism(rep_post_count) = 100*determinism_count/R;
        rep_post_reverse_trace(rep_post_count) = 100*reverse_trace/R;
        
        laminar = rep_recurence_map;
        laminar(tril(laminar) == 1) = 0;
        rep_post_laminar(rep_post_count) = 100*(sum(sum(laminar,1) > 1)+  sum(sum(laminar,2) > 1))/2/R;
        
        if sum(laminar) < 3
            rep_post_corm(:,rep_post_count)  = [NaN; NaN];
        else
            [xcm ycm] = centroid(laminar);
            rep_post_corm(:,rep_post_count) = [xcm;ycm];
        end
        
        rep_post_count = rep_post_count+1;

        
        nov_post_recurence_map = nov_post_recurence_map + nov_recurence_map;
        rep_post_recurence_map = rep_post_recurence_map + rep_recurence_map;
    end
end
%%

figure
%---For Pre-Novel---%
rm = nov_pre_recurence_map;
% rm = 100*rm/sum(rm(eye(size(rm,1))==1));
rm = 100*rm/rm(1);
rm(find(eye(size(rm)))) = NaN;
id = eye(size(rm));
id = [id(2:end,:); zeros(1,50)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [zeros(1,50);id(1:end-1,:)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [id(2:end,:); zeros(1,50)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [zeros(2,50);id(1:end-2,:)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [id(3:end,:); zeros(2,50)];
rm(find(id)) = NaN;

upperbound = prctile(rm(:),99);
lowerbound = prctile(rm(:),1.0);

subplot(2,2,1)
h = imagesc(rm(1:25,1:25));
axis xy, set(h,'alphadata',~isnan(rm(1:25,1:25))); 
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
box off
colorbar
colormap('viridis') 
title('Pre-Novel')

%---For Pre-Repeat---%
rm = rep_pre_recurence_map;
% rm = 100*rm/sum(rm(eye(size(rm,1))==1));
rm = 100*rm/rm(1);
rm(find(eye(size(rm)))) = NaN;
id = eye(size(rm));
id = [id(2:end,:); zeros(1,50)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [zeros(1,50);id(1:end-1,:)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [id(2:end,:); zeros(1,50)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [zeros(2,50);id(1:end-2,:)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [id(3:end,:); zeros(2,50)];
rm(find(id)) = NaN;

upperbound = prctile(rm(:),99);
lowerbound = prctile(rm(:),1.0);

subplot(2,2,2)
h = imagesc(rm(1:25,1:25));
axis xy, set(h,'alphadata',~isnan(rm(1:25,1:25))); 
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
colormap('viridis') 
colorbar
box off
title('Pre-Repeat')


%---For Pre-Repeat minus Pre-Novel---%
rm1 = nov_pre_recurence_map;
% rm1 = 100*rm1/sum(rm1(eye(size(rm1,1))==1));
rm1 = 100*rm1/rm1(1);

rm2 = rep_pre_recurence_map;
%rm2 = 100*rm2/sum(rm2(eye(size(rm2,1))==1));
rm2 = 100*rm2/rm2(1);

rm = rm2-rm1;
rm(find(eye(size(rm)))) = 0;

upperbound = prctile(rm(:),99);
lowerbound = prctile(rm(:),1.0);

subplot(2,2,3)
h = imagesc(rm(1:25,1:25));
axis xy, set(h,'alphadata',~isnan(rm(1:25,1:25))); 
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
colormap('viridis') 
box off
colorbar
title('Pre-Repeat minus Pre-Novel')


[nvpre,nvprei] = hist(nov_pre_recurence_measure,[5:2.5:50]);
[rppre,rpprei] = hist(rep_pre_recurence_measure,[5:2.5:50]);

subplot(2,2,4)
plot(nvprei,nvpre,'b')
hold on
plot(rpprei,rppre,'r')
hold off
legend('Novel','Repeat')
title('Distribution of Recurrence Measures')
xlabel('Recurrence in %')
ylabel('Image Count')
xlim([5 40])
box off

subtitle([pre_files{1}(1:2) ': Pre-lesion recurence'])
%%
figure
%---For post-Novel---%
rm = nov_post_recurence_map;
rm(find(eye(size(rm)))) = NaN;
id = eye(size(rm));
id = [id(2:end,:); zeros(1,50)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [zeros(1,50);id(1:end-1,:)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [id(2:end,:); zeros(1,50)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [zeros(2,50);id(1:end-2,:)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [id(3:end,:); zeros(2,50)];
rm(find(id)) = NaN;

upperbound = prctile(rm(:),99);
lowerbound = prctile(rm(:),1.0);

subplot(2,2,1)
h = imagesc(rm(1:25,1:25));
axis xy, set(h,'alphadata',~isnan(rm(1:25,1:25))); 
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
box off
colorbar
colormap('viridis') 
title('post-Novel')

%---For post-Repeat---%
rm = rep_post_recurence_map;
rm(find(eye(size(rm)))) = NaN;
id = eye(size(rm));
id = [id(2:end,:); zeros(1,50)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [zeros(1,50);id(1:end-1,:)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [id(2:end,:); zeros(1,50)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [zeros(2,50);id(1:end-2,:)];
rm(find(id)) = NaN;
id = eye(size(rm));
id = [id(3:end,:); zeros(2,50)];
rm(find(id)) = NaN;

upperbound = prctile(rm(:),99);
lowerbound = prctile(rm(:),1.0);

subplot(2,2,2)
h = imagesc(rm(1:25,1:25));
axis xy, set(h,'alphadata',~isnan(rm(1:25,1:25))); 
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
colormap('viridis') 
colorbar
box off
title('post-Repeat')


%---For post-Repeat minus post-Novel---%
rm = rep_post_recurence_map-nov_post_recurence_map;
rm(find(eye(size(rm)))) = NaN;

upperbound = prctile(rm(:),99);
lowerbound = prctile(rm(:),1.0);

subplot(2,2,3)
h = imagesc(rm(1:25,1:25));
axis xy, set(h,'alphadata',~isnan(rm(1:25,1:25))); 
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
colormap('viridis') 
box off
colorbar
title('post-Repeat minus post-Novel')


[nvpost,nvposti] = hist(nov_post_recurence_measure,[5:2.5:50]);
[rppost,rpposti] = hist(rep_post_recurence_measure,[5:2.5:50]);

subplot(2,2,4)
plot(nvposti,nvpost,'b')
hold on
plot(rpposti,rppost,'r')
hold off
legend('Novel','Repeat')
title('Distribution of Recurrence Measures')
xlabel('Recurrence in %')
ylabel('Image Count')
xlim([5 40])
box off

subtitle([post_files{1}(1:2) ': post-lesion recurence'])

%%

figure
%---For Novel---%
rm_post = nov_post_recurence_map;
rm_post(find(eye(size(rm_post)))) = NaN;
rm_pre = nov_pre_recurence_map;
rm_pre(find(eye(size(rm_pre)))) = NaN;

rm = rm_pre-rm_post;
upperbound = prctile(rm(:),99);
lowerbound = prctile(rm(:),1.0);

subplot(2,2,1)
h = imagesc(rm(1:25,1:25));
axis xy, set(h,'alphadata',~isnan(rm(1:25,1:25))); 
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
box off
colorbar
colormap('viridis') 
title('Novel pre-post')

%---For Repeat---%
rm_post = rep_post_recurence_map;
rm_post(find(eye(size(rm_post)))) = NaN;
rm_pre = rep_pre_recurence_map;
rm_pre(find(eye(size(rm_pre)))) = NaN;

rm = rm_pre-rm_post;
upperbound = prctile(rm(:),99);
lowerbound = prctile(rm(:),1.0);

subplot(2,2,2)
h = imagesc(rm(1:25,1:25));
axis xy, set(h,'alphadata',~isnan(rm(1:25,1:25))); 
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
box off
colorbar
colormap('viridis') 
title('Repeat pre-post')


%---For Repeat-Novel ---%
rm_post = rep_post_recurence_map-nov_post_recurence_map;
rm_post(find(eye(size(rm_post)))) = NaN;
rm_pre = rep_pre_recurence_map-nov_pre_recurence_map;
rm_pre(find(eye(size(rm_pre)))) = NaN;

rm = rm_pre-rm_post;
upperbound = prctile(rm(:),99);
lowerbound = prctile(rm(:),1.0);

subplot(2,2,3)
h = imagesc(rm(1:25,1:25));
axis xy, set(h,'alphadata',~isnan(rm(1:25,1:25))); 
caxis([lowerbound upperbound])
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
box off
colorbar
colormap('viridis') 
title('Repeat-Novel pre-post')

subplot(2,2,4)
hold on
plot(nvprei,nvpre,'b')
plot(rpprei,rppre,'r')
plot(nvposti,nvpost,'g')
plot(rpposti,rppost,'m')
hold off
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat')
title('Distribution of Recurrence Measures')
xlabel('Recurrence in %')
ylabel('Image Count')
xlim([5 40])
box off


subtitle([post_files{1}(1:2) ': pre vs post lesion recurrence'])

%%
diag4 = [];

for d = 2:5
id = eye(size(rm));
id = [id(d:end,:); zeros(d-1,50)];
diag4 = [diag4; find(id == 1)];

end
   
cluster1 = zeros(50,50);
cluster1(diag4) = 1;
cluster1(:,1:5) = 0;
cluster1(1:5,:) = 0;
cluster1(:,21:end) = 0;
cluster1(1,21:end) = 0;

diag7 = [];
for d = 2:7
id = eye(size(rm));
id = [id(d:end,:); zeros(d-1,50)];
diag7 = [diag7; find(id == 1)];

end

cluster2 = zeros(50,50);
cluster2(diag7) = 2;
cluster2(cluster2 == 0) = 1;
cluster2(cluster2 == 2) = 0;
cluster2(:,1:5) = 0;
cluster2(1:5,:) = 0;
cluster2(:,21:end) = 0;
cluster2(1,21:end) = 0;
cluster2(tril(cluster2)== 1) = 0;

%%
% rep_pre_recurence_measure = [];
% nov_pre_recurence_measure = [];
% nov_pre_determinism = [];
% rep_pre_determinism = [];
% nov_pre_reverse_trace = [];
% rep_pre_reverse_trace = [];
% nov_pre_laminar = [];
% rep_pre_laminar = [];
% nov_pre_corm = [];
% rep_pre_corm = [];
% nov_pre_count = 1;
% rep_pre_count = 1;
%%
pre_image_count = length(nov_pre_recurence_measure);
post_image_count = length(nov_post_recurence_measure);

figure
subplot(3,3,1)
hold on
[x,xi] = hist(nov_pre_recurence_measure,0:0.25:100);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/pre_image_count));
[x,xi] = hist(rep_pre_recurence_measure,0:0.25:100);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/pre_image_count));
[x,xi] = hist(nov_post_recurence_measure,0:0.25:100);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/post_image_count));
[x,xi] = hist(rep_post_recurence_measure,0:0.25:100);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/post_image_count));
hold off
xlabel('Recurrence (%)')
ylabel('Image Count')
legend('Nov. Pre','Rep. Pre','Nov. Post','Rep Post')
title('Recurrence Rate')
xlim([0 20])

subplot(3,3,2)
hold on
[x,xi] = hist(nov_pre_determinism,0:0.25:100);
x(1) = max(x(2:end))/2;
plot(xi,filtfilt(ones(1,8)*1/8,1,x/pre_image_count));
[x,xi] = hist(rep_pre_determinism,0:0.25:100);
x(1) = max(x(2:end))/2;
plot(xi,filtfilt(ones(1,8)*1/8,1,x/pre_image_count));
[x,xi] = hist(nov_post_determinism,0:0.25:100);
x(1) = max(x(2:end))/2;
plot(xi,filtfilt(ones(1,8)*1/8,1,x/post_image_count));
[x,xi] = hist(rep_post_determinism,0:0.25:100);
x(1) = max(x(2:end))/2;
plot(xi,filtfilt(ones(1,8)*1/8,1,x/post_image_count));
hold off
xlabel('Determinism (%)')
ylabel('Image Count')
title('Retrace Rate')
xlim([0 40])

subplot(3,3,3)
hold on
[x,xi] = hist(nov_pre_reverse_trace,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/pre_image_count));
[x,xi] = hist(rep_pre_reverse_trace,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/pre_image_count));
[x,xi] = hist(nov_post_reverse_trace,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/post_image_count));
[x,xi] = hist(rep_post_reverse_trace,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/post_image_count));
hold off
xlabel('Determinism (%)')
ylabel('Image Count')
title('Reverse Retrace Rate')
xlim([0 60])

subplot(3,3,4)
hold on
[x,xi] = hist(nov_pre_laminar,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/pre_image_count));
[x,xi] = hist(rep_pre_laminar,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/pre_image_count));
[x,xi] = hist(nov_post_laminar,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/post_image_count));
[x,xi] = hist(rep_post_laminar,0:0.25:100);
x(1) = max(x(2:end))/4;
plot(xi,filtfilt(ones(1,12)*1/12,1,x/post_image_count));
hold off
xlabel('Laminarity (%)')
ylabel('Image Count')
title('Laminar Re-Re-Exploration')
xlim([0 50])

subplot(3,3,5)
hold on
[x,xi] = hist(nov_pre_corm(1,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/pre_image_count));
[x,xi] = hist(rep_pre_corm(1,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/pre_image_count));
[x,xi] = hist(nov_post_corm(1,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/post_image_count));
[x,xi] = hist(rep_post_corm(1,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/post_image_count));
hold off
xlabel('Ordinal Fixation #')
ylabel('Image Count')
title('Horizontal Center of Mass')
xlim([0 30])

center_nov_pre_count = sum(~isnan(nov_pre_corm(1,:)));
center_rep_pre_count = sum(~isnan(rep_pre_corm(1,:)));
center_nov_post_count = sum(~isnan(nov_post_corm(1,:)));
center_rep_post_count = sum(~isnan(rep_post_corm(1,:)));

subplot(3,3,6)
hold on
[x,xi] = hist(nov_pre_corm(2,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/center_nov_pre_count));
[x,xi] = hist(rep_pre_corm(2,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/center_rep_pre_count));
[x,xi] = hist(nov_post_corm(2,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/center_nov_post_count));
[x,xi] = hist(rep_post_corm(2,:),0:1:30);
plot(xi,filtfilt(ones(1,2)*1/2,1,x/center_rep_post_count));
hold off
xlabel('Ordinal Fixation #')
ylabel('Image Count')
title('Vertical Center of Mass')
xlim([0 30])


subplot(3,3,7)
hold on
[x,xi] = hist(abs(nov_pre_corm(1,:)-nov_pre_corm(2,:)),0:0.25:20);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/center_nov_pre_count));
[x,xi] = hist(abs(rep_pre_corm(1,:)-rep_pre_corm(2,:)),0:0.25:20);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/center_rep_pre_count));
[x,xi] = hist(abs(nov_post_corm(1,:)-nov_post_corm(2,:)),0:0.25:20);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/center_nov_post_count));
[x,xi] = hist(abs(rep_post_corm(1,:)-rep_post_corm(2,:)),0:0.25:20);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/center_rep_post_count));
hold off
xlabel('Local --> Global')
ylabel('Image Count')
title('Distance from Center: Local vs Global')


pre_nov = abs(nov_pre_corm(1,:)-nov_pre_corm(2,:));
pre_rep = abs(rep_pre_corm(1,:)-rep_pre_corm(2,:));
post_nov = abs(nov_post_corm(1,:)-nov_post_corm(2,:));
post_rep = abs(rep_post_corm(1,:)-rep_post_corm(2,:));

subplot(3,3,8)
hold on
[x,xi] = hist(pre_rep-pre_nov,-10:0.25:10);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/sum(~isnan(pre_nov-pre_rep))));
[x,xi] = hist(post_rep-post_nov,-10:0.25:10);
plot(xi,filtfilt(ones(1,4)*1/4,1,x/sum(~isnan(post_nov-post_rep))));
hold off
xlabel('Change in Local --> Global')
ylabel('Image Count')
title('Change from Global to Local')
legend('Pre Rep-Nov','Post Rep-Nov')



subtitle(['All Recurrence Measurse for ' pre_files{1}(1:2)])
%%
