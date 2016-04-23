% fixation duration by image #

all_nov_fix_durs = [];
all_rep_fix_durs = [];
for set = 1:size(first_fix,2)
    all_nov_fix_durs = [all_nov_fix_durs; first_fix{set}(1,:)];
    all_rep_fix_durs = [all_rep_fix_durs; first_fix{set}(2,:)]; 
end

figure
hold on
plot(nanmean(all_nov_fix_durs))
plot(nanmean(all_rep_fix_durs),'r')
plot([1 1],[130 220],'k--')   %first image block 1
plot([17 17],[130 220],'k--') %first image block 2
plot([33 33],[130 220],'k--') %first image block 3
plot([49 49],[130 220],'k--') %first image block 4
plot([65 65],[130 220],'k--') %first image block 5
plot([81 81],[130 220],'k--') %first image block 6
hold off
xlabel('Image #')
ylabel('Fixation Duration (ms)')
legend('Novel','Repeat','Block Start')

%% fixation duration averaged over 16 images by block #
an = [];
ar = [];

for block = 1:6;
    an = [an; nanmean(nanmean(all_nov_fix_durs(:,16*(block-1)+1:block*16)))];
    ar = [ar; nanmean(nanmean(all_rep_fix_durs(:,16*(block-1)+1:block*16)))];
end

figure
hold on
plot(an)
plot(ar,'r')
hold off
xlabel('Block #')
ylabel('Fixation Duration (ms)')

%% fixation duration by image within block
anb = cell(1,16);
arb = cell(1,16);

for img = 1:16;
    anb{img} = [anb{img}; all_nov_fix_durs(:,img:16:end)];
    arb{img} = [arb{img}; all_rep_fix_durs(:,img:16:end)];
end

anbm = NaN(1,16);
arbm = NaN(1,16);
for img = 1:16
    anbm(img) = nanmean(anb{img}(1:end));
    arbm(img) = nanmean(arb{img}(1:end));
end

figure
hold on
plot(anbm)
plot(arbm,'r')
hold off
xlabel('Image #')
ylabel('Fixation Duration (ms)')


%% by Set number
set_means = [];
for set = 1:length(first_fix)
    set_means(1,set) = nanmean(first_fix{set}(1,:));
    set_means(2,set) = nanmean(first_fix{set}(2,:));
end
figure
plot(set_means')
xlabel('Data Session')
ylabel('Fixation Duration (ms)')

%% get fixation pdfs for first fixation
nov_map = zeros(600,800);
rep_map = zeros(600,800);
for set = 1:length(first_fix_location);
    for imgnum = 1:96
        if ~isnan(first_fix_location{1,set}(1,imgnum))
            xy = floor(first_fix_location{1,set}(:,imgnum));
            xy(xy < 1) = 1;
            xy(2,xy(2) > 600) = 600;
            xy(1,xy(1) > 800) = 800;
            nov_map(xy(2),xy(1)) = nov_map(xy(2),xy(1))+1;
            
        end
        
        if ~isnan(first_fix_location{2,set}(1,imgnum))
            xy = floor(first_fix_location{2,set}(:,imgnum));
            xy(xy < 1) = 1;
            xy(2,xy(2) > 600) = 600;
            xy(1,xy(1) > 800) = 800;
            rep_map(xy(2),xy(1)) = rep_map(xy(2),xy(1))+1;
        end
    end
end
