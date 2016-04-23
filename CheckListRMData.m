function CheckListRMData(data_dir,fixationfile)
%written by Seth Konig 5/8/15
%code displays the number of images viewed and the number of images viewed
%by the amount of time the images were displayed for/how much the monkey
%looked away

load([data_dir fixationfile]);

disp('______________________________')
disp(['Data for ' fixationfile(1:10)])

viewed = zeros(2,96);
img_dur = zeros(2,96);
for trial = 1:length(per);
    [r,c] = find(imgs == per(trial).cnd-1000);
    if viewed(r,c) ~= 0
        error('This should be the first viewing of this image, but it is not!')
    else
        viewed(r,c) = 1;
        img_dur(r,c) = per(trial).alltim(per(trial).allval == 24)-...
            per(trial).alltim(per(trial).allval == 23);
    end
end

%1st image trial was displayed with the wrong timing file, then when fixed
%started the task over again. The first image therefore is not novel and thus
%we will not be analyzing data for this image on this session. 
if strcmpi('PW150416_2',fixationfile(1:10))
    viewed(:,1) = 0;
    imgdur(:,1) = 0;
end

total_images_viewed = sum(sum(viewed));
disp(['Viewed ' num2str(total_images_viewed) ' images'])

total_nov_repeat_pairs = length(find(viewed(1,:) == viewed(2,:))); 
disp(['Viewed ' num2str(total_nov_repeat_pairs) ' novel and repeat pairs'])

great_views = NaN(2,96);
great_views(viewed & img_dur < 8000) = 1; %1000 ms buffer for accidental lookaways and cortex timing

total__great_nov_repeat_pairs = length(find(great_views(1,:) == great_views(2,:))); 
disp([num2str(total__great_nov_repeat_pairs) ' novel and repeat pairs with < 1 sec of looking away (great)'])

ok_views = NaN(2,96);
ok_views(viewed & img_dur < 10000) = 1; %looked away a bit 

total__ok_nov_repeat_pairs = length(find(ok_views(1,:) == ok_views(2,:))); 
disp([num2str(total__ok_nov_repeat_pairs) ' novel and repeat pairs with < 3 secs of looking away (ok)'])