function DetermineListSQCortexKLdiverence(ListSQdatafile)
%code based on KLdivergence.m in ListSQ folder. Written/modified 3/18/15 by
%Seth Konig.
%Code calculates KLdivergence values aka similarity in fixation locations
%during novel and familiar presentation of an image.

load(ListSQdatafile)

[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);

which_image = NaN(1,115);%cnd to image number
which_image(20:35) = 1:16;
which_image(36:51) = 17:32;
which_image(52:67) = 33:48;
which_image(68:83) = 49:64;
which_image(84:99) = 65:80;
which_image(100:115) = 81:96;

all_images = NaN(2,192);%row 1 trial #, row 2 img #
img_count = 1;
for t = 1:length(fixationstats);
    if ~any(itmlist(per(t).cnd) == sequence_items) %only want image trials
        img = which_image(itmlist(per(t).cnd));
        all_images(1,img_count) = t;%trial #
        all_images(2,img_count) = img; %img #
        img_count =  img_count+1;
    end
end

KLshuff = NaN(1,5);
KLnorm = NaN(1,5);
imageX = 800;
imageY = 600;


for img = 1:max(all_images(2,:));
    
    img_ind = find(all_images(2,:) == img);
    if length(img_ind) == 2 %so novel and repeat
        %these will contain PDFs (probability distribution functions) of fixation locations
        novelfixations = cell(1,5);
        repeatfixations = cell(1,5);
        shuffled_novelfixations = cell(1,5);
        shuffled_repeatfixations = cell(1,5);
        for i = 1:5; %fill with a matrix of zeros
            novelfixations{i} = zeros(imageY,imageX);
            repeatfixations{i}=zeros(imageY,imageX);
            shuffled_novelfixations{i} = zeros(imageY,imageX);
            shuffled_repeatfixations{i} = zeros(imageY,imageX);
        end
        
        nov_fixations = fixationstats{all_images(1,img_ind(1))}.fixations; %fixation locations from novel presentation
        rep_fixations = fixationstats{all_images(1,img_ind(2))}.fixations; %fixation locations from repeat presentation
        
        %for plotting and debugging purposes
        nov_fixtimes = fixationstats{all_images(1,img_ind(1))}.fixationtimes; %fixation times from novel presentation
        rep_fixtimes = fixationstats{all_images(1,img_ind(2))}.fixationtimes; %fixation times from repeat presentation
        nov_xy = fixationstats{all_images(1,img_ind(1))}.XY; %eye trace from novel presentation
        rep_xy = fixationstats{all_images(1,img_ind(2))}.XY; %eye trace from repeat presntations
        
        %remove fixations before the crosshair
        crossfixation = find(nov_fixations(1,:) > -2.5 & nov_fixations(1,:) < 2.5 ...
            & nov_fixations(2,:) > -2.5 & nov_fixations(2,:) < 2.5);
        crossfixation = crossfixation(1);
        nov_fixations(:,1:crossfixation) = [];
        nov_fixtimes(:,1:crossfixation) = [];
        
        crossfixation = find(rep_fixations(1,:) > -2.5 & rep_fixations(1,:) < 2.5 ...
            & rep_fixations(2,:) > -2.5 & rep_fixations(2,:) < 2.5);
        crossfixation = crossfixation(1);
        rep_fixations(:,1:crossfixation) = [];
        rep_fixtimes(:,1:crossfixation) = [];
        
        nov_xy = nov_xy(:,nov_fixtimes(1):end);
        rep_xy = rep_xy(:,rep_fixtimes(1):end);
        
        %convert from dva to pixels
        nov_fixations(1,:) = round(24*nov_fixations(1,:)+400);
        nov_fixations(2,:) = round(24*nov_fixations(2,:)+300);
        rep_fixations(1,:) = round(24*rep_fixations(1,:)+400);
        rep_fixations(2,:) = round(24*rep_fixations(2,:)+300);
        
        nov_xy(1,:) = 24*nov_xy(1,:)+400;
        nov_xy(2,:) = 24*nov_xy(2,:)+300;
        rep_xy(1,:) = 24*rep_xy(1,:)+400;
        rep_xy(2,:) = 24*rep_xy(2,:)+300;
        
        %we want to take the same number of fixations from the novel and repeat trials
        maxfixations = min(size(nov_fixations,2),size(rep_fixations,2));
        
        for fixation = 1:maxfixations;
            %for all fixations 1-all
            
            nov_fix_x = nov_fixations(1,fixation);%horizontal fixation position for novel presentation
            nov_fix_y = nov_fixations(2,fixation);%vertical fixation position for novel presentation
            rep_fix_x = rep_fixations(1,fixation);%horizonal fixation position for repeat presentation
            rep_fix_y = rep_fixations(2,fixation);%vertical fixation position for repeat presentation
            shuff_nov_fix_x = randi(imageX);
            shuff_nov_fix_y = randi(imageY);
            shuff_rep_fix_x = randi(imageX);
            shuff_rep_fix_y = randi(imageY);
            
            %make sure fixations are within image borders
            nov_fix_x(nov_fix_x < 1) = 1;
            nov_fix_x(nov_fix_x > imageX) = imageX;
            nov_fix_y(nov_fix_y < 1) = 1;
            nov_fix_y(nov_fix_y > imageY) = imageY;
            rep_fix_x(rep_fix_x < 1) = 1;
            rep_fix_x(rep_fix_x > imageX) = imageX;
            rep_fix_y(rep_fix_y < 1) = 1;
            rep_fix_y(rep_fix_y > imageY) = imageY;
            %shuffled may be 0 but not greater than imageX or imageY
            shuff_nov_fix_x(shuff_nov_fix_x < 1) = 1;
            shuff_nov_fix_y(shuff_nov_fix_y < 1) = 1;
            shuff_rep_fix_x(shuff_rep_fix_x < 1) = 1;
            shuff_rep_fix_y(shuff_rep_fix_y < 1) = 1;
            
            novelfixations{5}(nov_fix_y,nov_fix_x) = novelfixations{5}(nov_fix_y,nov_fix_x)+1;
            repeatfixations{5}(rep_fix_y,nov_fix_x) = repeatfixations{5}(rep_fix_y,nov_fix_x)+1;
            shuffled_novelfixations{5}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{5}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
            shuffled_repeatfixations{5}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{5}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
        end
        
        maxfixations(maxfixations > 20) = 20; %don't care if there are more than 20 fixations
        %since puting into groups of five remove the remainder of #/5,
        %median is 17-19 for TT and PW for novel and repeats
        maxfixations = maxfixations-rem(maxfixations,5);
        
        if maxfixations >=5
            for fixation = 1:maxfixations;
                nov_fix_x = nov_fixations(1,fixation);%horizontal fixation position for novel presentation
                nov_fix_y = nov_fixations(2,fixation);%vertical fixation position for novel presentation
                rep_fix_x = rep_fixations(1,fixation);%horizonal fixation position for repeat presentation
                rep_fix_y = rep_fixations(2,fixation);%vertical fixation position for repeat presentation
                shuff_nov_fix_x = randi(imageX);
                shuff_nov_fix_y = randi(imageY);
                shuff_rep_fix_x = randi(imageX);
                shuff_rep_fix_y = randi(imageY);
                
                %make sure fixations are within image borders
                nov_fix_x(nov_fix_x < 1) = 1;
                nov_fix_x(nov_fix_x > imageX) = imageX;
                nov_fix_y(nov_fix_y < 1) = 1;
                nov_fix_y(nov_fix_y > imageY) = imageY;
                rep_fix_x(rep_fix_x < 1) = 1;
                rep_fix_x(rep_fix_x > imageX) = imageX;
                rep_fix_y(rep_fix_y < 1) = 1;
                rep_fix_y(rep_fix_y > imageY) = imageY;
                %shuffled may be 0 but not greater than imageX or imageY
                shuff_nov_fix_x(shuff_nov_fix_x < 1) = 1;
                shuff_nov_fix_y(shuff_nov_fix_y < 1) = 1;
                shuff_rep_fix_x(shuff_rep_fix_x < 1) = 1;
                shuff_rep_fix_y(shuff_rep_fix_y < 1) = 1;
                
                %put fixations in their appropriate PDF. Mark matrix with a
                %1 where there was a fixation
                if fixation <= 5
                    novelfixations{1}(nov_fix_y,nov_fix_x) = novelfixations{1}(nov_fix_y,nov_fix_x)+1;
                    repeatfixations{1}(rep_fix_y,nov_fix_x) = repeatfixations{1}(rep_fix_y,nov_fix_x)+1;
                    shuffled_novelfixations{1}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{1}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                    shuffled_repeatfixations{1}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{1}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                elseif fixation <= 10
                    novelfixations{2}(nov_fix_y,nov_fix_x) = novelfixations{2}(nov_fix_y,nov_fix_x)+1;
                    repeatfixations{2}(rep_fix_y,nov_fix_x) = repeatfixations{2}(rep_fix_y,nov_fix_x)+1;
                    shuffled_novelfixations{2}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{2}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                    shuffled_repeatfixations{2}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{2}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                elseif fixation <=15
                    novelfixations{3}(nov_fix_y,nov_fix_x) = novelfixations{3}(nov_fix_y,nov_fix_x)+1;
                    repeatfixations{3}(rep_fix_y,nov_fix_x) = repeatfixations{3}(rep_fix_y,nov_fix_x)+1;
                    shuffled_novelfixations{3}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{3}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                    shuffled_repeatfixations{3}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{3}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                elseif fixation <=20
                    novelfixations{4}(nov_fix_y,nov_fix_x) = novelfixations{4}(nov_fix_y,nov_fix_x)+1;
                    repeatfixations{4}(rep_fix_y,nov_fix_x) = repeatfixations{4}(rep_fix_y,nov_fix_x)+1;
                    shuffled_novelfixations{4}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{4}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                    shuffled_repeatfixations{4}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{4}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                end
            end
        end
        for i = 1:5
            if i == 5
                if sum(sum(novelfixations{i})) >= 5 &&  sum(sum(novelfixations{i})) >= 5
                    Distance = KL_Divergence(novelfixations{i},repeatfixations{i});
                    KLnorm(img,i) = Distance;
                    Distance = KL_Divergence(shuffled_novelfixations{i},shuffled_repeatfixations{i});
                    KLshuff(img,i) = Distance;
                else
                    KLnorm(img,i) = NaN;
                    KLshuff(img,i) = NaN;
                end
            else
                
                if sum(sum(novelfixations{i})) == 5 &&  sum(sum(novelfixations{i})) == 5
                    Distance = KL_Divergence(novelfixations{i},repeatfixations{i});
                    KLnorm(img,i) = Distance;
                    Distance = KL_Divergence(shuffled_novelfixations{i},shuffled_repeatfixations{i});
                    KLshuff(img,i) = Distance;
                else
                    KLnorm(img,i) = NaN;
                    KLshuff(img,i) = NaN;
                end
            end
        end
    end
end

save([ListSQdatafile(1:end-13) '-KLdivergence.mat'],'KLshuff','KLnorm')
