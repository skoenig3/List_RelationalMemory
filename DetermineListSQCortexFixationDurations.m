function DetermineListSQCortexFixationDurations(ListSQdatafile)
% written by Seth Konig December, 2014
% function imports ListSQ cortex data and determines fixation durations for
% novel and repeated stimuli
% modified 3/18/2015 by Seth Konig to also calculate saccade amplitudes


load(ListSQdatafile)

fixation_durations{1} = NaN(96,40);
fixation_durations{2} = NaN(96,40);
saccade_amplitudes{1} = NaN(96,40); 
saccade_amplitudes{2} = NaN(96,40); 

[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);

which_image = NaN(1,115);%cnd to image number
which_image(20:35) = 1:16;
which_image(36:51) = 17:32;
which_image(52:67) = 33:48;
which_image(68:83) = 49:64;
which_image(84:99) = 65:80;
which_image(100:115) = 81:96;

all_images = ones(1,96);
for t = 1:length(fixationstats);
    if ~any(itmlist(per(t).cnd) == sequence_items) %only want image trials
        img = which_image(itmlist(per(t).cnd));
        if all_images(img) == 1%novel image
            all_images(img) = 0;
            novrow = 1;
        else %repeated image
            novrow = 2;
        end
            
        fixationtimes = fixationstats{t}.fixationtimes;
        fixations = fixationstats{t}.fixations;
        xy = fixationstats{t}.XY;
        saccadetimes = fixationstats{t}.saccadetimes; 
        if ~isempty(fixationtimes)
            crossfixation = find(fixations(1,:) > -2.5 & fixations(1,:) < 2.5 ...
                & fixations(2,:) > -2.5 & fixations(2,:) < 2.5);
            crossfixation = crossfixation(1);
            fixationtimes(:,1:crossfixation) = [];
            
            if fixationtimes(1,1) < saccadetimes(1,1);%if fixations occur before saccades
                saccadetimes(:,1:crossfixation-1) = [];
            else
                saccadetimes(:,1:crossfixation) = [];
            end

            fixdurs = diff(fixationtimes)+1;
            if length(fixdurs) > 40
                fixdurs = fixdurs(1:40);
            end
            
            sacamps = [];
            for s = 1:size(saccadetimes,2);
                sacx = xy(1,saccadetimes(2,s))-xy(1,saccadetimes(1,s));
                sacy = xy(2,saccadetimes(2,s))-xy(2,saccadetimes(1,s));
                sacamps(s) = sqrt(sacx^2+sacy^2);
            end
            
            fixation_durations{novrow}(img,1:length(fixdurs))=fixdurs;
            saccade_amplitudes{novrow}(img,1:length(sacamps)) = sacamps; 
            
        end
    end
end

save([ListSQdatafile(1:end-13) '-fixdurs.mat'],'fixation_durations','saccade_amplitudes')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not official code but for plotting on the fly
% leave commented out

% dd = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Cortex Data\';
% a = what(dd);
% a = a.mat; 
% 
% all_fixdurs = cell(1,2);
% for f = 1:size(a,1);
%    if ~isempty(strfind(a(f,:),'fixdurs'))
%        load(a{f,:})
%        all_fixdurs{1} = [all_fixdurs{1}; fixation_durations{1}(:,1:20)];
%        all_fixdurs{2} = [all_fixdurs{2}; fixation_durations{2}(:,1:20)];
%    end
% end
% %%
% figure
% hold on
% plot(nanmean(all_fixdurs{1}(:,1:20)),'blue')
% errorb(nanmean(all_fixdurs{1}(:,1:20)),...
%     nanstd(all_fixdurs{1}(:,1:20))./sqrt(sum(~isnan(all_fixdurs{1}(:,1:20)))),'color','blue')
% plot(nanmean(all_fixdurs{2}(:,1:20)),'red')
% errorb(nanmean(all_fixdurs{2}(:,1:20)),...
%     nanstd(all_fixdurs{2}(:,1:20))./sqrt(sum(~isnan(all_fixdurs{2}(:,1:20)))),'color','red')
% xlabel('Ordinal Fixation #')
% ylabel('Fixtaion Duration (ms)')
% legend('Novel','Repeat')