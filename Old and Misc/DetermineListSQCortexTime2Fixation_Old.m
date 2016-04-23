function DetermineListSQCortexTime2Fixation(ListSQdatafile)
%written by Seth Konig December, 2014
%code determines how long it takes to fixation an item in the sequence
%portion of the ListSQ task. For cortex/behavioral only days/
load(ListSQdatafile)

[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);

num_trials = length(fixationstats);
fixwin = 5; 
reaction_time = NaN(length(fixationstats),4); %how long after a stimulus disappears do they saccade
time_to_fixation = NaN(length(fixationstats),4); %how fast did they get to it the first time around
fixation_accuracy = NaN(length(fixationstats),4); %how far off
fixation_duration = NaN(length(fixationstats),4); %fixation duration
fixation_on_test = NaN(length(fixationstats),4); %how long the stimulus was on while fixating
extrafixations = NaN(length(fixationstats),4); %how many noncross hair fixations do they make
for trial = 1:num_trials
    if any(itmlist(per(trial).cnd) == sequence_items)
        % x and y data needs to be converted back into pixels
        fixations = 24*fixationstats{trial}.fixations;
        fixations(1,:) = fixations(1,:) + 400;
        fixations(2,:) = fixations(2,:) + 300;
        fixationtimes = fixationstats{trial}.fixationtimes;
        saccadetimes = fixationstats{trial}.saccadetimes;
        xy = 24*fixationstats{trial}.XY;
        xy(1,:) = xy(1,:)+400;
        xy(2,:) = xy(2,:)+300;
        
        % Find saccades that were small in amplitude and determine if these
        % smaller saccades were corrective and I should essentially count
        % the pre and post saccade fixations as both on the cross hair.
        saccade_amplitudes = NaN(1,size(saccadetimes,2));
        for sacc = 1:size(saccadetimes,2);
            sacx = xy(1,saccadetimes(2,sacc))-xy(1,saccadetimes(1,sacc));
            sacy = xy(2,saccadetimes(2,sacc))-xy(2,saccadetimes(1,sacc));
            saccade_amplitudes(sacc) = sqrt(sacx.^2+sacy.^2);
        end
        if any(saccade_amplitudes(2:end) <= fixwin/2) %going to generalize  these as corrective saccades. Ignore 1st saccade
            [fixationtimes,fixations,saccadetimes] = remove_corrective_saccades(...
                xy,saccadetimes,saccade_amplitudes,fixwin);
        end
        
        if itmlist(per(trial).cnd) == sequence_items(1)
            locs = sequence_locations{1};
        else
           locs = sequence_locations{2};
        end
           
        
%         
%             figure
%             hold on
%             plot(xy(1,:),xy(2,:),'g');
%             %         for f = 1:length(fixationtimes);
%             %             plot(fixations(1,f),fixations(2,f),'k*','markersize',5);
%             %         end
%             for c = 1:size(locs,2);
%                 plot(locs(1,c),locs(2,c),'kx','markersize',12)
%             end
%             xlim([0 800])
%             ylim([0 600])
        
        valid_fixations = ones(1,size(fixations,2));
        for c = 1:size(locs,2)
            if c == 1
                valid_window = [per(trial).test(c,1) per(trial).test(c,2)];
            else
                valid_window = [per(trial).test(c-1,2) per(trial).test(c,2)];
            end
            valid_ind = valid_window(1):valid_window(2);
            valid_ind(valid_ind > length(xy)) = [];
            
            potential_fix = [];
            for f = 1:size(fixationtimes,2);
                if valid_fixations(f)
                    fixtimes = fixationtimes(1,f):fixationtimes(2,f);
                    C = intersect(fixtimes,valid_ind);
                    if length(C) >= 50;
                        potential_fix = [potential_fix f];
                    end
                end
            end
            if ~isempty(potential_fix)
                if length(potential_fix)  == 1;
                    fixstart = fixationtimes(1,potential_fix);
                    valid_fixations(1:potential_fix) = 0;
                    dist = sqrt((fixations(1,potential_fix)-locs(1,c)).^2+...
                        (fixations(2,potential_fix)-locs(2,c)).^2);
                elseif length(potential_fix) > 1;
                    dist = sqrt((fixations(1,potential_fix)-locs(1,c)).^2+...
                        (fixations(2,potential_fix)-locs(2,c)).^2);
                    [~,thefix] = min(dist);
                    fixstart = fixationtimes(1,potential_fix(thefix));
                    valid_fixations(1:potential_fix(thefix)) = 0;
                    potential_fix = potential_fix(thefix);
                    dist = dist(thefix);
                    if thefix ~= 1
                        extrafixations(trial,c) = thefix-1;
                    end
                end
                
                time_to_fixation(trial,c) = fixstart-per(trial).test(c,1);
                
                fixation_duration(trial,c) = fixationtimes(2,potential_fix)-fixationtimes(1,potential_fix)+1;
                if c < size(locs,2)
                    reaction_time(trial,c) = fixationtimes(2,potential_fix)-per(trial).test(c,2);
                end
                
                timeON = per(trial).test(c,1):per(trial).test(c,2);
                fixON = fixationtimes(1,potential_fix):fixationtimes(2,potential_fix);
                fixation_on_test(trial,c) = length(intersect(timeON,fixON));
                
                fixation_accuracy(trial,c) = dist;
                
                
%                             plot(xy(1,fixationtimes(1,potential_fix):fixationtimes(2,potential_fix)),...
%                                 xy(2,fixationtimes(1,potential_fix):fixationtimes(2,potential_fix)),'c');
%                             plot(xy(1,fixstart),xy(2,fixstart),'m*','markersize',10)
                
            end
        end
%         
%             pause(0.5)
%             close all
    end
end


save([ListSQdatafile(1:end-13) '-SQRTs.mat'],'time_to_fixation','reaction_time',...
   'fixation_accuracy','fixation_duration','fixation_on_test','extrafixations')

end


function [fixationtimes,fixations,saccadetimes] = remove_corrective_saccades(...
    xy,saccadetimes,saccade_amplitudes,fixwin)
% written by Seth Konig on June 26, 2014
% function removes "corrective saccades" determine as saccades (other than
% the 1st recorded saccade) that are shorter than a defined threshold.
% This threshold is set by the "saccade_amplitude" variable.

corrective_saccades = find(saccade_amplitudes <= fixwin/2);
corrective_saccades(corrective_saccades == 1) = []; %ignore these
saccadetimes(:,corrective_saccades) =[];%remove corrective saccades

saccade_indexes = [];
for sac = 1:size(saccadetimes,2)
    saccade_indexes = [saccade_indexes saccadetimes(1,sac):saccadetimes(2,sac)];
end

fixation_indexes = 1:size(xy,2);
[~ , ia, ~] = intersect(fixation_indexes,saccade_indexes);
fixation_indexes(ia) = [];

[fixationtimes,fixations] = BehavioralIndexXY(fixation_indexes,xy(1,:),xy(2,:));
end

function [behaviortime, behaviormean] = BehavioralIndexXY(behavind,x,y)
%function is the same as above but also calculates mean fixation position
dind = diff(behavind);
gaps =find(dind > 1);
behaveind = zeros(length(gaps),50);
if ~isempty(gaps)
    for gapind = 1:length(gaps)+1;
        if gapind == 1;
            temp = behavind(1:gaps(gapind));
        elseif gapind == length(gaps)+1
            temp = behavind(gaps(gapind-1)+1:end);
        else
            temp = behavind(gaps(gapind-1)+1:gaps(gapind));
        end
        behaveind(gapind,1:length(temp)) = temp;
    end
else
    behaveind =  behavind;
end
behaviortime = zeros(2,size(behaveind,1));
behaviormean = zeros(2,size(behaveind,1));
for index=1:size(behaveind,1)
    rowfixind = behaveind(index,:);
    rowfixind(rowfixind == 0) = [];
    behaviortime(:,index) = [rowfixind(1);rowfixind(end)];
    behaviormean(:,index) = [mean(x(rowfixind));...
        mean(y(rowfixind))];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for plotting data on fly keep commented out

% dd = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Cortex Data\';
% 
% a = what(dd);
% a = a.mat;
% 
% s = [];
% all_t2f = [];%cell(1,5);
% for f = 1:size(a,1);
%     if ~isempty(strfind(a{f,:},'SQRTs.mat'))
%         load([dd a{f,:}],'time_to_fixation')
%         all_t2f = [all_t2f; time_to_fixation];   
% 
% %         all_t2f{f} = time_to_fixation;
%     end
% end
% 
% %%
% count= [];
% for f = 1:size(all_t2f,2);
%     if ~isempty(all_t2f{f})
%         pct = sum(all_t2f{f} < 50)./sum(~isnan(all_t2f{f}))*100;
%         pct(pct <10 ) = [];
%         count = [count pct];
%     end
% end