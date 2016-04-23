% written by Seth Konig December 2014
% code runs all the other code preprocess then process behavioral only data 
% when task was run on cortex and not during a recording

datadir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ_RelationalMemory\Cortex Data\';

cch25_files = {'TT141020.1','TT141024.2','TT141027.2','TT141028.2','TT141104.2'};
listsq_files = {'TT141020.2','TT141024.1','TT141027.1','TT141028.1','TT141104.1'};
item_sets =  {'ListSQ07.itm','ListSQ09.itm','ListSQ10.itm','ListSQ11.itm','ListSQ13.itm'};

cch25_files = {'PW140805.2','PW140806.2','PW140829.1','PW140910.1','PW140915.1',...
               'PW140917.1','PW141007.1','PW141008.1','PW141010.1','PW141013.1',...
               'PW141014.1','PW141015.1','PW141016.1'};
listsq_files = {'PW140805.3','PW140806.3','PW140829.3','PW140910.3','PW140915.3',...
                'PW140917.3','PW141007.3','PW141008.3','PW141010.3','PW141013.3',...
                'PW141014.3','PW141015.3','PW141016.3'};
item_sets = {'ListSQ07.itm','ListSQ08.itm','ListSQ12.itm','ListSQ14.itm','ListSQ15.itm',...
             'ListSQ16.itm','ListSQ20.itm','ListSQ21.itm','ListSQ22.itm','ListSQ23.itm',...
             'ListSQ24.itm','ListSQ25.itm','ListSQ26.itm'};

% cch25_files = {'PW150306.1','PW150311.1','PW150313.1','PW150316.1'};
% listsq_files = {'PW150306.2','PW150311.2','PW150313.2','PW150316.2'};
% item_sets =  {'ListSQVR.itm','ListSQVS.itm','ListSQVT.itm','ListSQVU.itm'};


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %---Preprocess all the ListSQ data---%%%
% for ls =1:length(listsq_files)
%     ImportListSqData(cch25_files{ls},listsq_files{ls},item_sets{ls})
%     close all
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ---Autmoatically analyze Reaction times in Sequence Task---%%%
% for ls =1:length(listsq_files)
%     DetermineListSQCortexTime2Fixation([datadir listsq_files{ls}(1:8) '_'  listsq_files{ls}(end) '-fixation.mat'])
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ---Autmoatically analyze Fixation Durations in List task---%%%
% for ls =1:length(listsq_files)
%     DetermineListSQCortexFixationDurations([datadir listsq_files{ls}(1:8) '_'  listsq_files{ls}(end) '-fixation.mat'])
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%---Automatically analyze KL Divergence Data for List task---%%%
% for ls = 1:length(listsq_files)
%    DetermineListSQCortexKLdiverence([datadir listsq_files{ls}(1:8) '_'  listsq_files{ls}(end) '-fixation.mat']);
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot mean reaction time across sessions---%%%
block_prop = NaN(length(listsq_files),21);
means = NaN(1,4);
for ls =1:length(listsq_files)
    load([datadir listsq_files{ls}(1:8) '_' listsq_files{ls}(end) '-SQRTs.mat'],'reaction_time');
    means(ls,:) =100*sum(reaction_time(21:end,:) < 150)./sum(~isnan(reaction_time(21:end,:)));
    nb = floor(size(reaction_time,1)/20);
    for n = 1:nb
        temp = reaction_time(20*(n-1)+1:20*n,2:4);
        temp = temp(1:end);
        block_prop(ls,n) = sum(temp < 150)./sum(~isnan(temp));
    end
end

figure
hold on
bar(nanmean(means));
errorb(nanmean(means),nanstd(means)./sqrt(sum(~isnan(means))))
set(gca,'Xtick',[1:4])
xlabel('Item Number')
ylabel('Mean Reacion times (time to fixation) (ms)')
ylabel('Percentageof Reaction times < 150 ms')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Plot mean Fixation Durations and Saccade Amplitudes Across Sets---%%%
% sac_means =cell(1,2); 
% means = cell(1,2);
% n{1} = [];
% n{2} = [];
% for ls =1:length(listsq_files)
%     load([datadir listsq_files{ls}(1:8) '_' listsq_files{ls}(end) '-fixdurs.mat'],...
%         'fixation_durations','saccade_amplitudes');
%     means{1} = [means{1}; nanmean(fixation_durations{1}(:,2:21))];
%     means{2} = [means{2}; nanmean(fixation_durations{2}(:,2:21))];
%     n{1} = [n{1} sum(~isnan(fixation_durations{1}'))];
%     n{2} = [n{2} sum(~isnan(fixation_durations{2}'))];
%     sac_means{1} = [sac_means{1}; nanmean(saccade_amplitudes{1}(:,1:20))];
%     sac_means{2} = [sac_means{2}; nanmean(saccade_amplitudes{2}(:,1:20))];
% end
% figure
% hold on
% errorbar(nanmean(means{1}),nanstd(means{1})./sqrt(sum(~isnan(means{1}))),'b')
% errorbar(nanmean(means{2}),nanstd(means{2})./sqrt(sum(~isnan(means{2}))),'r')
% xlabel('Ordinal Fixation #')
% ylabel('Fixation Duration (ms)')
% legend('Novel','Repeat')
% 
% nn_90 = [];
% nn_80 = [];
% nr_90 = [];
% nr_80 = [];
% p = [];
% for f = 1:20;
%     nn_80(f) = sampsizepwr('t',[mean(means{1}(:,f)) std(means{1}(:,f))], mean(means{2}(:,f)),0.8);
%     nn_90(f) = sampsizepwr('t',[mean(means{1}(:,f)) std(means{1}(:,f))], mean(means{2}(:,f)),0.9);
%     nr_80(f) = sampsizepwr('t',[mean(means{2}(:,f)) std(means{2}(:,f))], mean(means{1}(:,f)),0.8);
%     nr_90(f) = sampsizepwr('t',[mean(means{2}(:,f)) std(means{2}(:,f))], mean(means{1}(:,f)),0.9);
%    [~,p(f)] = ttest2( means{1}(:,f),means{2}(:,f));
% end
% 
% figure
% hold on
% errorbar(nanmean(sac_means{1}),nanstd(sac_means{1})./sqrt(sum(~isnan(sac_means{1}))),'b')
% errorbar(nanmean(sac_means{2}),nanstd(sac_means{2})./sqrt(sum(~isnan(sac_means{2}))),'r')
% xlabel('Ordinal Sacccade #')
% ylabel('Saccade Amplitude (dva)')
% legend('Novel','Repeat')
%% Plot KL divergence 

KLnorm_mean = NaN(length(listsq_files),5);
KLshuff_mean = NaN(length(listsq_files),5);
means = NaN(1,4);
for ls =1:length(listsq_files)
   load([datadir listsq_files{ls}(1:8) '_' listsq_files{ls}(end) '-KLdivergence.mat'],'KLnorm','KLshuff');
   KLnorm_mean(ls,:) = nanmean(KLnorm);
   KLshuff_mean(ls,:) = nanmean(KLshuff);
end

figure
plot(nanmean(KLnorm_mean))
hold on
plot(nanmean(KLshuff_mean),'r')