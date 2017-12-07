%repeated measures test of ListRM fixation durations
%run ListRM_pre_vs_post.m first
%%


subject_id = [ones(17,1); 2*ones(20,1); ones(15,1); 2*ones(20,1);...
             ones(17,1); 2*ones(20,1); ones(15,1); 2*ones(20,1)];

pre_post_id = [ones(17,1); 2*ones(15,1); ones(17,1); 2*ones(15,1);];

nov_rep_id = [ones(17,1); ones(15,1);
             2*ones(17,1); 2*ones(15,1)];

data = [mean(pre_novel_fix_dur(:,2:10),2); mean(post_novel_fix_dur(:,2:10),2);...
        mean(pre_repeat_fix_dur(:,2:10),2); mean(post_repeat_fix_dur(:,2:10),2)];

%%
rmSTATS =rm_anova2(data,subject_id,nov_rep_id,pre_post_id,{'Nov/Rep','Pre/Post'})
STATS =anovan(data,[nov_rep_id,pre_post_id],'model','interaction',...
    'varnames',strvcat('Nov/Rep', 'Pre/Post'))
%%
%% Levene's test for equal variance
p_pre_post = vartestn(data,[pre_post_id],'TestType','LeveneAbsolute')
p_nov_rep = vartestn(data,[nov_rep_id],'TestType','LeveneAbsolute')
%%
%% Create table
Factor = rmSTATS(2:end,1);
Sum_Sq = cell2mat(rmSTATS(2:end,2));
df = cell2mat(rmSTATS(2:end,3));
Mean_Sq = cell2mat(rmSTATS(2:end,4));
F = [cell2mat(rmSTATS(2:end,5)); NaN(3,1)];
p = [cell2mat(rmSTATS(2:end,6)); NaN(3,1)];
T = table(Sum_Sq,df,Mean_Sq,F,p,'RowNames',Factor);

uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);