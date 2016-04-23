function combinedSimilarity(data_dir,cortex_files)

all_KLdivergence = [];
all_shuff_KLdivergence = [];
all_shuff_KLdivergence2 = [];
all_ObservedSimilarity = [];
all_ShuffledSimilarity = [];
all_ShuffledSimilarity2 = [];
for file = 1:length(cortex_files)
    load([data_dir cortex_files{file}(1:8) '_' cortex_files{file}(end) '-KLdivergence.mat']);
    all_KLdivergence = [all_KLdivergence; nanmean(KLnorm)];
    all_shuff_KLdivergence = [all_shuff_KLdivergence; nanmean(KLshuff)];
    all_shuff_KLdivergence2 = [all_shuff_KLdivergence2; nanmean(KLshuff2)];
    all_ObservedSimilarity(file) = 100*nanmean(ObservedSimilarity)/(600*800);
    all_ShuffledSimilarity(file) = 100*nanmean(ShuffledSimilarity/(600*800));
    all_ShuffledSimilarity2(file) = 100*nanmean(ShuffledSimilarity2/(600*800));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% KL divergence stuff %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---For random chance levels---%
%chance levels for fixations in groups of 5
chance = all_shuff_KLdivergence(:,1:5);
chance = chance(1:end);

%chance levels for all fixations ungrouped
chance12 = all_shuff_KLdivergence(:,6);
chance12 = chance12(1:end);

cIs = [];
pvalus = [];
for fix = 1:5
    [~,p,cI] = ztest(chance,mean(all_KLdivergence(:,fix)),std(chance));
    cIs(fix) = cI(1);
    pvalues(fix) = p;
end

[~,pvalues2,cI2] = ztest(chance12,mean(all_KLdivergence(:,6)),std(chance12));

%---For chance levels calculated by random pairs of novel and repeat images--%
pvaluesp = [];
for fix = 1:6
    [~,p] = ttest(all_shuff_KLdivergence2(:,fix),all_KLdivergence(:,fix));
    pvaluesp(fix) = p;
end


figure
subplot(1,3,[1:2])
bar(mean(all_KLdivergence(:,1:5)))
hold on
errorb(mean(all_KLdivergence(:,1:5)),std(all_KLdivergence(:,1:5))...
    ./sqrt(size(all_KLdivergence,2)))
plot([0.5 5.5],[cIs(1) cIs(1)],'k--')
for fix = 1:5
   if pvalues(fix) < 0.05
      plot(fix,cI(1)+1,'k*')
   end
end
errorbar(mean(all_shuff_KLdivergence2(:,1:5)),std(all_shuff_KLdivergence2(:,1:5))./...
    sqrt(size(all_shuff_KLdivergence2,1)),'r')
for fix = 1:5
   if pvaluesp(fix) < 0.05
      plot(fix,mean(all_shuff_KLdivergence2(:,fix))+std(all_shuff_KLdivergence2(:,fix)),'r*')
   end
end
set(gca,'Xtick',[1:5])
set(gca,'XtickLabel',{'1-5','6-10','11-15','15-20','21-25'})
xlabel('Fixation #s')
ylabel('KL divergence (Bits)')

subplot(1,3,3)
bar(mean(all_KLdivergence(:,6)))
hold on
errorb(mean(all_KLdivergence(:,6)),std(all_KLdivergence(:,6))...
    ./sqrt(size(all_KLdivergence,2)))
plot([0.5 1.5],[cI2(1) cI2(1)],'k--')
if pvalues2 < 0.05
      plot(1,cI2(1)+1,'k*')
end
errorbar(mean(all_shuff_KLdivergence2(:,6)),std(all_shuff_KLdivergence2(:,6))./...
    sqrt(size(all_shuff_KLdivergence2,1)),'r')
   if pvaluesp(6) < 0.05
      plot(1,mean(all_shuff_KLdivergence2(:,6))+std(all_shuff_KLdivergence2(:,6)),'r*')
   end
set(gca,'Xtick',[])
xlabel('All fixations')
subtitle(cortex_files{1}(1:2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spatial Similarity Stuff %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,p_sim]= ttest2(all_ObservedSimilarity,all_ShuffledSimilarity);
[~,p_sim2]= ttest2(all_ObservedSimilarity,all_ShuffledSimilarity2);

figure
hold on
bar([mean(all_ObservedSimilarity) mean(all_ShuffledSimilarity)  mean(all_ShuffledSimilarity2)])
errorb([mean(all_ObservedSimilarity) mean(all_ShuffledSimilarity) mean(all_ShuffledSimilarity2)],...
    [std(all_ObservedSimilarity) std(all_ShuffledSimilarity)  std(all_ShuffledSimilarity2)]...
    ./sqrt(length(all_ObservedSimilarity)))
if p_sim < 0.05
    plot([2],mean(all_ShuffledSimilarity)+2*std(all_ShuffledSimilarity),'k*')
end
if p_sim2 < 0.05
    plot([3],mean(all_ShuffledSimilarity2)+2*std(all_ShuffledSimilarity2),'r*')
end
hold off
set(gca,'Xtick',[1:3])
set(gca,'XtickLabel',{'Observed','Shuffled','Shuffled Paired'})
hold off
ylabel('Normalized Similarity (a.u.)')
title(cortex_files{1}(1:2))