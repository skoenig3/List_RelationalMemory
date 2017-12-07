x = 1:27;
nov = nanmean(all_novel_fix_dur);
y = nov;
myfit = fittype('a + b*log(x)',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','b'});
[FO, G]  = fit(x',y',myfit)
%%
yp = FO.a+FO.b*log(x);

