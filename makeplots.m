function x = makeplots(tt,SmallStepsI,hI,hII)
%PLots full wavefunctions plus narrow regions around t0
figure;
plot(tt,hI);
hold on;
figure;
plot(tt,hII);
hold on;
narrow_range = 3000;
tplot_narrow = tt(SmallStepsI - narrow_range: SmallStepsI + narrow_range);
hI_narrow = hI(SmallStepsI - narrow_range: SmallStepsI + narrow_range);
figure;
plot(tplot_narrow,hI_narrow);
hold on;
figure;
hII_narrow = hII(SmallStepsI - narrow_range: SmallStepsI + narrow_range);
plot(tplot_narrow,hII_narrow);
hold on;
%older version:
%midpt = SmallStepsI - dN1;
%narrow_range = 1000;
%tplot_narrow = tplot(midpt - narrow_range: midpt + narrow_range);
%hI_narrow = hI(midpt - narrow_range: midpt+narrow_range);
%figure;
%plot(tplot_narrow,hI_narrow);
%hold on;
%figure;
%hII_narrow = hII(midpt - narrow_range: midpt+narrow_range);
%plot(tplot_narrow,hII_narrow);
%plot(tplot,hII);
x=0;

