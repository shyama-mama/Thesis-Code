function [initContactStartInd,initContactEndInd,slipStartInd,slipEndInd] = getTrackingStartAndEndInds(synchFx, synchFy, synchFz, synchContactPixelsArray)

FN = -synchFz';
FT = sqrt(synchFx'.^2 + synchFy'.^2);
FTonFN = FT./FN;

figure; set(gcf,'Color',[1 1 1]);
nPlot = 4;
ax(1) = subplot(nPlot,1,1);
plot(synchContactPixelsArray,'k');
xlabel('Frame Number'); ylabel('Pixels');
ax(2) = subplot(nPlot,1,2);
plot(FT,'k');
xlabel('Frame Number'); ylabel('F_T (N)');
ax(3) = subplot(nPlot,1,3);
plot(FN,'k');
xlabel('Frame Number'); ylabel('F_N (N)');
ax(4) = subplot(nPlot,1,4);
plot(FTonFN,'k');
xlabel('Frame Number'); ylabel('F_T/F_N');
linkaxes(ax,'x');



% User input to select the start and end of initial contact, and slip
minLen = min(length(synchContactPixelsArray),length(FT));
temp = [synchContactPixelsArray(1:minLen), FT(1:minLen), FN(1:minLen), FTonFN(1:minLen)]';
subplot(nPlot,1,1); hold on;
legendStrings = {'Contact Area Size'};
xlim([0 length(synchContactPixelsArray)]);
legend(legendStrings);
title('Please, select the start of initial contact');
[x,y] = ginput(1);
initContactStartInd = round(x);
for i = 1:nPlot
    subplot(nPlot,1,i); hold on;
    plot([initContactStartInd initContactStartInd],[0,temp(i,initContactStartInd)],'r-');
end
legendStrings = [legendStrings, {'Start Init Contact'}];
subplot(nPlot,1,1); hold on;
legend(legendStrings);
xlim([0 length(synchContactPixelsArray)]);
title('Please, select the end of initial contact');
[x,y] = ginput(1);
initContactEndInd = round(x);
for i = 1:nPlot
    subplot(nPlot,1,i); hold on;
    plot([initContactEndInd initContactEndInd],[0,temp(i,initContactEndInd)],'y-');
end
legendStrings = [legendStrings, {'End Init Contact'}];
subplot(nPlot,1,1); hold on;
legend(legendStrings);
xlim([0 length(synchContactPixelsArray)]);
title('Please select the start of tangential force loading');
[x,y] = ginput(1);
slipStartInd = round(x);
for i = 1:nPlot
    subplot(nPlot,1,i); hold on;
    plot([slipStartInd slipStartInd],[0,temp(i,slipStartInd)],'g-');
end
legendStrings = [legendStrings, {'Start F_T Loading'}];
subplot(nPlot,1,1); hold on;
legend(legendStrings);
xlim([0 length(synchContactPixelsArray)]);
title('Please select the point of slip');
[x,y] = ginput(1);
slipEndInd = round(x);
for i = 1:nPlot
    subplot(nPlot,1,i); hold on;
    plot([slipEndInd slipEndInd],[0,temp(i,slipEndInd)],'b-');
end
subplot(nPlot,1,1); hold on;
legendStrings = [legendStrings, {'Slip'}];
xlim([0 length(synchContactPixelsArray)]);
legend(legendStrings);