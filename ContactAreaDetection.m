clc;
clear;
fprintf('Moar Testing\n');

% Contact Area Params
se = strel('disk',8);% Structure element for top-hat filter
se2 = strel('disk', 5);% Structure element for the closing area
A =[20 20];% Matrix to make the border flat
otsuThreshold = 0.015;%0.025;% Otsu level threshold
pixelsThreshold = 1000; %10000;% Threshold in number of pixels to delete minor areas

synchedFile = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Split data\TestSynch.mat';
cropAndMarkerFilename = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Split data\cropAndMarker.mat';


fprintf('Loading Synched File...');
load(synchedFile);
fprintf('Done\n');

figure; imshow(squeeze(synchVideoArray(midInd,:,:,:)));
title('Please, select the top left corner of cropped image');
[x1,y1] = ginput(1); xmin_crop = uint16(round(x1)); ymin_crop = uint16(round(y1));
title('Please, select the bottom right corner of cropped image');
[x2,y2] = ginput(1); xmax_crop = uint16(round(x2)); ymax_crop = uint16(round(y2));

[initContactStartInd,initContactEndInd,slipStartInd,slipEndInd] = getTrackingStartAndEndInds(synchFx, synchFy, synchFz, zeros(size(synchFz))');


%{
 Gets the start and end indexes of initial contact 
synchContactPixelsArray = zeros(size(synchFz)); 
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

subplot(nPlot,1,1); hold on;
legendStrings = [legendStrings, {'Slip'}];
xlim([0 length(synchContactPixelsArray)]);
legend(legendStrings);
%}

fprintf('Start of intitial contact is %d and end of initial contact is %d\n', initContactStartInd, initContactEndInd);


