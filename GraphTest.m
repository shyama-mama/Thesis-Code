clc;
clear;


load('C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Force Files\ForcesPlateModAndSynchTrigger.mat');
load('C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Tracked Points\4_12FPS_PlateOn_TrackedPoints.mat');
saveFilename = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Triangles\AreaStress\Graphs\';
fs = samplerate(25);
FxAll = data(datastart(25):dataend(25)); 
FyAll = data(datastart(26):dataend(26)); 
FzAll = data(datastart(27):dataend(27));
TxAll = data(datastart(28):dataend(28));
TyAll = data(datastart(29):dataend(29));
TzAll = data(datastart(30):dataend(30));
plateModAll = data(datastart(31):dataend(31));
forceTrigSig = data(datastart(32):dataend(32));
fps = 120;

forceTrigThresh = 5;
forceTrigStartEndSig = diff(forceTrigSig>forceTrigThresh);
forceTrigStartInds = find(forceTrigStartEndSig==1);
forceTrigEndInds = find(forceTrigStartEndSig==-1);
forceTrigDurations = forceTrigEndInds - forceTrigStartInds;
forceTriggerMidInds = forceTrigStartInds+round(forceTrigDurations./2);

%figure; set(gcf,'Color',[1 1 1]);
%        ax(1) = subplot(2,1,1); plot(forceTrigSig); xlabel('Sample #'); ylabel('Trigger Signal');
%        hold on; plot([0 length(forceTrigSig)],[forceTrigThresh forceTrigThresh],'r');
%        ax(2) = subplot(2,1,2); plot(forceTrigStartEndSig); xlabel('Sample #'); ylabel('Logic Trigger Signal');
%        linkaxes(ax,'x');

[b,a] = butter(3,50./(fs/2),'low');
FxAllFiltered = filtfilt(b,a,FxAll);
FyAllFiltered = filtfilt(b,a,FyAll);
FzAllFiltered = filtfilt(b,a,FzAll);
TxAllFiltered = filtfilt(b,a,TxAll);
TyAllFiltered = filtfilt(b,a,TyAll);
TzAllFiltered = filtfilt(b,a,TzAll);
fprintf('Done loading forces\n\n');


thisStimForceStartSamp = forceTriggerMidInds(25);
thisStimForceEndSamp = forceTriggerMidInds(26);
thisStimForceDurationSamps = (thisStimForceEndSamp - thisStimForceStartSamp);

thisStimFx = FxAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
thisStimFy = FyAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
thisStimFz = FzAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
thisStimTx = TxAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
thisStimTy = TyAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
thisStimTz = TzAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
thisStimPlateMod = plateModAll(thisStimForceStartSamp:thisStimForceEndSamp);

forceSampInds = (1:round(fs./fps):thisStimForceDurationSamps);
synchFps = fps;
b = ones(1,round(fs./fps))./round((fs./fps)); a = [1 zeros(1,round(fs./fps)-1)]; % moving average filter
temp = filter(b,a,thisStimFx);  
synchFx = temp(forceSampInds); % thisStimFx(forceSampInds); 
temp = filter(b,a,thisStimFy);  synchFy = temp(forceSampInds); % thisStimFy(forceSampInds); 
temp = filter(b,a,thisStimFz);  synchFz = temp(forceSampInds); % thisStimFz(forceSampInds); 
temp = filter(b,a,thisStimTx);  synchTx = temp(forceSampInds); % thisStimTx(forceSampInds); 
temp = filter(b,a,thisStimTy);  synchTy = temp(forceSampInds); % thisStimTy(forceSampInds); 
temp = filter(b,a,thisStimTz);  synchTz = temp(forceSampInds); % thisStimTz(forceSampInds); 
temp = filter(b,a,thisStimPlateMod);  synchPlateMod = temp(forceSampInds); % thisStimPlateMod(forceSampInds);
synchFs = synchFps;


FN = -synchFz';
FT = sqrt(synchFx'.^2 + synchFy'.^2);
FTonFN = FT./FN;
synchContactPixelsArray = zeros(size(synchFz))';

figure; set(gcf,'Color',[1 1 1]);
nPlot = 3;
ax(1) = subplot(nPlot,1,1);
plot(FT,'k');
xlabel('Frame Number'); ylabel('F_T (N)');
ax(2) = subplot(nPlot,1,2);
plot(FN,'k');
xlabel('Frame Number'); ylabel('F_N (N)');
ax(3) = subplot(nPlot,1,3);
plot(FTonFN,'k');
xlabel('Frame Number'); ylabel('F_T/F_N');
linkaxes(ax,'x');

save = sprintf('%s11_PlateOff', saveFilename);
print(save,'-dpng');