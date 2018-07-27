% SplitAndSynchForceAndVideo

DataFilename = 'Video 1';%'SampleForSynchTestWithFingerAndPlate2'; %sprintf('%s_%s_%s_offonX4',subject,finger,direction);
RootPath = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Video'; %'C:\HEBA_RESEARCH\TACTILE\HUMAN_FRICTION\CAMERA_FRICTION_PLATE\Victor_RawData';
WritePath = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Split data'; %'C:\HEBA_RESEARCH\TACTILE\HUMAN_FRICTION\CAMERA_FRICTION_PLATE\Victor_SplitData';

readVideoFilename = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Video\Video 1.MP4';
readForceFilename = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Force Files\ForcesPlateModAndSynchTrigger.mat';

% GET STIM START AND END FROM VIDEO

[videoTrigSig,fps] = readVideoForTrigger(readVideoFilename);

[b,a] = butter(3,0.1/(fps/2),'high');
videoTrigSigFilt = filtfilt(b,a,videoTrigSig);
figure; plot(videoTrigSigFilt)
title(sprintf('mean: %.2f, std: %.2f',mean(videoTrigSigFilt),std(videoTrigSigFilt)));
% videoTrigSigFilt = medfilt1(videoTrigSig,2*fps);
% videoTrigSigFilt = videoTrigSig-videoTrigSigFilt;
% videoTrigSigFilt(1) = videoTrigSigFilt(2);
videoTrigThresh = mean(videoTrigSigFilt)+2.5*std(videoTrigSigFilt);
videoTrigStartEndSig = diff(videoTrigSigFilt>videoTrigThresh);
videoTrigStartInds = find(videoTrigStartEndSig==1);
videoTrigEndInds = find(videoTrigStartEndSig==-1);
videoTrigDurations = videoTrigEndInds - videoTrigStartInds;
videoTriggerMidInds = videoTrigStartInds+round(videoTrigDurations./2);

figure; set(gcf,'Color',[1 1 1]);
        ax(1) = subplot(2,1,1); plot(videoTrigSigFilt); xlabel('Frame #'); ylabel('mean Frame intensity');
        hold on; plot([0 length(videoTrigSigFilt)],[videoTrigThresh videoTrigThresh],'r');
        ax(2) = subplot(2,1,2); plot(videoTrigStartEndSig); xlabel('Frame #'); ylabel('Logic Trigger Signal');
        linkaxes(ax,'x');

% GET STIM START AND END FROM FORCES

[FxAll, FyAll, FzAll, TxAll, TyAll, TzAll, plateModAll, forceTrigSig, fs] = readForces(readForceFilename,true);
forceTrigThresh = 5;
forceTrigStartEndSig = diff(forceTrigSig>forceTrigThresh);
forceTrigStartInds = find(forceTrigStartEndSig==1);
forceTrigEndInds = find(forceTrigStartEndSig==-1);
forceTrigDurations = forceTrigEndInds - forceTrigStartInds;
forceTriggerMidInds = forceTrigStartInds+round(forceTrigDurations./2);

figure; set(gcf,'Color',[1 1 1]);
        ax(1) = subplot(2,1,1); plot(forceTrigSig); xlabel('Sample #'); ylabel('Trigger Signal');
        hold on; plot([0 length(forceTrigSig)],[forceTrigThresh forceTrigThresh],'r');
        ax(2) = subplot(2,1,2); plot(forceTrigStartEndSig); xlabel('Sample #'); ylabel('Logic Trigger Signal');
        linkaxes(ax,'x');
        
[b,a] = butter(3,50./(fs/2),'low');
FxAllFiltered = filtfilt(b,a,FxAll);
FyAllFiltered = filtfilt(b,a,FyAll);
FzAllFiltered = filtfilt(b,a,FzAll);
TxAllFiltered = filtfilt(b,a,TxAll);
TyAllFiltered = filtfilt(b,a,TyAll);
TzAllFiltered = filtfilt(b,a,TzAll);
        
videoReader = VideoReader(readVideoFilename);        
height = videoReader.Height;
width = videoReader.Width;

stimNumber = 1;
videoFrameInd = 1;
for i = 1:2:length(videoTriggerMidInds)
    thisStimVideoStartFrame = videoTriggerMidInds(i);
    thisStimVideoEndFrame = videoTriggerMidInds(i+1);
    thisStimVideoDurationFrames = (thisStimVideoEndFrame - thisStimVideoStartFrame+1);
    thisStimVideoDurationSecs = thisStimVideoDurationFrames./fps;
    
    thisStimForceStartSamp = forceTriggerMidInds(i);
    thisStimForceEndSamp = forceTriggerMidInds(i+1);
    thisStimForceDurationSamps = (thisStimForceEndSamp - thisStimForceStartSamp);
    thisStimForceDurationSecs = thisStimForceDurationSamps./fs;
    
    if abs(thisStimVideoDurationSecs - thisStimForceDurationSecs) > 2./fps
        fprintf('Stimulus duration mismatch!');
        keyboard
    end
    
    % Get this stimulus forces
    thisStimFx = FxAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
    thisStimFy = FyAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
    thisStimFz = FzAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
    thisStimTx = TxAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
    thisStimTy = TyAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
    thisStimTz = TzAllFiltered(thisStimForceStartSamp:thisStimForceEndSamp);
    thisStimPlateMod = plateModAll(thisStimForceStartSamp:thisStimForceEndSamp);
    
    % Skip to first frame of stimulus
    while hasFrame(videoReader) && videoFrameInd < thisStimVideoStartFrame
        thisFrame = readFrame(videoReader);
        videoFrameInd = videoFrameInd + 1;
    end
    
    % Read until last frame of stimulus
    stimFrameInd = 1;
    thisStimVideoArray = zeros(thisStimVideoDurationFrames,height,width,3,'uint8');
    while hasFrame(videoReader) && videoFrameInd <= thisStimVideoEndFrame
        thisFrame = readFrame(videoReader);
        thisStimVideoArray(stimFrameInd,:,:,:) = thisFrame;
        videoFrameInd = videoFrameInd + 1;
        stimFrameInd = stimFrameInd + 1;
    end
    
    % Write this stimulus video to avi - not necessary
%     videoWriter = VideoWriter(sprintf('%s\\%s\\%s_%d.avi',WritePath,DataFilename,DataFilename,stimNumber));
%     open(videoWriter);
%     for frameInd = 1:thisStimVideoDurationFrames
%         writeVideo(videoWriter,squeeze(thisStimVideoArray(frameInd,:,:,:)));
%     end
%     close(videoWriter);
    % Save this stimulus video arrays at original sampling rate - not necessary
%     save(sprintf('%s\\%s\\%s_%dvideo.mat',WritePath,DataFilename,DataFilename,stimNumber),'thisStimVideoArray','fps','thisStimVideoStartFrame','thisStimVideoEndFrame');
    % Save this stimulus force arrays at original sampling rate - not necessary
%     save(sprintf('%s\\%s\\%s_%dforce.mat',WritePath,DataFilename,DataFilename,stimNumber),'thisStimFx','thisStimFy','thisStimFz','thisStimTx','thisStimTy','thisStimTz','thisStimPlateMod','fs','thisStimForceStartSamp','thisStimForceEndSamp');
    
    % Synch video and force
    synchVideoArray = thisStimVideoArray;
    % Resample forces with moving average filter
    forceSampInds = (1:(fs./fps):thisStimForceDurationSamps);
    synchFps = fps;
    b = ones(1,fs./fps)./(fs./fps); a = [1 zeros(1,fs./fps-1)]; % moving average filter
    temp = filter(b,a,thisStimFx);  synchFx = temp(forceSampInds); % thisStimFx(forceSampInds); 
    temp = filter(b,a,thisStimFy);  synchFy = temp(forceSampInds); % thisStimFy(forceSampInds); 
    temp = filter(b,a,thisStimFz);  synchFz = temp(forceSampInds); % thisStimFz(forceSampInds); 
    temp = filter(b,a,thisStimTx);  synchTx = temp(forceSampInds); % thisStimTx(forceSampInds); 
    temp = filter(b,a,thisStimTy);  synchTy = temp(forceSampInds); % thisStimTy(forceSampInds); 
    temp = filter(b,a,thisStimTz);  synchTz = temp(forceSampInds); % thisStimTz(forceSampInds); 
    temp = filter(b,a,thisStimPlateMod);  synchPlateMod = temp(forceSampInds); % thisStimPlateMod(forceSampInds);
    synchFs = synchFps;
    
    % Save synched video and force
    if ~exist(sprintf('%s\\%s',WritePath,DataFilename),'dir')
        [SUCCESS,MESSAGE,MESSAGEID] = mkdir(WritePath,DataFilename);
    end       
    save(sprintf('%s\\%s\\%s_%dsynched.mat',WritePath,DataFilename,DataFilename,stimNumber),'synchVideoArray',...
        'synchFps','synchFx','synchFy','synchFz','synchTx','synchTy','synchTz','synchPlateMod','synchFs',...
        'thisStimVideoStartFrame','thisStimVideoEndFrame','thisStimForceStartSamp','thisStimForceEndSamp', '-v7.3');
    
    % FOR MY SANITY
%     newLength = 800; newHeight = 400; % Dimensions of the cropped image
%     xShift = 0; yShift = 0; % Shifts the cropped area
%     [croppedVideoArray] = cropVideo(synchVideoArray,newLength,newHeight,xShift,yShift);
%     se = strel('disk',8);% Structure element for top-hat filter
%     se2 = strel('disk', 5);% Structure element for the closing area
%     A =[20 20];% Matrix to make the border flat
%     otsuThreshold = 0.025;% Otsu level threshold
%     pixelsThreshold = 1000; %10000;% Threshold in number of pixels to delete minor areas
%     [grayImageArray,filteredImageArray,bwImageArray,flatBorderImageArray,contactLogicalImageArray,contactPixelsArray,otsuLevelArray] = filterForContactArea(croppedVideoArray,A,otsuThreshold,pixelsThreshold,se,se2);
%     figure; set(gcf,'Color',[1 1 1]);
%         ax(1) = subplot(4,1,1); plot(contactPixelsArray); xlabel('Frame #'); ylabel('# Pixels');
%         ax(2) = subplot(4,1,2); plot(sqrt(synchFx.^2 + synchFy.^2)); xlabel('Frame #'); ylabel('F_T (N)');
%         ax(3) = subplot(4,1,3); plot(-synchFz); xlabel('Frame #'); ylabel('F_N (N)');
%         ax(4) = subplot(4,1,4); plot(-sqrt(synchFx.^2 + synchFy.^2)./synchFz); xlabel('Frame #'); ylabel('F_T/F_N');
%         linkaxes(ax,'x');
%     keyboard
    
    stimNumber = stimNumber + 1;
end
