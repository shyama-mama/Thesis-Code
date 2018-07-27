clc;
clear all;
close all;

%%  ==============================================================
%======================= CONSTANTS ===============================
%=================================================================

savePath = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis';
dataPath = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Split data\'; 

%%  ==============================================================
%======================= INPUTS ==================================
%=================================================================

subjectNum = 1;
stimOrderNum = 1;
stimulusName = 'Stuff';
dateNum = 131115;
fileName = sprintf('Subject%d_StimOrder%d_%d',subjectNum,stimOrderNum,dateNum);
cropAndMarkerFilename = sprintf('%scropAndMarker.mat',dataPath); 

if stimOrderNum == 1 || stimOrderNum == 2
    nStim = 1;
else
    nStim = 48;
end
for stimulusNumber = 1:nStim
%%  ==============================================================
%===================== PARAMETERS ================================
%=================================================================

% Input Filenames
synchedFilename = sprintf('%sMovingPingas.mat',dataPath);
startEndContactSlipFilename = sprintf('%scontactInds.mat',dataPath); 

% Output Filenames
backwardsMatlabRoute = sprintf('%s\\%s_BackwardsTrackingData.mat',savePath,stimulusName);% Route of the backwards tracking data
forwardsMatlabRoute = sprintf('%s\\%s_ForwardsTrackingData.mat',savePath,stimulusName);% Route of the forwards tracking data
plateauMatlabRoute = sprintf('%s\\%s_PlateauTrackingData.mat',savePath,stimulusName);% Route of the backwards tracking data
initContOnOffCompMatlabRoute = sprintf('%s\\%s_InitContOnOffCompTrackingData.mat',savePath,stimulusName);% Route of the backwards tracking data

% Contact Area Params
se = strel('disk',8);% Structure element for top-hat filter
se2 = strel('disk', 5);% Structure element for the closing area
A =[20 20];% Matrix to make the border flat
otsuThreshold = 0.015;%0.025;% Otsu level threshold
pixelsThreshold = 1000; %10000;% Threshold in number of pixels to delete minor areas

% Point Tracking Params
artefactsThreshold = 0.5;% Maximum bidirectional error in pixels
distanceThreshold = 1;% Distance threshold between marker and tracked point to determine slip

%%  ==============================================================
%============ READ RAW VIDEO, UNDISTORT and CROP =================
%=================================================================

% READ RAW VIDEO
if exist(synchedFilename,'file')
    fprintf('Loading synched data...');
    load(synchedFilename);
    fprintf('Done.\n');
else
    fprintf('ERROR: Synched .mat file does not exist. Did you run SplitAndSynchForceAndVideo.mat?\n');
    return;
end

% ESTIMATE CAMERA PARAMETERS AND PIXELS PER MM
checkerboardFiles = {sprintf('%sang1bottomright.jpg',dataPath);...
    sprintf('%sang2topleft.jpg',dataPath);...
    sprintf('%sang4bottomleft.jpg',dataPath);...
    sprintf('%sbottomleft.jpg',dataPath);...
    sprintf('%sbottomright.jpg',dataPath);...
    sprintf('%sleft2.jpg',dataPath);...
    sprintf('%stop.jpg',dataPath);...
    sprintf('%stop2.jpg',dataPath);...
    sprintf('%stopleft.jpg',dataPath);...
    sprintf('%stopright.jpg',dataPath);...
    };
%[cameraParams, pixelsSquareMm] = estimateCameraParams(checkerboardFiles,checkerboardOnPlate,squareSize);

% UNDISTORT VIDEO
% [undistortedVideoArray] = undistortVideo(rawVideoArray,cameraParams);
undistortedVideoArray = synchVideoArray;

% CROP VIDEO
% Cropping Params
if exist(cropAndMarkerFilename,'file') %
    load(cropAndMarkerFilename);
elseif stimulusNumber == 1
    midInd = uint16(round(size(undistortedVideoArray,1)./2));
    figure; imshow(squeeze(undistortedVideoArray(midInd,:,:,:)));
    title('Please, select the top left corner of cropped image');
    [x1,y1] = ginput(1); xmin_crop = uint16(round(x1)); ymin_crop = uint16(round(y1));
    title('Please, select the bottom right corner of cropped image');
    [x2,y2] = ginput(1); xmax_crop = uint16(round(x2)); ymax_crop = uint16(round(y2));
end
% newLength = 500; newHeight = 400;%400; % Dimensions of the cropped image
% xShift = 100; yShift = -185; % Shifts the cropped area
[croppedVideoArray] = cropVideo(undistortedVideoArray,ymin_crop,ymax_crop,xmin_crop,xmax_crop);%newLength,newHeight,xShift,yShift);

if ~exist(cropAndMarkerFilename,'file') && stimulusNumber == 1
    figure; imshow(squeeze(croppedVideoArray(midInd,:,:,:)));
    title('Please, select the marker');
    [x_marker,y_marker] = ginput(1);
    save(cropAndMarkerFilename,'xmin_crop','xmax_crop','ymin_crop','ymax_crop','x_marker','y_marker');
end


if stimOrderNum == 1 || stimOrderNum == 2
    %%  ==============================================================
    %======== GET START & END OF INIT CONTACT AND SLIP ===============
    %=================================================================

    if exist(startEndContactSlipFilename,'file')
        load(startEndContactSlipFilename);
    else
        [initContactStartInd,initContactEndInd,slipStartInd,slipEndInd] = getTrackingStartAndEndInds(synchFx, synchFy, synchFz, zeros(size(synchFz))');
        save(startEndContactSlipFilename,'initContactStartInd','initContactEndInd','slipStartInd','slipEndInd');
    end

    %%  ==============================================================
    %================= POINT TRACKING BACKWARDS ======================
    %=================================================================

    offsetFx = mean(synchFx((initContactStartInd-4):initContactStartInd));
    offsetFy = mean(synchFy((initContactStartInd-4):initContactStartInd));
    offsetFz = mean(synchFz((initContactStartInd-4):initContactStartInd));

    % Track from slipStartInd forwards to slipEndInd
    forwardInds = slipStartInd:slipEndInd;
    Fx = synchFx(forwardInds)-offsetFx; Fy = synchFy(forwardInds)-offsetFy; Fz = synchFz(forwardInds)-offsetFz;
    [grayTrackingArray,filteredTrackingArray,bwTrackingArray,flatBorderTrackingArray,contactAreaTrackingArray,otsuLevelStructArray,...
        pointLocsPerFrame,pointValidityPerFrame,xDistToMarker,yDistToMarker,movementDistance,hasChanged,hasSlipped,hasBeenOutOfContactArea,trackingImageArray]...
        =trackPoints(artefactsThreshold,distanceThreshold,croppedVideoArray(forwardInds,:,:,:),x_marker,y_marker,A,otsuThreshold,pixelsThreshold,se,se2,true,false);

    fprintf('Saving slip tracking...');
    save(forwardsMatlabRoute,'Fx','Fy','Fz','grayTrackingArray','filteredTrackingArray','bwTrackingArray','flatBorderTrackingArray','contactAreaTrackingArray','trackingImageArray','pointLocsPerFrame','pointValidityPerFrame','otsuLevelStructArray', '-v7.3');
    fprintf('Done.\n');

    % Track from initContactEndInd backwards to initContactStartInd-1
    otsuLevelStruct = otsuLevelStructArray(1);
    backwardInds = initContactEndInd:-1:(initContactStartInd);
    Fx = synchFx(backwardInds)-offsetFx; Fy = synchFy(backwardInds)-offsetFy; Fz = synchFz(backwardInds)-offsetFz;
    [grayTrackingArray,filteredTrackingArray,bwTrackingArray,flatBorderTrackingArray,contactAreaTrackingArray,otsuLevelStructArray,...
        pointLocsPerFrame,pointValidityPerFrame,xDistToMarker,yDistToMarker,movementDistance,hasChanged,hasSlipped,hasBeenOutOfContactArea,trackingImageArray]...
        =trackPoints(artefactsThreshold,distanceThreshold,croppedVideoArray(backwardInds,:,:,:),x_marker,y_marker,A,otsuThreshold,pixelsThreshold,se,se2,true,false,otsuLevelStruct);

    fprintf('Saving initial contact tracking...');
    save(backwardsMatlabRoute,'Fx','Fy','Fz','grayTrackingArray','filteredTrackingArray','bwTrackingArray','flatBorderTrackingArray','contactAreaTrackingArray','trackingImageArray','pointLocsPerFrame','pointValidityPerFrame','otsuLevelStructArray');
    fprintf('Done.\n');
    

    % Image alignment and comparison init contact: final contact compared
    % to very first contact area
%     [grayFrame0,filteredFrame0,bwFrame0,flatBorderFrame0,contactAreaFrame0,contactPixels0,otsuLevelStruct] = filterForContactArea(croppedVideoArray(initContactEndInd,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,0);
%     grayFrame0 = squeeze(grayFrame0);
%     filteredFrame0 = squeeze(filteredFrame0);
%     bwFrame0 = squeeze(bwFrame0);
%     flatBorderFrame0 = squeeze(flatBorderFrame0);
%     contactAreaFrame0 = squeeze(contactAreaFrame0);
%     contactAreaFrame0 = ignoreMultipleContact(contactAreaFrame0,grayFrame0,flatBorderFrame0);
%     
%     fixed = (bwFrame0).*((contactAreaFrame0));
%     fixed1 = (grayFrame0).*uint8((contactAreaFrame0));
%     
%     f1 = figure; 
%     for frameInd=(initContactEndInd-1):-1:(initContactStartInd-1)
%         [grayFrame1,filteredFrame1,bwFrame1,flatBorderFrame1,contactAreaFrame1,contactPixels1,tempOtsuLevelStruct1] = filterForContactArea(croppedVideoArray(frameInd,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,otsuLevelStruct);
%         grayFrame1 = squeeze(grayFrame1);
%         contactAreaFrame1 = squeeze(contactAreaFrame1);
%         contactAreaFrame1 = ignoreMultipleContact(contactAreaFrame1,grayFrame1,flatBorderFrame1);
%         
%         imshow(contactAreaFrame1); pause(0.1);
%         if sum(sum(contactAreaFrame1 & contactAreaFrame0)) == 0
%             break;
%         end
%     end    
%     close(f1);
%     frameInd = frameInd + 1;
%         
%     [grayFrame1,filteredFrame1,bwFrame1,flatBorderFrame1,contactAreaFrame1,contactPixels1,tempOtsuLevelStruct1] = filterForContactArea(croppedVideoArray(frameInd,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,otsuLevelStruct);
%     grayFrame1 = squeeze(grayFrame1);
%     filteredFrame1 = squeeze(filteredFrame1);
%     bwFrame1 = squeeze(bwFrame1);
%     flatBorderFrame1 = squeeze(flatBorderFrame1);
%     contactAreaFrame1 = squeeze(contactAreaFrame1);
%     contactAreaFrame1 = ignoreMultipleContact(contactAreaFrame1,grayFrame1,flateBorderFrame1);
%     moving = (bwFrame1).*((contactAreaFrame1));
%     moving1 = (grayFrame1).*uint8((contactAreaFrame1));
% 
% 
%     tempFixed1 = uint8(~(bwFrame0.*contactAreaFrame0));%bwmorph(~(bwFrame0.*contactAreaFrame0),'thin',inf).*contactAreaFrame0; %thins the fingerprint rigde to a single pixel
%     tempMoving1 = uint8(~(bwFrame1.*contactAreaFrame1));%bwmorph(~(bwFrame1.*contactAreaFrame1),'thin',inf).*contactAreaFrame1; %thins the fingerprint rigde to a single pixel
%     [optimizer,metric] = imregconfig('monomodal');%'monomodal;);
%     tempTform = imregtform(tempMoving1, tempFixed1, 'rigid', optimizer, metric); % rigid = rotation and translation, no scaling
%     tempMovingRegistered = imwarp(tempMoving1,tempTform,'OutputView',imref2d(size(tempFixed1)));
%     figure; subplot(1,2,1);
%     h = imshowpair(tempMoving1, tempFixed1, 'falsecolor');
%     colImage = get(h,'CData');
%     imshow(colImage);%y);
%     subplot(1,2,2);
%     h = imshowpair(tempMovingRegistered, tempFixed1, 'falsecolor');
%     colImage = get(h,'CData');
%     imshow(colImage);%y);
%     
%     [optimizer,metric] = imregconfig('monomodal');%'monomodal;);
% %     optimizer.InitialRadius = optimizer.InitialRadius/3.5;
% %     optimizer.MaximumIterations = 300;
%     tform = imregtform(moving1, fixed1, 'rigid', optimizer, metric); % rigid = rotation and translation, no scaling
%     movingRegisteredAdjustedInitRad300 = imwarp(moving1,tform,'OutputView',imref2d(size(fixed1)));
%     figure; subplot(1,2,1);
%     h = imshowpair(moving1, fixed1, 'falsecolor');
%     colImage = get(h,'CData');
%     hgamma = vision.GammaCorrector(10.0,'Correction','De-gamma');%'Gamma');%
%     y = step(hgamma, colImage);
%     imshow(y);
%     subplot(1,2,2);
%     h = imshowpair(movingRegisteredAdjustedInitRad300, fixed1, 'falsecolor');
%     colImage = get(h,'CData');
%     hgamma = vision.GammaCorrector(10.0,'Correction','De-gamma');%'Gamma');%
%     y = step(hgamma, colImage);
%     imshow(y);
% 
%     keyboard
    
     % Image alignment and comparison init contact: final contact compared
    % to previous frames
%     if stimulusNumber == 1
%         [grayFrame0,filteredFrame0,bwFrame0,flatBorderFrame0,contactAreaFrame0,contactPixels0,otsuLevelStruct] = filterForContactArea(croppedVideoArray(initContactEndInd,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,0);
%     else
%         [grayFrame0,filteredFrame0,bwFrame0,flatBorderFrame0,contactAreaFrame0,contactPixels0,otsuLevelStruct] = filterForContactArea(croppedVideoArray(initContactEndInd,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,otsuLevelStruct);
%     end
%     grayFrame0 = squeeze(grayFrame0);
%     filteredFrame0 = squeeze(filteredFrame0);
%     bwFrame0 = squeeze(bwFrame0);
%     flatBorderFrame0 = squeeze(flatBorderFrame0);
%     contactAreaFrame0 = squeeze(contactAreaFrame0);
%     contactAreaFrame0 = ignoreMultipleContact(contactAreaFrame0,grayFrame0,flatBorderFrame0);
%     
%     s = regionprops(contactAreaFrame0,'Centroid');
%     c0 = s.Centroid;
%     
%     fixed = (bwFrame0).*((contactAreaFrame0));
%     fixed1 = (grayFrame0).*uint8((contactAreaFrame0));
%     
%     
%     backwardInds = initContactEndInd:-1:(initContactStartInd);
%     Fx = synchFx-offsetFx; Fy = synchFy-offsetFy; Fz = synchFz-offsetFz;
%     
%     centroidDist{stimulusNumber,1} = [];
%     centroidPos{stimulusNumber,1} = [];
%     centroidPos{stimulusNumber,1} = [centroidPos{stimulusNumber,1}; c0(:)'];
%     contactSize{stimulusNumber,1} = [];
%     contactSize{stimulusNumber,1} = [contactSize{stimulusNumber,1},contactPixels0];
%     ft{stimulusNumber,1} = [];
%     fn{stimulusNumber,1} = [];
% %     f2 = figure;
%     for frameInd=(initContactEndInd-1):-1:(initContactStartInd-1)
%         
%         [grayFrame1,filteredFrame1,bwFrame1,flatBorderFrame1,contactAreaFrame1,contactPixels1,tempOtsuLevelStruct1] = filterForContactArea(croppedVideoArray(frameInd,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,otsuLevelStruct);
%         grayFrame1 = squeeze(grayFrame1);
%         filteredFrame1 = squeeze(filteredFrame1);
%         bwFrame1 = squeeze(bwFrame1);
%         flatBorderFrame1 = squeeze(flatBorderFrame1);
%         contactAreaFrame1 = squeeze(contactAreaFrame1);
%         contactAreaFrame1 = ignoreMultipleContact(contactAreaFrame1,grayFrame1,flatBorderFrame1);
%         
%         if sum(sum(contactAreaFrame1 & contactAreaFrame0)) == 0
%             break;
%         end
%         
%         s = regionprops(contactAreaFrame1,'Centroid');
%         c1 = s.Centroid; 
%         
%         d = sqrt((c0(1)-c1(1)).^2+(c0(2)-c1(2)).^2);
%         centroidPos{stimulusNumber,1} = [centroidPos{stimulusNumber,1}; c1(:)'];
%         centroidDist{stimulusNumber,1} = [centroidDist{stimulusNumber,1}, d];
%         contactSize{stimulusNumber,1} = [contactSize{stimulusNumber,1},contactPixels1];
%         ft{stimulusNumber,1} = [ft{stimulusNumber,1}, sqrt(Fx(frameInd).^2+Fy(frameInd).^2)];
%         fn{stimulusNumber,1} = [fn{stimulusNumber,1}, -Fz(frameInd)];
%         
%         moving = (bwFrame1).*((contactAreaFrame1));
%         moving1 = (grayFrame1).*uint8((contactAreaFrame1));
        
        % Plot frame by frame centroids
%         figure(f2);
%         imshowpair(moving1,fixed1,'falsecolor');
%         hold on;
%         line([c0(1), c1(1)], [c0(2), c1(2)], 'LineWidth', 2, 'Color', [0 0 1]); 
%         pause(0.1);

%         fixedFFT = fft(fixed1);
%         movingFFT = fft(moving1);
%         
%         fixedFFT2 = fft2(fixed1);
%         movingFFT2 = fft2(moving1);
%         
%         figure; ax(1) = subplot(1,2,1); imshow(abs(fixedFFT2)./max(max(abs(fixedFFT2))));
%             ax(2) = subplot(1,2,2); imshow(abs(movingFFT2)./max(max(abs(movingFFT2))));
           
%         [optimizer,metric] = imregconfig('multimodal');%'monomodal;);
%         optimizer.InitialRadius = optimizer.InitialRadius/3.5;
%         optimizer.MaximumIterations = 300;
%         tform = imregtform(moving1, fixed1, 'rigid', optimizer, metric); % rigid = rotation and translation, no scaling
%         movingRegisteredAdjustedInitRad300 = imwarp(moving1,tform,'OutputView',imref2d(size(fixed1)));
% 
%         figure;
%         h = imshowpair(movingRegisteredAdjustedInitRad300, fixed1, 'falsecolor');
%         colImage = get(h,'CData');
%         hgamma = vision.GammaCorrector(10.0,'Correction','De-gamma');%'Gamma');%
%         y = step(hgamma, colImage);
%         imshow(y);
%     
%         keyboard
%     end
    
    % Plot centroid distance and Ft and Fn
    
%     figure; subplot(3,1,1); plot(centroidDist{stimulusNumber,1}); xlabel('sample #'); ylabel('Centroid distance (Pi)');
%     hold on; subplot(3,1,2); plot(ft{stimulusNumber,1}); xlabel('sample #'); ylabel('F_T (N)');
%     hold on; subplot(3,1,3); plot(fn{stimulusNumber,1}); xlabel('sample #'); ylabel('F_N (N)');
%     if ~isempty(centroidDist{stimulusNumber,1})
%         centroidFinalSubInitDist(stimulusNumber,1) = centroidDist{stimulusNumber,1}(end);
%         centroidTrajectoryDist(stimulusNumber,1) = sum(centroidDist{stimulusNumber,1});
%         ftFinal(stimulusNumber,1) = ft{stimulusNumber,1}(1);
%         fnFinal(stimulusNumber,1) = fn{stimulusNumber,1}(1);  
%     else
%         centroidFinalSubInitDist(stimulusNumber,1) = nan;
%         centroidTrajectoryDist(stimulusNumber,1) = nan;
%         ftFinal(stimulusNumber,1) = nan;
%         fnFinal(stimulusNumber,1) = nan;
%     end
    %%  ==============================================================
    %======= POINT TRACKING INIT CONTACT ON/OFF COMPARISON ===========
    %=================================================================
    
%     isTrack = false; % false = use imshowpair, true = use point tracking
%     
%     croppedAccum(stimulusNumber,:,:,:) = croppedVideoArray(initContactEndInd,:,:,:);
%     x_markerAccum(stimulusNumber,1) = x_marker;
%     y_markerAccum(stimulusNumber,1) = y_marker;
% 
%     save(sprintf('%s_endInitContact_OnVsOff_FrameByFrame.mat',fileName),'ftFinal','fnFinal','croppedAccum','x_markerAccum','y_markerAccum');
%     
%     if stimulusNumber == nStim
%                 
%         figure; 
%         subplot(2,1,1);
%         plot(centroidFinalDist);
%         subplot(2,1,2);
%         plot(ftFinal);
%                 
%         figure; 
%         subplot(2,1,1);
%         plot(centroidFinalDist(1:2:end),'r'); hold on;
%         plot(centroidFinalDist(2:2:end),'b');
%         subplot(2,1,2);
%         plot(ftFinal(1:2:end),'r'); hold on;
%         plot(ftFinal(2:2:end),'b');
%         
%         figure; 
%         for stimulusNumber = 1:96
%             centroidStep{stimulusNumber,1} = diff(centroidDist{stimulusNumber,1});
%             if mod(stimulusNumber,2) == 0
%                 plot(centroidStep{stimulusNumber-1,1},'r');
%                 hold on;
%                 plot(centroidStep{stimulusNumber,1},'b');
%                 hold off;
%                 keyboard
%             end
%         end
%         
%         save(sprintf('%s_%s.mat',fileName,'centroids'),'centroidPos','ft','fn','contactSize');
%         keyboard;
        
%         for i = 1:1:nStim
%             
%             if isTrack
%                 [grayTrackingArray,filteredTrackingArray,bwTrackingArray,flatBorderTrackingArray,contactAreaTrackingArray,otsuLevelStructArray,...
%                 pointLocsPerFrame,pointValidityPerFrame,xDistToMarker,yDistToMarker,movementDistance,hasChanged,hasSlipped,hasBeenOutOfContactArea,trackingImageArray]...
%                 =trackPoints(inf,distanceThreshold,croppedAccum(i:(i+1),:,:,:),x_markerAccum(1,1),y_markerAccum(1,1),A,otsuThreshold,pixelsThreshold,se,se2,true,false);
% 
%                 [xDistToMarker,yDistToMarker,movementDistance,isInContact,hasChanged,hasBeenOutOfContactArea,hasSlipped,videoArray] =  ...
%                 measureSlip(grayTrackingArray,filteredTrackingArray,bwTrackingArray,flatBorderTrackingArray,contactAreaTrackingArray,trackingImageArray,pointLocsPerFrame,pointValidityPerFrame,distanceThreshold);
% 
%                 figure; subplot(1,2,1); imshow(squeeze(videoArray(1,:,:,:))); 
%                 subplot(1,2,2); imshow(squeeze(videoArray(2,:,:,:)));
%             else
%                 [grayFrame0,filteredFrame0,bwFrame0,flatBorderFrame0,contactAreaFrame0,contactPixels0,otsuLevelStruct] = filterForContactArea(croppedAccum(i,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,0);
%                 grayFrame0 = squeeze(grayFrame0);
%                 filteredFrame0 = squeeze(filteredFrame0);
%                 bwFrame0 = squeeze(bwFrame0);
%                 flatBorderFrame0 = squeeze(flatBorderFrame0);
%                 contactAreaFrame0 = squeeze(contactAreaFrame0);
%                 contactAreaFrame0 = ignoreMultipleContact(contactAreaFrame0,grayFrame0);
%                 fixed = (bwFrame0).*((contactAreaFrame0));
%                 fixed1 = (grayFrame0).*uint8((contactAreaFrame0));
%                 
%                 [grayFrame1,filteredFrame1,bwFrame1,flatBorderFrame1,contactAreaFrame1,contactPixels1,tempOtsuLevelStruct1] = filterForContactArea(croppedAccum(i+1,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,otsuLevelStruct);
%                 grayFrame1 = squeeze(grayFrame1);
%                 filteredFrame1 = squeeze(filteredFrame1);
%                 bwFrame1 = squeeze(bwFrame1);
%                 flatBorderFrame1 = squeeze(flatBorderFrame1);
%                 contactAreaFrame1 = squeeze(contactAreaFrame1);
%                 contactAreaFrame1 = ignoreMultipleContact(contactAreaFrame1,grayFrame1);
%                 moving = (bwFrame1).*((contactAreaFrame1));
%                 moving1 = (grayFrame1).*uint8((contactAreaFrame1));
%                 
%                 if mod(i,2) == 0
%                     temp = fixed1;
%                     fixed1 = moving1;
%                     moving1 = temp;
%                 end
%                 fixed1(fixed1 == 0) = nan;
%                 moving1(moving1 == 0) = nan;
%                 [optimizer,metric] = imregconfig('multimodal');%'monomodal;);
%                 optimizer.InitialRadius = optimizer.InitialRadius/1000;%3.5;
%                 optimizer.GrowthFactor = 1.001;
%                 optimizer.MaximumIterations = 1000;
%                 tform = imregtform(moving1, fixed1, 'rigid', optimizer, metric); % rigid = rotation and translation, no scaling
%                 movingRegistered = imwarp(moving1,tform,'OutputView',imref2d(size(fixed1)));
%       
%                 allImages(i,1,:,:) = fixed1;
%                 allImages(i,2,:,:) = moving1;
%                 allImages(i,3,:,:) = movingRegistered;
%                 
%             end
%         
%         end
%         
%         keyboard;
%         
%         save(sprintf('%s_endInitContact_OnVsOff_LastFrame.mat',fileName),'allImages');
%                 
%         fprintf('Saving initial contact On/Off comparison tracking...');
%         save(initContOnOffCompMatlabRoute,'Fx','Fy','Fz','grayTrackingArray','filteredTrackingArray','bwTrackingArray','flatBorderTrackingArray','contactAreaTrackingArray','trackingImageArray','pointLocsPerFrame','pointValidityPerFrame');
%         fprintf('Done.\n');
%     end
    
else
    plateOnSig = synchPlateMod > 5/2;
    
    plateOnStart = find(diff(plateOnSig)==1);
    plateOffStart = find(diff(plateOnSig)==-1);
    if length(plateOnStart) ~= 2 || length(plateOffStart) ~= 2
        fprintf('ERROR: Expected 2 repetitions of plate switching during plateau.\n');
        keyboard
    end
    offWidth = plateOnStart(2) - plateOffStart(1);
    for i = 1:length(plateOnStart)
        onWidth(i) = plateOffStart(i) - plateOnStart(i);
    end
    offFrameInds = [plateOnStart(1) - 5, plateOnStart(2) - uint16(round(offWidth./2)), plateOffStart(2) + 5];
    onFrameInds = [plateOnStart(1) + uint16(round(onWidth(1)./2)), plateOnStart(2) + uint16(round(onWidth(2)./2))];
    
    offsetFx = synchFx(min(onWidth));
    offsetFy = synchFy(min(onWidth));
    offsetFz = synchFz(min(onWidth));
    for i = 1:length(offFrameInds)
        inds = (-2:2) + double(offFrameInds(i));%(-round((onWidth(1)/2)):round((onWidth(1)/2))) + double(offFrameInds(i));
        Fx(i) = median(synchFx(inds)) - offsetFx;
        Fy(i) = median(synchFy(inds)) - offsetFy;
        Fz(i) = median(synchFz(inds)) - offsetFz;
    end
    [grayTrackingArray,filteredTrackingArray,bwTrackingArray,flatBorderTrackingArray,contactAreaTrackingArray,otsuLevelStructArray,...
        pointLocsPerFrame,pointValidityPerFrame,xDistToMarker,yDistToMarker,movementDistance,hasChanged,hasSlipped,hasBeenOutOfContactArea,trackingImageArray]...
        =trackPoints(artefactsThreshold,distanceThreshold,croppedVideoArray(offFrameInds,:,:,:),x_marker,y_marker,A,otsuThreshold,pixelsThreshold,se,se2,false,true);
    
    onFrameInds = [plateOnStart(1) + uint16(round(onWidth(1)./2)), plateOnStart(2) + uint16(round(onWidth(2)./2))];
    for i = 1:length(onFrameInds)
        inds = (-2:2) + double(onFrameInds(i));%(-round((onWidth(1)/2)):round((onWidth(1)/2))) + double(offFrameInds(i));
        onFx(i) = median(synchFx(inds)) - offsetFx;
        onFy(i) = median(synchFy(inds)) - offsetFy;
        onFz(i) = median(synchFz(inds)) - offsetFz;
    end
    otsuLevelStruct = otsuLevelStructArray(1);
    otsuLevelStruct.initContactArea = squeeze(contactAreaTrackingArray(1,:,:));
    [onGrayTrackingArray,onFilteredTrackingArray,onBwTrackingArray,onFlatBorderTrackingArray,onContactAreaTrackingArray,onOtsuLevelStructArray,...
        onPointLocsPerFrame,onPointValidityPerFrame,onXDistToMarker,onYDistToMarker,onMovementDistance,onHasChanged,onHasSlipped,onHasBeenOutOfContactArea,onTrackingImageArray]...
        =trackPoints(artefactsThreshold,distanceThreshold,croppedVideoArray(onFrameInds,:,:,:),x_marker,y_marker,A,otsuThreshold,pixelsThreshold,se,se2,false,true,otsuLevelStruct);
    
    fprintf('Saving plateau tracking...');
    save(plateauMatlabRoute,'Fx','Fy','Fz','grayTrackingArray','filteredTrackingArray','bwTrackingArray','flatBorderTrackingArray','contactAreaTrackingArray','trackingImageArray','pointLocsPerFrame','pointValidityPerFrame',...
        'onFx','onFy','onFz','onGrayTrackingArray','onFilteredTrackingArray','onBwTrackingArray','onFlatBorderTrackingArray','onContactAreaTrackingArray','onTrackingImageArray','onPointLocsPerFrame','onPointValidityPerFrame','otsuLevelStructArray');
    fprintf('Done.\n');
end

end