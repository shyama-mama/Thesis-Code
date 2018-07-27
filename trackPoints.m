function [grayTrackingArray,filteredTrackingArray,bwTrackingArray,flatBorderTrackingArray,contactAreaTrackingArray,otsuLevelStructArray,...
    pointLocsPerFrame,pointValidityPerFrame,xDistToMarker,yDistToMarker,movementDistance,hasChanged,hasSlipped,hasBeenOutOfContactArea,trackedImageArray]...
    = trackPoints(artefactsThreshold,distanceThreshold,croppedVideoArray,markerX,markerY,A,otsuThreshold,pixelsThreshold,se,se2,isIgnoreMultipleContact,isFitEllipse,otsuLevelStruct0)

isPlotIntermediate = false;
isPlot = false;

fprintf('Tracking Points...');

nFrames = size(croppedVideoArray,1);
imageLength = size(croppedVideoArray,2);
imageHeight = size(croppedVideoArray,3);

grayTrackingArray = uint8(zeros(nFrames,imageLength,imageHeight));
filteredTrackingArray = uint8(zeros(nFrames,imageLength,imageHeight));
bwTrackingArray = false(nFrames,imageLength,imageHeight);
flatBorderTrackingArray = false(nFrames,imageLength,imageHeight);
contactAreaTrackingArray = false(nFrames,imageLength,imageHeight);

if exist('otsuLevelStruct0','var')
    [grayFrame0,filteredFrame0,bwFrame0,flatBorderFrame0,contactAreaFrame0,contactPixels0,tempOtsuLevelStruct0] = filterForContactArea(croppedVideoArray(1,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,otsuLevelStruct0);
    otsuLevelStruct = otsuLevelStruct0;
else
    [grayFrame0,filteredFrame0,bwFrame0,flatBorderFrame0,contactAreaFrame0,contactPixels0,otsuLevelStruct0] = filterForContactArea(croppedVideoArray(1,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,0);
end
grayFrame0 = squeeze(grayFrame0);
filteredFrame0 = squeeze(filteredFrame0);
bwFrame0 = squeeze(bwFrame0);
flatBorderFrame0 = squeeze(flatBorderFrame0);
contactAreaFrame0 = squeeze(contactAreaFrame0);

% POINTS SELECTION

if isIgnoreMultipleContact
    % if multiple contact regions, pick the one with the biggest contrast
    s = regionprops(contactAreaFrame0,'all');
    nRegion = length(s);
    if nRegion == 0
        contactAreaFrame0 = flatBorderFrame0;
        s = regionprops(contactAreaFrame0,'all');
        nRegion = length(s);
        if nRegion == 0
            fprintf('Cannot detect a contact area\n');
            keyboard;
        end
    end

    if isPlotIntermediate
        figure; set(gcf,'Color',[1 1 1]);
        subplot((nRegion+1),3,1); imshow(grayFrame0);
        subplot((nRegion+1),3,2); imshow(contactAreaFrame0);
    end
    for i = 1:nRegion
        pixIdx = s(i).PixelIdxList;
        thisRegionContactArea = s(i).FilledImage;
        thisRegionContactArea = false(size(contactAreaFrame0)); thisRegionContactArea(pixIdx) = true;   
        temp = double(grayFrame0); temp(~thisRegionContactArea)=nan;
        intensity(i) = nanmedian(temp(:));
        if isPlotIntermediate
            subplot((nRegion+1),3,(i)*3+2); 
                imshow(thisRegionContactArea); 
                title(sprintf('Region %d Contact Area',i));
            subplot((nRegion+1),3,(i)*3+3);
                hist(temp(:)); 
                title(sprintf('Region %d Intensity Dist.',i));
                legend(sprintf('Median: %d',intensity(i)));
        end
    end
    [maxVal,regionInd] = max(intensity);
    pixIdx = s(regionInd).PixelIdxList;
    contactAreaFrame0 = false(size(contactAreaFrame0)); contactAreaFrame0(pixIdx) = true;
end

if isFitEllipse
    contactAreaFrame0 = myFitEllipseToContact(contactAreaFrame0);
end

if exist('otsuLevelStruct0','var') && isfield(otsuLevelStruct0,'initContactArea')
    contactAreaFrame0 = otsuLevelStruct0.initContactArea;
%     figure; hparent = imshow(grayFrame0);
%     h = imellipse;
%     vert = wait(h);
%     p = h.getPosition();
%     
%     [x,y] = addEllipseToPlot([p(1)+p(3)./2,p(2)+p(4)./2],max(p(3:4)),min(p(3:4)),90,false);
%     ellipseIm = false(size(contactAreaFrame0)); 
%     for i = 1:length(x) 
%         xi = min(max(round(y(i)),1),size(ellipseIm,1));
%         yi = min(max(round(x(i)),1),size(ellipseIm,2));
%         ellipseIm(xi,yi) = true; 
%     end
%     ellipseContactFrame = imfill(ellipseIm,'holes');
%     contactAreaFrame0 = ellipseContactFrame;
end
trackingArea = imdilate(contactAreaFrame0,strel('disk',50));
grayTrackingArray(1,:,:) = grayFrame0;
filteredTrackingArray(1,:,:) = filteredFrame0;
bwTrackingArray(1,:,:) = bwFrame0;
flatBorderTrackingArray(1,:,:) = flatBorderFrame0;
contactAreaTrackingArray(1,:,:) = contactAreaFrame0;

% Points to track
[pointsToTrackX,pointsToTrackY] = selectPointsToTrack(trackingArea,grayFrame0);
% Marker
% [markerX,markerY] = selectMarker(imageHeight/2,imageLength - 25,grayFrame0);
x = [markerX, pointsToTrackX];
y = [markerY, pointsToTrackY];
nPoints = length(x);

% INITIALISING ARRAYS

pointLocsPerFrame = zeros(nFrames,nPoints,2);
pointValidityPerFrame = zeros(nFrames,nPoints);

xDistToMarker = zeros(nFrames,nPoints);
yDistToMarker = zeros(nFrames,nPoints);
movementDistance = zeros(nFrames,nPoints);
hasChanged = zeros(nPoints,1);
hasBeenOutOfContactArea = zeros(nPoints,1);
hasSlipped = zeros(nPoints,1);
fprintf('\nnFrames %d, ImLenght %d, imgHeight %d\n', nFrames,imageLength,imageHeight);
trackedImageArray = zeros(nFrames,imageLength,imageHeight,3);

% INITIALISE POINT TRACKER
fprintf('Initilize Point Tracker...');
pointTracker = vision.PointTracker('MaxBidirectionalError',artefactsThreshold);
POINT = [x' y'];
initialize(pointTracker, POINT, grayFrame0);
fprintf('Done\n');

% TRACKING LOOP 
if isPlot
    figure; set(gcf,'Color',[1 1 1]);
end

fprintf('Starting loop...\n');
for i = 1:nFrames

    fprintf('First If statement of %dth out of %d Frames', i, nFrames);
    if i == 1
        grayFrame = grayFrame0;
        filteredFrame = filteredFrame0;
        bwFrame = bwFrame0;
        flatBorderFrame = flatBorderFrame0;
        contactAreaFrame = contactAreaFrame0;
        otsuLevelStruct = otsuLevelStruct0;
    else % i>1
        searchArea = imdilate(squeeze(contactAreaTrackingArray(i-1,:,:)),strel('disk',10));
        [grayFrame,filteredFrame,bwFrame,flatBorderFrame,contactAreaFrame,contactPixels,tempOtsuLevelStruct] = filterForContactArea(croppedVideoArray(i,:,:,:),A,otsuThreshold,pixelsThreshold,se,se2,otsuLevelStruct);%otsuLevel0);
        grayFrame = squeeze(grayFrame);
        filteredFrame = squeeze(filteredFrame);
        bwFrame = squeeze(bwFrame);
        flatBorderFrame = squeeze(flatBorderFrame);
        contactAreaFrame = squeeze(contactAreaFrame).*searchArea;
        if isFitEllipse
            contactAreaFrame = myFitEllipseToContact(contactAreaFrame);
        end
        if exist('otsuLevelStruct0','var') && isfield(otsuLevelStruct0,'initContactArea')
            contactAreaFrame = otsuLevelStruct0.initContactArea;
%             contactAreaFrame = contactAreaFrame0;
        end
        if isIgnoreMultipleContact
            contactAreaFrame = ignoreMultipleContact(contactAreaFrame,grayFrame,flatBorderFrame);
        end
    end   
    
    fprintf('Done\n');
     
    fprintf('Saving the Different Frames...');
    grayTrackingArray(i,:,:) = grayFrame;
    filteredTrackingArray(i,:,:) = filteredFrame;
    bwTrackingArray(i,:,:) = bwFrame;
    flatBorderTrackingArray(i,:,:) = flatBorderFrame;
    contactAreaTrackingArray(i,:,:) = contactAreaFrame;
    otsuLevelStructArray(i,1) = otsuLevelStruct;
    fprintf('Done\n');
    % Checks if the algorithm has detected any area
    fprintf('Check if contact Area is detected...\n');
    isContactAreaDetected = getIsContactAreaDetected(contactAreaFrame);
    fprintf('Done\n');
    %%  =============== CONTOUR DETECTION ========================
    
    fprintf('Contour detectionn...');
    % Create contour around contact area
    contourFrame = createContour(contactAreaFrame); 

    % Add contour to image
    videoFrame = mat2gray(filteredFrame); % HEBA % videoFrame = filteredFrame;
    videoFrame(contourFrame==1)=1;
    fprintf('Done\n');
    %%  ====================== TRACKING =========================

    % Tracks all the points in the gray frame(the first point is the marker)
    fprintf('Track All points in gray frame...');
    [pointLocs, pointValidity] = step(pointTracker, grayFrame);
    pointLocsPerFrame(i,:,:) = pointLocs;
    pointValidityPerFrame(i,:) = pointValidity;
    visiblePoints = pointLocs(pointValidity,:);
    MarkerDot = visiblePoints(1,:);
    fprintf('Done\n');
    
    
    % Determines when all the points are slipping
    slippingPoints = 0;
    trackedPointsInContact = 0;
    totalPoints = length(pointLocs);

    %%  ==================== TRACKING ==========================
    fprintf('Looping through %d points...', nPoints);
    markerSize = 3;
    for j = 2:nPoints
        try
            x = min(max(round(pointLocs(j,2)),1),imageLength);
            y = min(max(round(pointLocs(j,1)),1),imageHeight);
            isPointInContact = contactAreaFrame(x,y);
        catch
            keyboard
        end

        % If the tracked point is in contact and has always been in contact
        if  (isPointInContact == 1) && (hasBeenOutOfContactArea(j,1) == 0)
           % fprintf('Hiii\n');
            trackedPointsInContact = trackedPointsInContact + 1;

            % RED:          Marker
            % MAGENTA:      Slipped
            % BLUE:         Stuck
            % GREEN:        Lost

            if (pointValidity(j) == 1)
                
                xDistToMarker(i,j) = abs(pointLocs(j,1)-MarkerDot(1,1));
                yDistToMarker(i,j) = abs(pointLocs(j,2)-MarkerDot(1,2));
                if i == 1
                    videoFrame = insertMarker(videoFrame, pointLocs(j,:),'*','Color','blue','Size',markerSize);
                    continue;
                end
                xMovement = abs(xDistToMarker(i,j)-xDistToMarker(i-1,j));
                yMovement = abs(yDistToMarker(i,j)-yDistToMarker(i-1,j));

                movementDistance(i,j) = sqrt((xMovement^2) + (yMovement^2));

                % The distance between the marker and the point has changed
                if (movementDistance(i,j) >= distanceThreshold)
                    %fprintf('Moved too much\n');
                    videoFrame = insertMarker(videoFrame, pointLocs(j,:),'x','Color','magenta','Size',markerSize);
                    hasChanged(j,1) = 1;
                    hasSlipped(j,1) = 1;
                    slippingPoints = slippingPoints + 1;
                % The point has not moved
                else
                    %fprintf('Getting some points\n');
                    videoFrame = insertMarker(videoFrame, pointLocs(j,:),'*','Color','blue','Size',markerSize);
                end

            % The tracker lost the point
            else
                fprintf('Lost\n');
                videoFrame = insertMarker(videoFrame, pointLocs(j,:),'o','Color','green','Size',markerSize);
                hasSlipped(j,1) = 2;
            end
        else
            hasBeenOutOfContactArea(j,1) = 1;
        end  
    end

    % Marks the marker and updates it for comparison in the next frame
    videoFrame = insertMarker(videoFrame, MarkerDot,'Color','red','Size',10);
    trackedImageArray(i,:,:, :) = videoFrame;

    nTrackedPointsInContactArray(i) = trackedPointsInContact;
    nSlippingPointsArray(i) = slippingPoints;

    if i == 1
        tempFlatBorderFrame = flatBorderFrame;
    else %i > 1
        prevAreaContour = createContour(searchArea); 
        tempFlatBorderFrame = mat2gray(flatBorderFrame);
        tempFlatBorderFrame = repmat(tempFlatBorderFrame,1,1,3);
        for x = 1:size(tempFlatBorderFrame,1)
            for y = 1:size(tempFlatBorderFrame,2)
                if prevAreaContour(x,y) == 1
                    tempFlatBorderFrame(x,y,:) = [1 0 0];
                end
            end
        end
    end
    
    fprintf('Done Looping through points\n');
    if isPlot
        subplot(3,2,1); imshow(grayFrame); title('Gray Frame');
        subplot(3,2,3); imshow(uint8(double(filteredFrame)./double(max(max(filteredFrame))).*255)); title('Scaled Filtered Frame');
        subplot(3,2,5); imshow(bwFrame); title('BW Frame');
        subplot(3,2,2); imshow(tempFlatBorderFrame); title('Flat Border');
        subplot(3,2,4); imshow(contactAreaFrame); title('Contact Area');
        subplot(3,2,6); imshow(videoFrame); title('Tracking');
        keyboard%pause(0.1);
    end
end

release(pointTracker);

fprintf('Done.\n');

% lostPoints = 0;
% slippedPoints = 0;
% for i = 1:length(hasSlipped)
%     if hasSlipped(i) == 1
%         slippedPoints = slippedPoints + 1;
%     elseif hasSlipped(i) == 2
%         lostPoints = lostPoints + 1;
%     end
% end
% 
end


%% 
function [x,y] = selectPointsToTrack(contactAreaFrame0,grayFrame0)
    
    % Finds smallest rectangle around contact area
    fprintf('Selecting points to track...');
    s = regionprops(contactAreaFrame0,'basic');
    ROI = round(s.BoundingBox);
    
    % Finds good features to track
    eigenPoints = detectMinEigenFeatures(grayFrame0, 'ROI', ROI);

    % Saves the points coordinates in x and y
    trackedPoints = 0;
    for i = 1:length(eigenPoints)
        pointX = eigenPoints.Location(i,1);
        pointY = eigenPoints.Location(i,2);
        if contactAreaFrame0(round(pointY),round(pointX)) == 1
            trackedPoints = trackedPoints + 1;
            x(trackedPoints) = pointX;
            y(trackedPoints) = pointY;
        end
    end
    
    fprintf('Done selecting points to track...');

end


%%
function [x,y] = selectMarker(selectedX,selectedY,grayFrame0)

    % Finds the closest good-to-track feature to the selected point and sets it
    % as the marker point
    squareSide = 2;
    markerFound = 0;
    while markerFound == 0
        markerROIxmin = selectedX - squareSide/2;
        markerROIymin = selectedY - squareSide/2;
        markerROIwidth = squareSide;
        markerROIheight = squareSide;

        markerROI = round([markerROIxmin,markerROIymin,markerROIwidth,markerROIheight]);

        markerEigenPoints = detectMinEigenFeatures(grayFrame0, 'ROI', markerROI);

        numberOfMarkerEigenPoints = size(markerEigenPoints.Location,1);
        if numberOfMarkerEigenPoints == 0
            squareSide = 2*squareSide;
        else
            minDistanceToMarker = squareSide;
            for i = 1:numberOfMarkerEigenPoints
                eigenPointX = markerEigenPoints.Location(i,1);
                eigenPointY = markerEigenPoints.Location(i,2);

                % Calculates distance from the selected point to the
                % good-to-track point
                distanceToMarker = sqrt((selectedX-eigenPointX)^2 + (selectedY-eigenPointY)^2);

                % Refreshes closest good-to-track point
                if distanceToMarker < minDistanceToMarker
                    x = eigenPointX;
                    y = eigenPointY;
                end
            end
            markerFound = 1;
        end
    end
    
end