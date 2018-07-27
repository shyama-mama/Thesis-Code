function [xDistToMarker,yDistToMarker,movementDistance,isInContact,hasChanged,hasBeenOutOfContactArea,hasSlipped,videoArray] =  ...
    measureSlip(grayTrackingArray,filteredTrackingArray,bwTrackingArray,flatBorderTrackingArray,contactAreaTrackingArray,trackingImageArray,pointLocsPerFrame,pointValidityPerFrame,distanceThreshold)

imageLength = size(grayTrackingArray,2);
imageHeight = size(grayTrackingArray,3);

nFrames = size(grayTrackingArray,1);
nPoints = size(pointLocsPerFrame,2);

videoArray = zeros(nFrames,imageLength,imageHeight,3);

xDistToMarker = zeros(nFrames,nPoints);
yDistToMarker = zeros(nFrames,nPoints);
movementDistance = zeros(nFrames,nPoints);
isInContact = false(nFrames,nPoints);
hasChanged = zeros(nPoints,1);
hasBeenOutOfContactArea = zeros(nPoints,1);
hasSlipped = zeros(nPoints,1);

fprintf('\nMeasuring Slip...\n');
for frameInd = 1:nFrames
    fprintf('Measuring for frame %d of %d\n', frameInd, nFrames);
    fprintf('\tSqueezing: Gray Frame...');
    grayFrame = squeeze(grayTrackingArray(frameInd,:,:));
    fprintf('FilteredFrame...');
    filteredFrame = squeeze(filteredTrackingArray(frameInd,:,:));
    fprintf('BWFrame...');
    bwFrame = squeeze(bwTrackingArray(frameInd,:,:));
    fprintf('Flat Border Frame...');
    flatBorderFrame = squeeze(flatBorderTrackingArray(frameInd,:,:));
    fprintf('Contact Area Frame...');
    contactAreaFrame = squeeze(contactAreaTrackingArray(frameInd,:,:));
    fprintf('Tracking Frame...');
    trackingFrame = squeeze(trackingImageArray(frameInd,:,:,:));
    fprintf('Point locations...');
    pointLocs = squeeze(pointLocsPerFrame(frameInd,:,:));
    fprintf('point valdity...');
    pointValidity = squeeze(pointValidityPerFrame(frameInd,:));
    fprintf('Done Squeezing\n');
    
    fprintf('\tCreating Contour...');
    contourFrame = createContour(contactAreaFrame); 
    videoFrame = mat2gray(grayFrame); % HEBA % videoFrame = filteredFrame;
    videoFrame(contourFrame==1)=1;
    videoFrame = repmat(videoFrame,1,1,3);
    MarkerDot = pointLocs(1,:);
    markerSize = 1;
    fprintf('Done\n');
    
    fprintf('\tLooping Through %d Points...', nPoints);
    for pointInd = 2:nPoints
        try
            x = min(max(round(pointLocs(pointInd,2)),1),imageLength);
            y = min(max(round(pointLocs(pointInd,1)),1),imageHeight);
            isPointInContact = contactAreaFrame(x,y);
            isInContact(frameInd,pointInd) = isPointInContact;
        catch
            keyboard
        end

        % If the tracked point is in contact and has always been in contact
        if  (isPointInContact == 1) %&& (hasBeenOutOfContactArea(pointInd,1) == 0)

            % RED:          Marker
            % MAGENTA:      Slipped
            % BLUE:         Stuck
            % GREEN:        Lost

            if (pointValidity(pointInd) == 1)

                xDistToMarker(frameInd,pointInd) = abs(pointLocs(pointInd,1)-MarkerDot(1,1));
                yDistToMarker(frameInd,pointInd) = abs(pointLocs(pointInd,2)-MarkerDot(1,2));

                if frameInd == 1
                    videoFrame = insertMarker(videoFrame, pointLocs(pointInd,:),'*','Color','blue','Size',markerSize);
                    continue;
                end
                xMovement = abs(xDistToMarker(frameInd,pointInd)-xDistToMarker(1,pointInd));
                yMovement = abs(yDistToMarker(frameInd,pointInd)-yDistToMarker(1,pointInd));

                movementDistance(frameInd,pointInd) = sqrt((xMovement^2) + (yMovement^2));

                % The distance between the marker and the point has changed
                if (movementDistance(frameInd,pointInd) >= distanceThreshold)
                    videoFrame = insertMarker(videoFrame, pointLocs(pointInd,:),'x','Color','magenta','Size',markerSize);
                    hasChanged(pointInd,1) = 1;
                    hasSlipped(pointInd,1) = 1;
%                         slippingPoints = slippingPoints + 1;
                % The point has not moved
                else
                    videoFrame = insertMarker(videoFrame, pointLocs(pointInd,:),'*','Color','blue','Size',markerSize);
                end

            % The tracker lost the point
            else
                videoFrame = insertMarker(videoFrame, pointLocs(pointInd,:),'o','Color','green','Size',markerSize);
                hasSlipped(pointInd,1) = 2;
            end
        else
            hasBeenOutOfContactArea(pointInd,1) = 1;
        end  
    end
    fprintf('Done\n');
    fprintf('\tSaving...');
    videoArray(frameInd,:,:,:) = videoFrame;
    fprintf('Done\n\n');
    
end