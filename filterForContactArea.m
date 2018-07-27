function [grayImageArray,filteredImageArray,bwImageArray,flatBorderImageArray,contactLogicalImageArray,contactPixelsArray,otsuLevelStructArray] = filterForContactArea(croppedVideoArray,A,otsuThreshold,pixelsThreshold,se,se2,otsuLevelStruct)

% fprintf('Filtering video for contact area...');

% Get dimensions
nFrames = size(croppedVideoArray,1);
frameHeight = size(croppedVideoArray,2);
frameWidth = size(croppedVideoArray,3);

% Initialise arrays
grayImageArray = uint8(zeros(nFrames,frameHeight,frameWidth));
filteredImageArray = uint8(zeros(nFrames,frameHeight,frameWidth));
bwImageArray = false(nFrames,frameHeight,frameWidth);
flatBorderImageArray = false(nFrames,frameHeight,frameWidth);
contactLogicalImageArray = false(nFrames,frameHeight,frameWidth);
contactPixelsArray = zeros(nFrames,1);
% otsuLevelStructArray = zeros(nFrames,1);

comparisonPercent = 10;
incrementPercent = 10;

for frameNumber = 1:nFrames

    % Update progress
%     percentComplete = 100.*(frameNumber)./nFrames;
%     if percentComplete > comparisonPercent
%         fprintf('%d%%...',comparisonPercent);
%         comparisonPercent = comparisonPercent + incrementPercent;
%     end
    
    % Selects the next frame
    croppedFrame = squeeze(croppedVideoArray(frameNumber,:,:,:));  
    
    
    if size(croppedFrame,3) == 3
        grayFrame = rgb2gray(croppedFrame);
    else
        % already gray
        grayFrame = croppedFrame;
    end
    grayImageArray(frameNumber,:,:) = grayFrame;

    % Measures the size of the cropped frame
    [croppedHeight,croppedLength] = size(grayFrame);

    
    filteredFrame = imtophat(grayFrame,se);
    filteredImageArray(frameNumber,:,:) = filteredFrame;
    
    
    % Calculates brightness threshold with Otsu method and show on screen
%     if isfield(otsuLevelStruct,'levels')
%         levels = otsuLevelStruct.levels;
%     else
%         levels = multithresh(filteredFrame,10);
%     end
%     otsuLevelStructArray(frameNumber,1).levels = levels;
%     
%     % Transforms frame to black and white format using the Otsu threshold if
%     % it is big enough to guarantee there is contact. It generates a black
%     % image if there is no contact.
%     segFrame = imquantize(filteredFrame,unique(levels));
%     RGBFrame = label2rgb(segFrame);
%     tempGrayFrame = rgb2gray(RGBFrame);
%         
%     if isfield(otsuLevelStruct,'level')
%         level = otsuLevelStruct.level;
%     else
%         level = graythresh(tempGrayFrame);
%     end
%     otsuLevelStructArray(frameNumber,1).level = level;
%     
%     BWFrame = im2bw(real(tempGrayFrame),level);
    
    if isfield(otsuLevelStruct,'level')
        level = otsuLevelStruct.level;
    else
        level = graythresh(filteredFrame);
    end
    otsuLevelStructArray(frameNumber,1).level = level;
    if level < otsuThreshold
        BWFrame = false(croppedHeight,croppedLength);
    else
        BWFrame = im2bw(real(filteredFrame),level);
    end
    bwImageArray(frameNumber,:,:) = BWFrame;
        
    % Morphological method to close fill the holes in the contact area
    closedFrame = imclose(BWFrame,se2);
    openFrame = imopen(closedFrame, se2);

    
    % Median filter to soften the border
    flatBorderFrame = medfilt2(openFrame,A);
    flatBorderImageArray(frameNumber,:,:) = flatBorderFrame;
    
        
    % Removes residual areas
    contactAreaFrame = bwareaopen(flatBorderFrame,pixelsThreshold);
    contactLogicalImageArray(frameNumber,:,:) = contactAreaFrame; % HEBA %
    
        
    % Writes the number of pixels in contact for this frame
    contactPixelsArray(frameNumber) = sum(sum(contactAreaFrame));
    
end

% fprintf('Done.\n');
