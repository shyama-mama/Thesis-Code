clc;
clear;

    
load('C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Video\12FPS_11.mat');
SavePath = ('C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\');
thisFrameArray = videoArray;%(5:end,:,:);

height = (0.8*1080+1)-(0.4*1080+1)+1;
width = (0.775*1920+1)-(0.15*1920+1)+1;
frames = 21;

%Contact Area Params
fprintf('Initiating Play Back...\n');
[frames, m, n]  = size(thisFrameArray);
for i = 1:frames
    imshow(squeeze(thisFrameArray(i,:,:)));
end

se = strel('disk',8);% Structure element for top-hat filter
se2 = strel('disk', 5);% Structure element for the closing area
A =[20 20];% Matrix to make the border flat
otsuThreshold = 0.015;%0.025;% Otsu level threshold
pixelsThreshold = 1000; %10000;% Threshold in number of pixels to delete minor areas

% Point Tracking Params
artefactsThreshold = 0.5;% Maximum bidirectional error in pixels
distanceThreshold = 1;% Distance threshold between marker and tracked point to determine slip

% Crop Image
fprintf('Cropping Image...');
figure; imshow(squeeze(thisFrameArray(1,:,:)));
    title('Please, select the top left corner of cropped image');
    [x1,y1] = ginput(1); cropxmin = uint16(round(x1)); cropymin = uint16(round(y1));
    title('Please, select the bottom right corner of cropped image');
    [x2,y2] = ginput(1); cropxmax = uint16(round(x2)); cropymax = uint16(round(y2));
    
fprintf('xmin: %d, xmax = %d, ymin = %d, ymax = %d\n', cropxmin, cropxmax, cropymin, cropymax);

croppedArray = thisFrameArray(:,cropymin:cropymax,cropxmin:cropxmax);
fprintf('Done\n');

fprintf('Getting Marker...');
figure; imshow(squeeze(croppedArray(1,:,:)));
    title('Please, select the marker');
    [x_marker,y_marker] = ginput(1);
fprintf('Done\n');


h1 = cropymax-cropymin+1;
w1 = cropxmax-cropxmin+1;

h2 = sprintf('%d',h1);
w2 = sprintf('%d', w1);

height = str2num(h2);
width = str2num(w2);

[grayTrackingArray,filteredTrackingArray,bwTrackingArray,flatBorderTrackingArray,contactAreaTrackingArray,otsuLevelStructArray,...
        pointLocsPerFrame,pointValidityPerFrame,xDistToMarker,yDistToMarker,movementDistance,hasChanged,hasSlipped,hasBeenOutOfContactArea,trackingImageArray]...
        =trackPoints(artefactsThreshold,distanceThreshold,thisFrameArray(:,:,:),x_marker,y_marker,A,otsuThreshold,pixelsThreshold,se,se2,true,false);
    %{
for i = 1:21
    figure; 
    hold on;
    
    subplot(imshow(squeeze(grayTrackingArray(i,:,:,:))));
    subplot(imshow(squeeze(filteredTrackingArray(i,:,:,:))));
    subplot(imshow(squeeze(bwTrackingArray(i,:,:,:))));
    subplot(imshow(squeeze(flatBorderTrackingArray(i,:,:,:))));
    subplot(imshow(squeeze(contactAreaTrackingArray(i,:,:,:))));
    subplot(imshow(squeeze(otsuLevelStructArray(i,:,:,:))));
    hold off;
end
        %}

save(sprintf('%s12FPS_PlateOff_TrackedPoints_11.mat',SavePath),'trackingImageArray', 'pointLocsPerFrame', 'pointValidityPerFrame', 'xDistToMarker', 'yDistToMarker', 'movementDistance', 'hasChanged', 'hasSlipped', 'hasBeenOutOfContactArea', 'grayTrackingArray', 'filteredTrackingArray', 'bwTrackingArray', 'flatBorderTrackingArray', 'contactAreaTrackingArray', 'otsuLevelStructArray');
