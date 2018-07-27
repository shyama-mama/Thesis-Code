function [croppedVideoArray] = cropVideo(videoArray,cropxmin,cropxmax,cropymin,cropymax)%newLength,newHeight,xShift,yShift)

fprintf('Cropping video...');

% xraw = size(videoArray,2);
% yraw = size(videoArray,3);
% 
% cropxmin = (xraw/2) - (newHeight/2) + xShift;
% cropxmax = (xraw/2) + (newHeight/2) + xShift;
% cropymin = (yraw/2) - (newLength/2) + yShift;
% cropymax = (yraw/2) + (newLength/2) + yShift;

croppedVideoArray = videoArray(:,cropxmin:cropxmax,cropymin:cropymax,:);

fprintf('Done.\n');
