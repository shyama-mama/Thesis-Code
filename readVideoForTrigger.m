function [trigSig,fps] = readVideoForTrigger(fullFilename)

fprintf('Retrieving video trigger signal...');

videoReader = VideoReader(fullFilename);
fps = videoReader.FrameRate;
height = videoReader.Height;
width = videoReader.Width;
duration = videoReader.Duration;
% bitsPerPixel = videoReader.BitsPerPixel;
% format = videoReader.VideoFormat;

comparisonPercent = 10;
incrementPercent = 10;

i = 1;
% meanIm = zeros(height,round(0.1*width));
% trigImArray = zeros(duration*fps,height,round(0.1*width),'uint8'); %uint8(zeros(duration*fps,height,round(0.2*width)));
while hasFrame(videoReader)

    
    % Selects the next frame
    thisFrame = readFrame(videoReader);%step(videoReader);
    
    thisTrigFrame = thisFrame(:,(0.9*width+1):end,:);
    thisTrigFrameGray = rgb2gray(thisTrigFrame);
%     trigImArray(i,:,:) = thisTrigFrameGray;
    trigSig(i,1) = mean2(thisTrigFrameGray);
%     if i < fps
%         meanIm = meanIm + double(thisTrigFrameGray);
%     elseif i == fps
%         meanIm = meanIm./fps;
%         trigThresh = mean2(meanIm).*1.5;
%     else
%         if mean2(thisTrigFrameGray) > trigThresh
%             keyboard
%         end
%     end
    
    % Update progress
    percentComplete = 100.*(i./fps)./duration;
    if percentComplete > comparisonPercent
        fprintf('%d%%...',comparisonPercent);
        comparisonPercent = comparisonPercent + incrementPercent;
    end
    
    i = i+1;
    
end

fprintf('Done.\n');
