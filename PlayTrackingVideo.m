
savePath = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis';

subjectNum = 1;
stimOrderNum = 1;
dateNum = 131115;
fileName = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Stuff_ForwardsTrackingData.MAT';
distanceThreshold = 1;

trackingType = 1;                                                                                            
if stimOrderNum == 4
    trackingType = 3;
end
if trackingType == 1
    s2 = 'Backwards';
elseif trackingType == 2
    s2 = 'Forwards';
elseif trackingType == 3
    s2 = 'Plateau';
end

figure; set(gcf,'Color',[1 1 1]);
for stimInd=1:2%1:96
    
    %matlabRoute = sprintf('%s\\%s_%d_%sTrackingData.mat',savePath,fileName,stimInd,s2);% Route of the tracking data 
    fprintf('Loading File...');
    load(fileName);
    fprintf('Done\n');

    if trackingType == 3
        tempFx = [Fx(1), onFx(1), Fx(2), onFx(2), Fx(3)];
        tempFy = [Fy(1), onFy(1), Fy(2), onFy(2), Fy(3)];
        tempFz = [Fz(1), onFz(1), Fz(2), onFz(2), Fz(3)];
        Ft = sqrt(tempFx.*tempFx + tempFy.*tempFy);
        Fn = -tempFz;
        
        tempGrayTrackingArray = uint8(zeros(5,size(grayTrackingArray,2),size(grayTrackingArray,3)));
        tempGrayTrackingArray(1:2:5,:,:) = grayTrackingArray; tempGrayTrackingArray(2:2:4,:,:) = onGrayTrackingArray;
        grayTrackingArray = tempGrayTrackingArray;
        tempFilteredTrackingArray = uint8(zeros(5,size(grayTrackingArray,2),size(grayTrackingArray,3)));
        tempFilteredTrackingArray(1:2:5,:,:) = filteredTrackingArray; tempFilteredTrackingArray(2:2:4,:,:) = onFilteredTrackingArray;
        filteredTrackingArray = tempFilteredTrackingArray;
        tempBwTrackingArray = false(5,size(grayTrackingArray,2),size(grayTrackingArray,3));
        tempBwTrackingArray(1:2:5,:,:) = bwTrackingArray; tempBwTrackingArray(2:2:4,:,:) = onBwTrackingArray;
        bwTrackingArray = tempBwTrackingArray;
        tempFlatBorderTrackingArray = false(5,size(grayTrackingArray,2),size(grayTrackingArray,3));
        tempFlatBorderTrackingArray(1:2:5,:,:) = flatBorderTrackingArray; tempFlatBorderTrackingArray(2:2:4,:,:) = onFlatBorderTrackingArray;
        flatBorderTrackingArray = tempFlatBorderTrackingArray;
        tempContactAreaTrackingArray = false(5,size(grayTrackingArray,2),size(grayTrackingArray,3));
        tempContactAreaTrackingarray(1:2:5,:,:) = contactAreaTrackingArray; tempContactAreaTrackingArray(2:2:4,:,:) = onContactAreaTrackingArray;
        contactAreaTrackingArray = tempContactAreaTrackingArray;
        tempTrackingImageArray = zeros(5,size(trackingImageArray,2),size(trackingImageArray,3),size(trackingImageArray,4));
        tempTrackingImageArray(1:2:5,:,:,:) = trackingImageArray; tempTrackingImageArray(2:2:4,:,:,:) = onTrackingImageArray;
        trackingImageArray = tempTrackingImageArray;
    else
        Ft = sqrt(Fx.*Fx + Fy.*Fy);
        Fn = -Fz;
    end
    
    fprintf('Copying Some variables...');
    [xDistToMarker,yDistToMarker,movementDistance,isInContact,hasChanged,hasBeenOutOfContactArea,hasSlipped,videoArray] = ...
        measureSlip(grayTrackingArray,filteredTrackingArray,bwTrackingArray,flatBorderTrackingArray,contactAreaTrackingArray,trackingImageArray,pointLocsPerFrame,pointValidityPerFrame,distanceThreshold);
    
    imageLength = size(grayTrackingArray,2);
    imageHeight = size(grayTrackingArray,3);
    
    nFrames = size(grayTrackingArray,1);
    nPoints = size(pointLocsPerFrame,2);
    fprintf('Done\n');
    
%     xDistToMarker = zeros(nFrames,nPoints);
%     yDistToMarker = zeros(nFrames,nPoints);
%     movementDistance = zeros(nFrames,nPoints);
%     hasChanged = zeros(nPoints,1);
%     hasBeenOutOfContactArea = zeros(nPoints,1);
%     hasSlipped = zeros(nPoints,1);
    
    frameInd = 1;
    fprintf('Starting while loop...\n');
    while frameInd <= nFrames
        fprintf('%d/%d...\n', frameInd, nFrames);
        grayFrame = squeeze(grayTrackingArray(frameInd,:,:));
        filteredFrame = squeeze(filteredTrackingArray(frameInd,:,:));
        bwFrame = squeeze(bwTrackingArray(frameInd,:,:));
        flatBorderFrame = squeeze(flatBorderTrackingArray(frameInd,:,:));
        contactAreaFrame = squeeze(contactAreaTrackingArray(frameInd,:,:));
        trackingFrame = squeeze(videoArray(frameInd,:,:,:));

        subplot(3,4,1); imshow(grayFrame); title({sprintf('stimInd: %d',stimInd),sprintf('Frame #%d of %d',frameInd,nFrames),'Grayscale'});
        subplot(3,4,5); imshow(uint8(double(filteredFrame)./double(max(max(filteredFrame))).*255)); title('Filtered');
        subplot(3,4,9); imshow(bwFrame); title('B+W');
        subplot(3,4,3); imshow(flatBorderFrame); 
            if frameInd == nFrames
                title({'L-key = next stim, A-key = prev frame','Flat Border'});
            elseif frameInd == 1
                title({'L-key = next frame','Flat Border'});
            else
                title({'L-key = next frame, A-key = prev frame','Flat Border'});
            end
        subplot(3,4,7); imshow(contactAreaFrame); title('Contact Area');
        subplot(3,4,11); imshow(trackingFrame); title('Tracking');
        if frameInd > 1
            temp = abs(grayFrame - squeeze(grayTrackingArray(frameInd-1,:,:)));            
            subplot(3,4,2); imshow(uint8(double(temp)./double(max(max(temp))).*255)); title({sprintf('F_T %.2f, F_N %.2f',Ft(frameInd),Fn(frameInd)),'Diff Grayscale'});
            temp = abs(filteredFrame - squeeze(filteredTrackingArray(frameInd-1,:,:)));
            subplot(3,4,6); imshow(uint8(double(temp)./double(max(max(temp))).*255)); title('Diff Filtered');
            subplot(3,4,10); imshow(abs(bwFrame - squeeze(bwTrackingArray(frameInd-1,:,:)))); title('Diff B+W');
            subplot(3,4,4); imshow(abs(flatBorderFrame - squeeze(flatBorderTrackingArray(frameInd-1,:,:)))); title('Diff Flat Border');
            subplot(3,4,8); imshow(abs(contactAreaFrame - squeeze(contactAreaTrackingArray(frameInd-1,:,:)))); title('Diff ContactArea');
        else
            subplot(3,4,2); imshow(1); subplot(3,4,6); imshow(1); subplot(3,4,10); imshow(1); subplot(3,4,4); imshow(1); subplot(3,4,8); imshow(1);
        end
        subplot(3,4,12); plot(Ft); xlabel('Frame #'); ylabel('F_T');
        yLim = get(gca,'ylim');
        hold on; plot([frameInd, frameInd], [0, yLim(2)],'r-');
        hold off;
        fig = gcf;
        w = 0; 
        while w==0
            w = waitforbuttonpress;
            if w==1
                c = fig.CurrentCharacter;
                if c == 'l' || c == 'a'
                    break;
                end
            end
        end
        
        switch c
            case 'l'
                frameInd = frameInd + 1;
            case 'a'
                if frameInd > 1
                    frameInd = frameInd - 1;
                end
        end
    end

end