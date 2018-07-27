function contactAreaFrame1 = ignoreMultipleContact(contactAreaFrame0,grayFrame0,flatBorderFrame0)

% if multiple contact regions, pick the one with the biggest contrast
s = regionprops(contactAreaFrame0,'all');
nRegion = length(s);
if nRegion == 0
    contactAreaFrame0 = flatBorderFrame0;
    s = regionprops(contactAreaFrame0,'all');
    nRegion = length(s);
    if nRegion == 0
        fprintf('Cannot detect a contact area\n');
        contactAreaFrame1 = contactAreaFrame0;
%         keyboard;
        return;
    end
end

%         if isPlotIntermediate
%             figure; set(gcf,'Color',[1 1 1]);
%             subplot((nRegion+1),3,1); imshow(grayFrame0);
%             subplot((nRegion+1),3,2); imshow(contactAreaFrame0);
%         end
for i = 1:nRegion
    pixIdx = s(i).PixelIdxList;
    thisRegionContactArea = s(i).FilledImage;
    thisRegionContactArea = false(size(contactAreaFrame0)); thisRegionContactArea(pixIdx) = true;   
    temp = double(grayFrame0); temp(~thisRegionContactArea)=nan;
    intensity(i) = nanmedian(temp(:));
%             if isPlotIntermediate
%                 subplot((nRegion+1),3,(i)*3+2); 
%                     imshow(thisRegionContactArea); 
%                     title(sprintf('Region %d Contact Area',i));
%                 subplot((nRegion+1),3,(i)*3+3);
%                     hist(temp(:)); 
%                     title(sprintf('Region %d Intensity Dist.',i));
%                     legend(sprintf('Median: %d',intensity(i)));
%             end
end
[maxVal,regionInd] = max(intensity);
pixIdx = s(regionInd).PixelIdxList;
contactAreaFrame1 = false(size(contactAreaFrame0)); 
contactAreaFrame1(pixIdx) = true;
