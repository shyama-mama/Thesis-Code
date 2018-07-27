function [ellipseContactFrame] = myFitEllipseToContact(contactAreaFrame)

isPlot = false;

temp1 = cumsum(contactAreaFrame,1)>0; 
temp2 = flipud(cumsum(flipud(contactAreaFrame),1)>0);
temp3 = cumsum(contactAreaFrame,2)>0; 
temp4 = fliplr(cumsum(fliplr(contactAreaFrame),2)>0);
connected = temp1 & temp2 & temp3 & temp4; 

if isPlot
    figure;
    subplot(2,3,1); imshow(temp1);
    subplot(2,3,2); imshow(temp2);
    subplot(2,3,4); imshow(temp3);
    subplot(2,3,5); imshow(temp4);
    subplot(2,3,3); imshow(contactAreaFrame);
    subplot(2,3,6); imshow(connected);
end

s = regionprops(connected,'all');
b = s.BoundingBox;

imToFit1 = false(size(contactAreaFrame)); 
try
    imToFit1(:,max(floor(b(1)),1):min(ceil((b(1)+b(3))),size(contactAreaFrame,2))) = true; 
catch
    keyboard
end
imToFit2 = false(size(contactAreaFrame));
try
    imToFit2(max(floor(b(2)),1):min(ceil((b(2)+b(4))),size(contactAreaFrame,1)),:) = true; 
catch
    keyboard
end
imToFit = imToFit1 & imToFit2; 

if isPlot
    figure; imshow(imToFit);
end


s = regionprops(imToFit,'all');
[x,y] = addEllipseToPlot(s.Centroid,s.MajorAxisLength,s.MinorAxisLength,s.Orientation,isPlot);
ellipseIm = false(size(contactAreaFrame)); 
for i = 1:length(x) 
    xi = min(max(round(y(i)),1),size(ellipseIm,1));
    yi = min(max(round(x(i)),1),size(ellipseIm,2));
    ellipseIm(xi,yi) = true; 
end
ellipseContactFrame = imfill(ellipseIm,'holes');
if isPlot
    figure; imshow(ellipseContactFrame); 
end