function [contourFlatFrame] = createContour(flatBorderFrame)

[imageLength,imageHeight] = size(flatBorderFrame);

% Calculates the contour of the contact area
C = contourc(double(flatBorderFrame));

% Extracts the number of points in the contour
contourPoints = size(C,2);

% Creates an image with a white line on the contour
contourFlatFrame = zeros(imageLength,imageHeight);
if contourPoints > 0
    xcontour = C(1,2:C(2,1)+1);
    ycontour = C(2,2:C(2,1)+1);
    for contourPointInd = 1:length(xcontour)
        contourFlatFrame(round(ycontour(contourPointInd)),round(xcontour(contourPointInd))) = 1;
    end
end

end