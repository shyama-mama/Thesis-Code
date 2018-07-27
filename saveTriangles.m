clc;
clear;

load('C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Tracked Points\10_12FPS_PlateOn_TrackedPoints.mat');
SavePath = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Triangles\';


[frames, points] = size(pointLocsPerFrame(1:end,1:end));
points = points/2;
A = 1;
B = 2;
C = 3;

for i = 1:frames
    imshow(squeeze(trackingImageArray(i,:,:,:)));
end

X = 1;
Y = 2;

markerPosition = zeros(frames,2);

for i = 1:frames
    pointsInContactArea = 0;
    for j = 2:points
        if hasBeenOutOfContactArea(j) == 0
            pointsInContactArea = pointsInContactArea + 1;
            trackedPointsPerFrame(i,pointsInContactArea,X) = pointLocsPerFrame(i,j,X);
            trackedPointsPerFrame(i,pointsInContactArea,Y) = pointLocsPerFrame(i,j,Y);       
        end
    end
    
    markerPosition(i,X) = pointLocsPerFrame(i,1,X);
    markerPosition(i,Y) = pointLocsPerFrame(i,1,Y);
    
end

points = pointsInContactArea;

for i = 1:frames
    x = trackedPointsPerFrame(i,1:points,X);
    y = trackedPointsPerFrame(i,1:points,Y);
    dt = delaunay(x', y');
   
    if i == 1
        triangles = size(dt,1);
    end
    
    numTris = size(dt,1);
    %fprintf('Frame %d and Triangles %d\n', i, numTris);
    for j = 1:numTris
        allTris(i,j,A) = dt(j,A);
        allTris(i,j,B) = dt(j,B);
        allTris(i,j,C) = dt(j,C);
    end
    
end  

copyAllTris = allTris;


found = 0;
trisNotFound = zeros(frames, round(triangles/2), 1);
for k = 2:frames
    counter = 0;
    for i = 1:triangles
        a = copyAllTris(1,i,A);
        b = copyAllTris(1,i,B);
        c = copyAllTris(1,i,C);
        for j = 1:triangles
            if sum(find(copyAllTris(k,j,:)== a))~=0
                if sum(find(copyAllTris(k,j,:)== b))~=0
                    if sum(find(copyAllTris(k,j,:)== c))~=0
                        found = 1;
                        allTris(k,i,A) = a;
                        allTris(k,i,B) = b;
                        allTris(k,i,C) = c;
                    end
                end
            end
        end
        if found == 1
            found = 0;
        else 
            counter = counter+1;
            trisNotFound(k,counter) = i;
        end
    end
end

triCoord = zeros(frames, triangles, 3, 2);
Xlengths = zeros(frames,triangles,3);
XYLengths = zeros(frames,triangles,3);
centroids = zeros(frames,triangles,2);
areas = zeros(frames,triangles,1);
angleValues = zeros(frames,triangles,3);

for j = 1:frames
    for i = 1:triangles
        a = allTris(j,i,A);
        b = allTris(j,i,B);
        c = allTris(j,i,C);
        
        triCoord(j,i,A,X) = trackedPointsPerFrame(j,a,X);
        triCoord(j,i,A,Y) = trackedPointsPerFrame(j,a,Y);
        triCoord(j,i,B,X) = trackedPointsPerFrame(j,b,X);
        triCoord(j,i,B,Y) = trackedPointsPerFrame(j,b,Y);
        triCoord(j,i,C,X) = trackedPointsPerFrame(j,c,X);
        triCoord(j,i,C,Y) = trackedPointsPerFrame(j,c,Y);
        
        x = [trackedPointsPerFrame(j,a,1), trackedPointsPerFrame(j,b,1), trackedPointsPerFrame(j,c,1)];
        y = [trackedPointsPerFrame(j,a,2), trackedPointsPerFrame(j,b,2), trackedPointsPerFrame(j,c,2)];
        areas(j,i,1) = polyarea(x,y);
        
        Xlengths(j,i,1) = abs(trackedPointsPerFrame(j,a,1) - trackedPointsPerFrame(j,b,1));
        Xlengths(j,i,2) = abs(trackedPointsPerFrame(j,a,1) - trackedPointsPerFrame(j,c,1));
        Xlengths(j,i,3) = abs(trackedPointsPerFrame(j,b,1) - trackedPointsPerFrame(j,c,1));
        
        p1 = [trackedPointsPerFrame(j,a,1), trackedPointsPerFrame(j,a,2)];
        p2 = [trackedPointsPerFrame(j,b,1), trackedPointsPerFrame(j,b,2)];
        p3 = [trackedPointsPerFrame(j,c,1), trackedPointsPerFrame(j,c,2)];
        
        angleValues(j,i,1) = atan2(2*areas(j,i),dot(p2-p1,p3-p1));
        angleValues(j,i,2) = atan2(2*areas(j,i),dot(p1-p2,p3-p2));
        angleValues(j,i,3) = atan2(2*areas(j,i),dot(p2-p3,p1-p3));
        
        AB = [p1;p2];
        AC = [p1;p3];
        BC = [p2;p3];
        
        XYLengths(j,i,1) = pdist(AB,'euclidean');
        XYLengths(j,i,2) = pdist(AC,'euclidean');
        XYLengths(j,i,3) = pdist(BC,'euclidean');      

        centroids(j,i,1) = sum(x)/3;
        centroids(j,i,2) = sum(y)/3;
    end
end

save(sprintf('%s10_PlateOn_Triangles_WithMarker.mat',SavePath),'triCoord', 'XYLengths', 'Xlengths', 'centroids', 'areas', 'grayTrackingArray', 'angleValues', 'trisNotFound', 'markerPosition');




