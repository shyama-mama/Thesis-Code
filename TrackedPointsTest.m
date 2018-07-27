clc;
clear;

load('C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Triangles\PlateOff_Trial_1_Triangles.mat');
%load('C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\12FPS_PlateOff_TrackedPoints.mat');
SavePath = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\StrainValues\';

[F, T, z, c] = size(triCoord);
extratests = 0;

if extratests == 0
    for h = 1:F-1
        strainValues = zeros(F-1T);
        XYstrainValues = zeros(T);
        lineStrainValues = zeros(T,3);
        
        %strainValues = zeros(T);
        for j = 1:F
            for i = 1:T
                a = allTris(j,i,1);
                b = allTris(j,i,2);
                c = allTris(j,i,3);
                triCoord(j,i,A,1) = trackedPointsPerFrame(j+h,a,1);
                triCoord(j,i,A,2) = trackedPointsPerFrame(j+h,a,2);
                triCoord(j,i,B,1) = trackedPointsPerFrame(j+h,b,1);
                triCoord(j,i,B,2) = trackedPointsPerFrame(j+h,b,2);
                triCoord(j,i,C,1) = trackedPointsPerFrame(j+h,c,1);
                triCoord(j,i,C,2) = trackedPointsPerFrame(j+h,c,2);
                x = [trackedPointsPerFrame(j+h,a,1), trackedPointsPerFrame(j+h,b,1), trackedPointsPerFrame(j+h,c,1)];
                y = [trackedPointsPerFrame(j+h,a,2), trackedPointsPerFrame(j+h,b,2), trackedPointsPerFrame(j+h,c,2)];
                Xlengths(j,i,1) = abs(trackedPointsPerFrame(j+h,a,1) - trackedPointsPerFrame(j+h,b,1));
                Xlengths(j,i,2) = abs(trackedPointsPerFrame(j+h,a,1) - trackedPointsPerFrame(j+h,c,1));
                Xlengths(j,i,3) = abs(trackedPointsPerFrame(j+h,b,1) - trackedPointsPerFrame(j+h,c,1));
                p1 = [trackedPointsPerFrame(j+h,a,1), trackedPointsPerFrame(j+h,a,2)];
                p2 = [trackedPointsPerFrame(j+h,b,1), trackedPointsPerFrame(j+h,b,2)];
                p3 = [trackedPointsPerFrame(j+h,c,1), trackedPointsPerFrame(j+h,c,2)];
                AB = [p1;p2];
                AC = [p1;p3];
                BC = [p2;p3];
                XYLengths(j,i,1) = pdist(AB,'euclidean');
                XYLengths(j,i,2) = pdist(AC,'euclidean');
                XYLengths(j,i,3) = pdist(BC,'euclidean');             
                areas(j,i,1) = polyarea(x,y);
                centroids(j,i,1) = sum(x)/3;
                centroids(j,i,2) = sum(y)/3;
                
                if j == 2
                    ABStrain = ((Xlengths(2,i,1)-Xlengths(1,i,1))/Xlengths(1,i,1))*100;
                    ACStrain = ((Xlengths(2,i,2)-Xlengths(1,i,2))/Xlengths(1,i,2))*100;
                    BCStrain = ((Xlengths(2,i,3)-Xlengths(1,i,3))/Xlengths(1,i,3))*100;
                    lineStrainValues(i,1) = ABStrain;
                    lineStrainValues(i,2) = ACStrain;
                    lineStrainValues(i,3) = BCStrain;
                    AvgStrain = (ABStrain+ACStrain+BCStrain)/3;
                    strainValues(i) = AvgStrain;               
                    ABStrain = ((XYLengths(2,i,1)-XYLengths(1,i,1))/XYLengths(1,i,1))*100;
                    ACStrain = ((XYLengths(2,i,2)-XYLengths(1,i,2))/XYLengths(1,i,2))*100;
                    BCStrain = ((XYLengths(2,i,3)-XYLengths(1,i,3))/XYLengths(1,i,3))*100;
                    AvgStrain = (ABStrain+ACStrain+BCStrain)/3;
                    XYstrainValues(i) = AvgStrain;
                end
                
            end
        end
        
        for i=1:T
            overallXAvgStrain(h+1,i) = strainValues(i);
            overallXYAvgStrain(h+1,i) = XYstrainValues(i);
            overallLineStrain(h+1,i,1) = lineStrainValues(i,1);
            overallLineStrain(h+1,i,2) = lineStrainValues(i,2);
            overallLineStrain(h+1,i,3) = lineStrainValues(i,3);
        end
         
        D = zeros(T,1);
        XDist = zeros(T,1);
        G = zeros(T,1);
        for i = 1:10
            p1 = [centroids(1,i,1), centroids(1,i,2)];
            p2 = [centroids(2,i,1), centroids(2,i,2)];
            X = [p1 ; p2];
            D(i) = pdist(X, 'euclidean');
            G(i) = (centroids(1,i,2) - centroids(2,i,2)) / (centroids(1,i,1) - centroids(2,i,1));
        end  
       
        figure; 
        imshow(squeeze(grayTrackingArray(u,:,:))); 
        hold on
        clear s
        %highest = max(XYstrainValues);
        %lowest = min(XYstrainValues);
        %tresh = (lowest+highest)/2;
        highest = max(lineStrainValues(:));
        
        lowest = min(lineStrainValues(:));
        
        tresh = mean(lineStrainValues(:));
        
        upperbound = (highest+tresh)/2;
        lowerbound = (lowest+tresh)/2;
        fprintf('H%d L%d T%d\n', highest, lowest, tresh);
        for j = 1:T
            
            if lineStrainValues(j,1) < lowerbound
                line([triCoord(2,j,1,1),triCoord(2,j,2,1)], [triCoord(2,j,1,2),triCoord(2,j,2,2)], 'Color', 'r');
            elseif lineStrainValues(j,1) > upperbound
                line([triCoord(2,j,1,1),triCoord(2,j,2,1)], [triCoord(2,j,1,2),triCoord(2,j,2,2)], 'Color', 'b');
            else
                line([triCoord(2,j,1,1),triCoord(2,j,2,1)], [triCoord(2,j,1,2),triCoord(2,j,2,2)], 'Color', 'g');
            end
            
            if lineStrainValues(j,2) < tresh
                line([triCoord(2,j,1,1),triCoord(2,j,3,1)], [triCoord(2,j,1,2),triCoord(2,j,3,2)], 'Color', 'r');
            elseif lineStrainValues(j,2) > tresh
                line([triCoord(2,j,1,1),triCoord(2,j,3,1)], [triCoord(2,j,1,2),triCoord(2,j,3,2)], 'Color', 'b');
            else
                line([triCoord(2,j,1,1),triCoord(2,j,3,1)], [triCoord(2,j,1,2),triCoord(2,j,3,2)], 'Color', 'g');
            end
            
            if lineStrainValues(j,3) < tresh
                line([triCoord(2,j,2,1),triCoord(2,j,3,1)], [triCoord(2,j,2,2),triCoord(2,j,3,2)], 'Color', 'r');
            elseif lineStrainValues(j,3) > tresh
                line([triCoord(2,j,2,1),triCoord(2,j,3,1)], [triCoord(2,j,2,2),triCoord(2,j,3,2)], 'Color', 'b');
            else
                line([triCoord(2,j,2,1),triCoord(2,j,3,1)], [triCoord(2,j,2,2),triCoord(2,j,3,2)], 'Color', 'g');
            end
            
            %{
            s.Vertices = [triCoord(2,j,1,1) triCoord(2,j,1,2); triCoord(2,j,2,1) triCoord(2,j,2,2); triCoord(2,j,3,1) triCoord(2,j,3,2)];
            s.Faces = [1 2 3];
            %x = [triCoord(2,j,1,1), triCoord(2,j,2,1), triCoord(2,j,3,1)];
            %y = [triCoord(2,j,1,2), triCoord(2,j,2,2), triCoord(2,j,3,2)];
            %if (areas(2,j,1) - areas(1,j,1))/areas(1,j,1) < 0
            if XYstrainValues(j) < tresh
            %if centroids(2,j,1) < centroids(1,j,1) 
                %patch(x,y,'red')
                s.EdgeColor = 'red';
                s.FaceColor = 'red';
            %elseif (areas(2,j,1) - areas(1,j,1))/areas(1,j,1) > 0
            %elseif centroids(2,j,1) > centroids(1,j,1)
            elseif XYstrainValues(j) > tresh 
            %patch(x,y,'blue');
                s.EdgeColor = 'blue';
                s.FaceColor = 'blue';
            else
                %patch(x,y,'green');
                s.EdgeColor = 'green';
                s.FaceColor = 'green';
            end
            patch(s);
            
            %}
        end
        save = sprintf('%s%d', saveFilename, h);
        print(save,'-dpng');
        hold off;
        %}

    end
    
    save(sprintf('%sStrainValuesTrial1.mat',SavePath),'overallXAvgStrain', 'overallXYAvgStrain', 'overallLineStrain');
    
    
   
else
    
    %{
    for h = 25:61
        
        %{
        fprintf('Loading Filter Play Back...\n');
        %{
        imshow(squeeze(grayTrackingArray(h,:,:,:)));
        path = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Testing Points\';
        savepath = sprintf('%s%d.jpg', path, h);
        imwrite(squeeze(grayTrackingArray(h,:,:,:)), savepath);
        I = imread(savepath);
        marked = insertMarker(I, [pointLocsPerFrame(h,1,1) pointLocsPerFrame(h,1,2)]);
        %}
        y = 1;
        g = 2;
        fprintf('The co-ordinates for frame %d are\n', h);
        while g <= 1761 && y < 10               
            if hasBeenOutOfContactArea(g) == 0
                %{
                marked = insertMarker(marked, [pointLocsPerFrame(h,g,1) pointLocsPerFrame(h,g,2)]); 
                imshow(squeeze(marked));
                trackedPointsPerFrame(h,y,1) = pointLocsPerFrame(h,g,1);
                trackedPointsPerFrame(h,y,2) = pointLocsPerFrame(h,g,2);
                %}
                fprintf('%d and %d\n', pointLocsPerFrame(h,g,1), pointLocsPerFrame(h,g,2)); 
                y = y + 1;

            end
            g = g+1;
        end
        fprintf('\n');
        %}
        
        
        figure;
        imshow(squeeze(grayTrackingArray(h,:,:,:)));
        hold on;
        P = 1761;
        array = zeros(P);
        for p = 1:P
            if hasBeenOutOfContactArea(p) == 0
                value = ((movementDistance(h+1,p) - movementDistance(h,p))/movementDistance(h,p))*100;
                array(p) = value;
            end
        end
        
        lowest = min(array);
        highest = max(array);
        tresh = (lowest+highest)/2;
        
        for p = 1:P
            if hasBeenOutOfContactArea(p) == 0
                value = ((movementDistance(h+1,p) - movementDistance(h,p))/movementDistance(h,p))*100;
                if value < tresh
                    scatter(pointLocsPerFrame(h+1,p,1),pointLocsPerFrame(h+1,p,2), 'filled', 'red'); 
                elseif  value > tresh
                    scatter(pointLocsPerFrame(h+1,p,1),pointLocsPerFrame(h+1,p,2), 'filled', 'blue'); 
                else
                    scatter(pointLocsPerFrame(h+1,p,1),pointLocsPerFrame(h+1,p,2), 'filled', 'green'); 
                end
            end
        end
        save = sprintf('%s%d', saveFilename, h);
        print(save,'-dpng');
        hold off;
        
    end
    %}
    
end


