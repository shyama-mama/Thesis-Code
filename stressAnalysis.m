clc;
clear;

load('C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Triangles\4_PlateOn_Triangles_WithMarker.mat');
load('C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Tracked Points\4_12FPS_PlateOn_TrackedPoints.mat');
saveFilename = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Triangles\AreaStress\';
[F, T, z, c] = size(triCoord);

areaStress = zeros(F, T);
XlineStrainValues = zeros(F,T,3);
angleStress = zeros(F,T,3);
centroidDistance = zeros(F,T);
markerDisplacement = zeros(F);

XAvgStrainValues = zeros(F,T);
XYAvgStrainValues = zeros(F,T);
%trisNotFound
for f = 2:F
    for t = 1:T
        if sum(find(trisNotFound(f,:,:)== t))==0
            areaStress(f,t) = ((areas(f,t)-areas(1,t))/areas(1,t))*100;
            ABStrain = ((Xlengths(f,t,1)-Xlengths(1,t,1))/Xlengths(1,t,1))*100;
            ACStrain = ((Xlengths(f,t,2)-Xlengths(1,t,2))/Xlengths(1,t,2))*100;
            BCStrain = ((Xlengths(f,t,3)-Xlengths(1,t,3))/Xlengths(1,t,3))*100;
            XlineStrainValues(f,t,1) = ABStrain;
            XlineStrainValues(f,t,2) = ACStrain;
            XlineStrainValues(f,t,3) = BCStrain;
            AvgStrain = (ABStrain+ACStrain+BCStrain)/3;
            XAvgStrainValues(f,t) = AvgStrain;               
            ABStrain = ((XYLengths(f,t,1)-XYLengths(1,t,1))/XYLengths(1,t,1))*100;
            ACStrain = ((XYLengths(f,t,2)-XYLengths(1,t,2))/XYLengths(1,t,2))*100;
            BCStrain = ((XYLengths(f,t,3)-XYLengths(1,t,3))/XYLengths(1,t,3))*100;
            AvgStrain = (ABStrain+ACStrain+BCStrain)/3;
            XYAvgStrainValues(f,t) = AvgStrain;
            
            p1 = [centroids(f,t,1), centroids(f,t,2)];
            p2 = [centroids(1,t,1), centroids(1,t,2)];
            displacement = [p1;p2];
            m1 = [markerPosition(f,1), markerPosition(f,2)];
            m2 = [markerPosition(1,1), markerPosition(1,2)];
            markerPoints = [m1;m2];
            markerDisplacement(f) = pdist(markerPoints,'euclidean');
            centroidDistance(f,t) = pdist(displacement,'euclidean');  
            centroidDistance(f,t) = centroidDistance(f,t);% - markerDisplacement(f);
            if centroids(f,t,1) - centroids(1,t,1) < 0
                centroidDistance(f,t) = centroidDistance(f,t)*-1;
            end
            if f ~= 1 
                angleStress(f,t,1) = ((angleValues(f,t,1)-angleValues(1,t,1))/angleValues(1,t,1))*100;
                angleStress(f,t,2) = ((angleValues(f,t,2)-angleValues(1,t,2))/angleValues(1,t,2))*100;
                angleStress(f,t,3) = ((angleValues(f,t,3)-angleValues(1,t,3))/angleValues(1,t,3))*100;
            end
        end
     end
end

isPlot = 1;
area = 0;
XYavg = 0;
XStrain = 0;
if isPlot 
    
    %{
    for f = 1:F
        figure; 
        hold on;
        h = histogram(centroidDistance(f,:));
        title('Centroid Distance');
        save = sprintf('%s%d', saveFilename, f);
        print(save,'-dpng');
        hold off;
    end
    %}
    h = histogram(markerDisplacement);
    
    
    
    figure; 
    h1 = histogram(areaStress(:,:));
    title('Area Histo');

    figure; 
    h2 = histogram(XAvgStrainValues(:,:));
    title('XAvgStrainValues');

    figure; 
    h3 = histogram(XYAvgStrainValues(:,:));
    title('XYAvgStrainValues');

    figure; 
    h4 = histogram(XlineStrainValues(:,:));
    title('Line');

    %figure; 
    %h5 = histogram(angleStress(25,:));
    %title('AngleStress');
    
    

    AREA = zeros(F,1);
    %fprintf('%d\n',AREA);

    for f = 1:F
        AREA(f) = (nnz(contactAreaTrackingArray(f,:,:)==1));  
    end

    %fprintf('%d\n', AREA);
    plot(AREA);
    %}
elseif area
    
    %fprintf('%d -> %d -> %d -> %d\n', mini, lowerbound, upperbound, maxi);
    for f = 30:F
        mini = (min(areaStress(f,:)));
        maxi = (max(areaStress(f,:)));
        avg = 0;
        
        lowerbound = (mini+avg)/2;
        upperbound = (maxi+avg)/2;
        lowermini = (lowerbound+mini)/2;
        uppermaxi = (lowerbound+maxi)/2;
        upperavg = (avg+lowerbound)/2;
        loweravg = (avg+upperbound)/2;
        figure; 
        imshow(squeeze(grayTrackingArray(f,:,:))); 
        hold on
        clear s
        for t = 1:T
            s.Vertices = [triCoord(f,t,1,1) triCoord(f,t,1,2); triCoord(f,t,2,1) triCoord(f,t,2,2); triCoord(f,t,3,1) triCoord(f,t,3,2)];
            s.Faces = [1 2 3];
            if (areaStress(f,t) == 0) || (areaStress(f,t) > avg && areaStress(f,t) < upperavg) || (areaStress(f,t) < avg && areaStress(f,t) > loweravg)
                s.EdgeColor = 'green';
                s.FaceColor = 'green';
            elseif areaStress(f,t) > uppermaxi
                s.EdgeColor = 'blue';
                s.FaceColor = 'blue';
            elseif areaStress(f,t) < upperbound && areaStress(f,t) > avg
                s.EdgeColor = 'green';
                s.FaceColor = 'green';
            elseif areaStress(f,t) > lowerbound && areaStress(f,t) < avg
                s.EdgeColor = 'green';
                s.FaceColor = 'green';
            elseif areaStress(f,t) < lowerbound && areaStress(f,t) > lowermini
                s.EdgeColor = 'yellow';
                s.FaceColor = 'yellow';
            elseif areaStress(f,t) > upperbound && areaStress(f,t) < uppermaxi
                s.EdgeColor = 'cyan';
                s.FaceColor = 'cyan';
            else
                s.EdgeColor = 'red';
                s.FaceColor = 'red';
            end
            patch(s);
        end
        save = sprintf('%s%d', saveFilename, f);
        print(save,'-dpng');
        hold off;
    end
elseif XStrain
        
    for f = 17:F
       
        figure; 
        imshow(squeeze(grayTrackingArray(f,:,:))); 
        hold on
        clear s
        for j = 1:T
             highest = max(XlineStrainValues(f,:));        
            lowest = min(XlineStrainValues(f,:));
            tresh = mean(XlineStrainValues(f,:));

            upperbound = (highest+tresh)/2;
            lowerbound = (lowest+tresh)/2;
            upperavg = (upperbound+tresh)/2;
            loweravg = (lowerbound+tresh)/2;
            highupper = (highest+upperbound)/2;
            lowlow = (lowest+lowerbound)/2;
            
            if XlineStrainValues(f,j,1) < lowerbound
                line([triCoord(f,j,1,1),triCoord(f,j,2,1)], [triCoord(f,j,1,2),triCoord(f,j,2,2)], 'Color', 'r');
            elseif XlineStrainValues(f,j,1) > upperbound
                line([triCoord(f,j,1,1),triCoord(f,j,2,1)], [triCoord(f,j,1,2),triCoord(f,j,2,2)], 'Color', 'b');
            elseif XlineStrainValues(f,j,1) > lowerbound && XlineStrainValues(f,j,1) < tresh
                line([triCoord(f,j,1,1),triCoord(f,j,2,1)], [triCoord(f,j,1,2),triCoord(f,j,2,2)], 'Color', 'y');
            elseif XlineStrainValues(f,j,1) < upperbound && XlineStrainValues(f,j,1) > tresh
                line([triCoord(f,j,1,1),triCoord(f,j,2,1)], [triCoord(f,j,1,2),triCoord(f,j,2,2)], 'Color', 'c');
            elseif (XlineStrainValues(f,j,1) < upperavg && XlineStrainValues(f,j,1) > tresh) || (XlineStrainValues(f,j,1) > lowlow && XlineStrainValues(f,j,1) < tresh)
                line([triCoord(f,j,1,1),triCoord(f,j,2,1)], [triCoord(f,j,1,2),triCoord(f,j,2,2)], 'Color', 'g');
            end
            
            if XlineStrainValues(f,j,2) < lowerbound
                line([triCoord(f,j,1,1),triCoord(f,j,3,1)], [triCoord(f,j,1,2),triCoord(f,j,3,2)], 'Color', 'r');
            elseif XlineStrainValues(f,j,2) > upperbound
                line([triCoord(f,j,1,1),triCoord(f,j,3,1)], [triCoord(f,j,1,2),triCoord(f,j,3,2)], 'Color', 'b');
            elseif XlineStrainValues(f,j,2) > lowerbound && XlineStrainValues(f,j,2) < tresh
                line([triCoord(f,j,1,1),triCoord(f,j,3,1)], [triCoord(f,j,1,2),triCoord(f,j,3,2)], 'Color', 'y');
            elseif XlineStrainValues(f,j,2) < upperbound && XlineStrainValues(f,j,2) > tresh
                line([triCoord(f,j,1,1),triCoord(f,j,3,1)], [triCoord(f,j,1,2),triCoord(f,j,3,2)], 'Color', 'c');
            elseif (XlineStrainValues(f,j,2) < upperavg && XlineStrainValues(f,j,2) > tresh) || (XlineStrainValues(f,j,2) > lowlow && XlineStrainValues(f,j,2) < tresh)
                line([triCoord(f,j,1,1),triCoord(f,j,3,1)], [triCoord(f,j,1,2),triCoord(f,j,3,2)], 'Color', 'g');
            end
            
            
            if XlineStrainValues(f,j,3) < lowerbound
                line([triCoord(f,j,2,1),triCoord(f,j,3,1)], [triCoord(f,j,2,2),triCoord(f,j,3,2)], 'Color', 'r');
            elseif XlineStrainValues(f,j,3) > upperbound
                line([triCoord(f,j,2,1),triCoord(f,j,3,1)], [triCoord(f,j,2,2),triCoord(f,j,3,2)], 'Color', 'b');
            elseif XlineStrainValues(f,j,3) > lowerbound && XlineStrainValues(f,j,3) < tresh
                line([triCoord(f,j,2,1),triCoord(f,j,3,1)], [triCoord(f,j,2,2),triCoord(f,j,3,2)], 'Color', 'y');
            elseif XlineStrainValues(f,j,3) < upperbound && XlineStrainValues(f,j,3) > tresh
                line([triCoord(f,j,2,1),triCoord(f,j,3,1)], [triCoord(f,j,2,2),triCoord(f,j,3,2)], 'Color', 'c');
            elseif (XlineStrainValues(f,j,3) < upperavg && XlineStrainValues(f,j,3) > tresh) || (XlineStrainValues(f,j,3) > lowlow && XlineStrainValues(f,j,3) < tresh)
                line([triCoord(f,j,2,1),triCoord(f,j,3,1)], [triCoord(f,j,2,2),triCoord(f,j,3,2)], 'Color', 'g');
            end
        end
        save = sprintf('%s%d', saveFilename, f);
        print(save,'-dpng');
        hold off;
    end
elseif XYavg
    
    for f = 17:F
       
        figure; 
        imshow(squeeze(grayTrackingArray(f,:,:))); 
        hold on
        clear s
        for j = 1:T
            s.Vertices = [triCoord(f,j,1,1) triCoord(f,j,1,2); triCoord(f,j,2,1) triCoord(f,j,2,2); triCoord(f,j,3,1) triCoord(f,j,3,2)];
            s.Faces = [1 2 3];
            highest = max(XYAvgStrainValues(f,:));        
            lowest = min(XYAvgStrainValues(f,:));
            tresh = mean(XYAvgStrainValues(f,:));
            fprintf('%d %d %d\n', highest, lowest, tresh);
            upperbound = (highest+tresh)/2;
            lowerbound = (lowest+tresh)/2;
            upperavg = (upperbound+tresh)/2;
            loweravg = (lowerbound+tresh)/2;
            highupper = (highest+upperbound)/2;
            lowlow = (lowest+lowerbound)/2;
            
            if XYAvgStrainValues(f,j) < lowerbound
                s.EdgeColor = 'red';
                s.FaceColor = 'red';
            elseif XYAvgStrainValues(f,j) > upperbound
                s.EdgeColor = 'blue';
                s.FaceColor = 'blue';
            elseif XYAvgStrainValues(f,j) > lowerbound && XYAvgStrainValues(f,j) < tresh
                s.EdgeColor = 'yellow';
                s.FaceColor = 'yellow';
            elseif XYAvgStrainValues(f,j) < upperbound && XYAvgStrainValues(f,j) > tresh
                s.EdgeColor = 'cyan';
                s.FaceColor = 'cyan';
            elseif (XYAvgStrainValues(f,j) < upperavg && XYAvgStrainValues(f,j) > tresh) || (XYAvgStrainValues(f,j) > lowlow && XYAvgStrainValues(f,j) < tresh)
                s.EdgeColor = 'green';
                s.FaceColor = 'green';
            end
            patch(s);
        end
        save = sprintf('%s%d', saveFilename, f);
        print(save,'-dpng');
        hold off;
    end
else
        maxi = max(centroidDistance(:));
        mini = min(centroidDistance(:));
        avg = mean(centroidDistance(:));
        lowerbound = (mini+avg)/2;
        upperbound = (maxi+avg)/2;
        upperavg = (upperbound+avg)/2;
        loweravg = (lowerbound+avg)/2;
    
    for f = 30:F
        %{
        maxi = max(centroidDistance(f,:));
        mini = min(centroidDistance(f,:));
        avg = mean(centroidDistance(f,:));
        lowerbound = (mini+avg)/2;
        upperbound = (maxi+avg)/2;
        upperavg = (upperbound+avg)/2;
        loweravg = (lowerbound+avg)/2;
        %}
        figure; 
        imshow(squeeze(grayTrackingArray(f,:,:))); 
        hold on
        for t = 1:T
            if (centroidDistance(f,t) > avg && centroidDistance(f,t) < upperavg) || (centroidDistance(f,t) < avg && centroidDistance(f,t) > loweravg)
                scatter(centroids(f,t,1),centroids(f,t,2), 'filled', 'green'); 
            elseif centroidDistance(f,t) > upperavg && centroidDistance(f,t) < upperbound
                scatter(centroids(f,t,1),centroids(f,t,2), 'filled', 'cyan'); 
            elseif centroidDistance(f,t) > upperbound
                scatter(centroids(f,t,1),centroids(f,t,2), 'filled', 'blue'); 
            elseif centroidDistance(f,t) < loweravg && centroidDistance(f,t) > lowerbound
                scatter(centroids(f,t,1),centroids(f,t,2), 'filled', 'yellow'); 
            elseif centroidDistance(f,t) < lowerbound
                scatter(centroids(f,t,1),centroids(f,t,2), 'filled', 'red'); 
            end
        end
        save = sprintf('%s%d', saveFilename, f);
        print(save,'-dpng');
        hold off;
    end
    
end