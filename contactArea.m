clc;
clear;
load('C:\Users\shyam\Desktop\Past Sem\Thesis\MAtlab tests\drive-download-20160930T040850Z\Analysis\Tracked Points\9_12FPS_PlateOff_TrackedPoints.mat');

[F, H, W, Z] = size(trackingImageArray);
%fprintf('%d %d %d %d\n', F, H, W, Z);
AREA = zeros(F);
saveFilename = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Analysis\Triangles\AreaStress\';

for f = 1:F
    %{
    figure;
    imshow(squeeze(contactAreaTrackingArray(f,100:end,450:900)));
    hold on;
    st = regionprops(squeeze(contactAreaTrackingArray(f,100:end,450:900)), 'BoundingBox' );
    rectangle('Position',[st.BoundingBox],'EdgeColor','r','LineWidth',2 );
    hold off;
    %}
    AREA(f) = (nnz(contactAreaTrackingArray(f,100:end,450:900)==1));
    
    %{
    figure; 
    imshow(squeeze(BW_filled));
    hold on;
    for k=1:size(boundaries)
       b = boundaries{k};
       plot(b(:,2),b(:,1),'g','LineWidth',3);
    end
    hold off;
    %}
  
end

figure; imshow(squeeze(trackingImageArray(1,:,:,:)));

%{
plot(AREA);
xlabel('Frame Number') % x-axis label
ylabel('Contact Area') % y-axis label
title('Contact Area VS Frames');
save = sprintf('%sContactAreaTrial10', saveFilename);
print(save,'-dpng');
%}

