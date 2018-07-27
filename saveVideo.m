clc;
clear;

savePath = 'C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Video\';

videoReader = VideoReader('C:\Users\shyam\Desktop\MAtlab tests\drive-download-20160930T040850Z\Video\Video 4.MP4');
fps = videoReader.FrameRate;
fprintf('frame rate: %d\n', fps);
height = videoReader.Height;
fprintf('height: %d\n', height);
width = videoReader.Width;
fprintf('width: %d\n', width);
duration = videoReader.Duration;
fprintf('duration: %d\n', duration);


starts = [22, 34, 48, 60, 73, 85, 98, 112, 125, 137, 150];
ends = [28, 40, 53, 66, 78, 91, 103, 116, 124, 142, 155];


i = 1;
j = 1;
trial = 1;
frames = ((ends(trial)*120)-(starts(trial)*120))/10;
videoArray = zeros(frames,(0.8*height+1)-(0.4*height+1)+1,(0.775*width+1)-(0.15*width+1)+1,'uint8');
fprintf('Reading in Frames for Trial %d...', trial);
while hasFrame(videoReader) && i <= 18600
    
    thisFrame = readFrame(videoReader);
    if i >= starts(trial)*120 && i <= ends(trial)*120 && mod(i, 10) == 0
        videoArray(j,:,:) = thisFrame(round((0.4*height+1)):round((0.8*height+1)),round((0.15*width+1)):round((0.775*width+1)));
        j = j + 1;  
    end
    i = i+1;
    
    if i >= ends(trial)*120
        fprintf('Saving...');
        save(sprintf('%s12FPS_%d.mat',savePath, trial),'videoArray');
        fprintf('Done\n');
        trial = trial + 1;
        j = 1;
        frames = ((ends(trial)*120)-(starts(trial)*120))/10;
        fprintf('Reading in Frames for Trial %d...', trial);
        videoArray = zeros(21,(0.8*height+1)-(0.4*height+1)+1,(0.775*width+1)-(0.15*width+1)+1,'uint8');
    end
end
