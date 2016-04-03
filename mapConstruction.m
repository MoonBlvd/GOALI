function [landmarkGT,landmarkLidar,landmarkCamera] = mapConstruction()
clear
if ~exist('plotFig','var') || isempty(plotFig)
    plotFig = false;
end
%load aa3_lsr2.mat
% for i = 1:size(LASER,1)
%     RR = double(LASER(i,:))/100;
%     x = detectTreesI16(RR,i)
% end
%% Get and save the x,y coordinates of landmarks based on lidar
% global State;
% numSteps = 500;
% run(numSteps,'vp',0.01,'nn');
% landmarkX = State.Ekf.mu(4:2:end-1);
% landmarkY = State.Ekf.mu(5:2:end);
% plot(landmarkX,landmarkY,'r*');
% landmarkLidar = [landmarkX landmarkY];
% save('landmarkLidar.mat','landmarkLidar');

%% Get the Lidar,ground truth, and camera x,y coordinates of landmarks
numL1 = randi([1,30],1,randi([1,10],1,1)); % number of the redundant observations of Lidar
numL2 = randi([1,30],1,randi([1,10],1,1)); % number of the missed observations of Lidar
numC1 = randi([1,30],1,randi([1,10],1,1)); % number of the redundant observations of Camera
numC2 = randi([1,30],1,randi([1,10],1,1)); % number of the missed observation of Camera
load ('landmarkLidar.mat');
landmarkGT = round(landmarkLidar);

% Build Lidar coordinates 
% id = randi([1,size(landmarkLidar,1)],1,length(numL2));
landmarkLidar(numL2,:) = [];
% j = 0;
for i = 1: length(numL1)
    index = numL1(i);
    ratio = [normrnd(1,0.01),normrnd(1,0.01)];
%     if index > j
    landmarkLidar = [landmarkLidar(1:index,:);
                    landmarkLidar(index,:).*ratio;
                    landmarkLidar(index+1:end,:)];
%     else
%         landmarkLidar = [landmarkLidar(1:index,:);
%                         landmarkLidar(index,:).*ratio;
%                         landmarkLidar(index+1:end,:)];
%     end
%     j = index;

end

% Build Camera coordinates 
% Manually introduce Gaussian noise to the camera data, 
%   with mu = 0, sigma = 1;
camNoise = normrnd(0,1,size(landmarkGT,1),size(landmarkGT,2));
landmarkCamera = landmarkGT+camNoise;
landmarkCamera(numC2,:) = [];
for i = 1: length(numC1)
    index = numC1(i);
    ratio = [normrnd(1,0.01),normrnd(1,0.01)];
    landmarkCamera = [landmarkCamera(1:index,:);
                    landmarkCamera(index,:).*ratio;
                    landmarkCamera(index+1:end,:)];
end
end
