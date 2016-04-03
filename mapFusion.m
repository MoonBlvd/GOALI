function mapFusion(plotFig)
if ~exist('plotFig','var') || isempty(plotFig)
    plotFig = false;
end
[landmarkGT,landmarkLidar,landmarkCamera] = mapConstruction();
% Plot the landmarks

if plotFig
    plot(landmarkCamera(:,1),landmarkCamera(:,2),'r*');
    hold on;
    plot(landmarkLidar(:,1),landmarkLidar(:,2),'go');
    hold on;
    plot(landmarkGT(:,1),landmarkGT(:,2),'bp');
end
numCam = size(landmarkCamera,1);
minL = 1.5;% Minimum length between two different features.
for i = 1:numCam
    temp = (landmarkCamera-repmat(landmarkCamera(i,:),numCam,1)).^2;
    dist = sum(temp,2);
    dist(i) = [];
    indexCam = find(dist < minL)
end
end
