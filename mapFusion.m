function mapFusion(plotFig)
if ~exist('plotFig','var') || isempty(plotFig)
    plotFig = false;
end
[landmarkGT,landmarkLidar,landmarkCamera] = mapConstruction();
% Plot the landmarks
originalLidar = landmarkLidar;
numCam = size(landmarkCamera,1);
numLidar = size(landmarkLidar,1);
minL = 2.25;% Minimum length between two different features.
newMap = [0 0];
numNewMap = 1;
for i = 1:numCam
    temp  =(newMap - repmat(landmarkCamera(i,:),numNewMap,1)).^2;
    dist = sum(temp,2);
    overlap = find(dist < minL);
    if overlap
        overlap
        i
        newMap(overlap,:) = mean([newMap(overlap,:);landmarkCamera(i,:)]);
    else
%     temp1 = (landmarkCamera-repmat(landmarkCamera(i,:),numCam,1)).^2;
%     dist1 = sum(temp1,2); % We may calculate Mahalanobis distance later
%     %dist1(i) = 1e10;
%     indexCam = find(dist1 < minL);
        temp2 = (landmarkLidar-repmat(landmarkCamera(i,:),numLidar,1)).^2;
        dist2 = sum(temp2,2); % We should calculate Mahalanobis distance later
        indexLidar = find(dist2 < minL);   

%     if indexCam
%         % Find the lidar observations that are close to the indexCamth
%         % camera observation
%         for j = 1:length(indexCam)
%             temp = (landmarkLidar-repmat(landmarkCamera(indexCam(j),:),numLidar,1)).^2;
%             dist = sum(temp,2); % We should calculate Mahalanobis distance later
%             indexLidar = find(dist < minL);            
%         end
%         x = mean([landmarkLidar(indexLidar,1);landmarkCam(indexCam,1)]);
%         y = mean([landmarkLidar(indexLidar,2);landmarkCam(indexCam,2)]);
%         newMap = [newMap;x,y];
%     end
        x = mean([landmarkLidar(indexLidar,1);landmarkCamera(i,1)]);
        y = mean([landmarkLidar(indexLidar,2);landmarkCamera(i,2)]);
        newMap = [newMap;x,y];
        landmarkLidar(indexLidar,:) = [];
        numLidar = size(landmarkLidar,1);
    end
    numNewMap = size(newMap,1);
%     plot(landmarkCamera(i,1),landmarkCamera(i,2),'r*');
%     plot(landmarkLidar(i,1),landmarkLidar(i,2),'go');

%     plot(x,y)

end
extraFreature = [0 0];
numExtra = 1;
for i = 1:numLidar
    temp1 = (extraFreature-repmat(landmarkLidar(i,:),numExtra,1)).^2;
    dist1 = sum(temp1,2); % We may calculate Mahalanobis distance later
    %dist1(i) = 1e10;
    indexLL = find(dist1 < minL);
    if indexLL
        extraFreature(indexLL,:) = mean([extraFreature(indexLL,:);landmarkLidar(i,:)]);
    else
        extraFreature = [extraFreature;landmarkLidar(i,:)];
    end
    numExtra = size(extraFreature,1);
end
newMap(1,:) = [];
extraFreature(1,:) = [];
newMap = [newMap;extraFreature];
if plotFig
    plot(landmarkCamera(:,1),landmarkCamera(:,2),'r*');
    hold on;
    legend('Camera map')
    k = waitforbuttonpress;
    plot(originalLidar(:,1),originalLidar(:,2),'go');
    hold on;
    legend('Camera map','Lidar map');
    k = waitforbuttonpress;
    plot(newMap(:,1),newMap(:,2),'kp');
    hold on;
    legend('Camera map','Lidar map','Fused Map');
    k = waitforbuttonpress;
    plot(landmarkGT(:,1),landmarkGT(:,2),'bs');
    legend('Camera map','Lidar map','Fused Map','Ground Truth');
end
end
