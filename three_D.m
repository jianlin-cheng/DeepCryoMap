%% Create dataSet of random particle interations in 3d space
for i=1:5000
    if i == 1
        dataSet = [rand()*100 rand()*100 rand()*100 rand() i];
    else
        dataSet(i,:) = [rand()*100 rand()*100 rand()*100 rand() i];
    end
end
% dataSet = [x y z interactionStrength imageNumber]

xLimits = [min(dataSet(:,1)) max(dataSet(:,1))];
yLimits = [min(dataSet(:,2)) max(dataSet(:,2))];
zLimits = [min(dataSet(:,3)) max(dataSet(:,3))];

binSize = 10; % Number of bins to split each spatial dimention into
binXInterval = (xLimits(2)-xLimits(1))/binSize;
binYInterval = (yLimits(2)-yLimits(1))/binSize;
binZInterval = (zLimits(2)-zLimits(1))/binSize;

histo = [];
for i=xLimits(1)+(binSize/2):binXInterval:xLimits(2) + (binSize/2)
    for j=yLimits(1)+(binSize/2):binYInterval:yLimits(2) + (binSize/2)
        for k=zLimits(1)+(binSize/2):binZInterval:zLimits(2) + (binSize/2)
            %% Filter out particle interactions found within the current spatial bin
            idx = find((dataSet(:,1) > (i - binSize)) .* (dataSet(:,1) < i));
            temp = dataSet(idx,:);
            idx = find((temp(:,2) > (j - binSize)) .* (temp(:,2) < j));
            temp = temp(idx,:);
            idx = find((temp(:,3) > (k - binSize)) .* (temp(:,3) < k));
            temp = temp(idx,:);
            %% Add up all interaction strengths found within this bin
            histo = [histo; i j k sum(temp(:,4))];
        end
    end
end
%% Remove bins with no particle interactions
idx = find(histo(:,4)>0);
histo = histo(idx,:);
numberOfImages = max(dataSet(:,5));
%% Plot result
PointSizeMultiplier = 100000;
scatter3(histo(:,1).*binXInterval + xLimits(1),histo(:,2).*binYInterval + yLimits(1),histo(:,3).*binZInterval + zLimits(1),(histo(:,4)/numberOfImages)*PointSizeMultiplier,(histo(:,4)/numberOfImages));
colormap hot;
%Size and color represent the average interaction intensity over time
