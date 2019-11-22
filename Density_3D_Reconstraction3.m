close all;
clear all;

% Binary_mask_images='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Images\Apoferritin\Binary Masks';
Particles_images='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\KLH\Side View\Localized Aligning\Original\Particle Images';
% Particles_images1='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Images\Original\Top_view';

cd(Particles_images);
D1 = dir('*.png');
D1=100;
X = zeros(350,350,D1);

   for i = 2:numel(D1)
       tic;
        close all;
        % start to read cancer image
        clc; 
        disp('======================================================================================');        
        disp('  P  A  R  T  I  C  L  E -  D  E  T  E  C  T  I  O  N  - A N D -  P  I  C  K  I  N  G ');
        disp('  ------------- T R A I N I N G - D A T A S E T - P R E P E R A T I O N ------------- ');
        disp('======================================================================================');
        disp(' ');
        fprintf('Cryo-EM Image No. : %d\n',i');  
        
        %% Read the Particles Masks
        particle_image = (imread(D1(i).name));
        X(:,:,i) = particle_image;     

   end

   montage(permute(X,[1 2 4 3]));

   
map = hsv(90);
XR =X;
Ds = smooth3(XR);
hiso = patch(isosurface(Ds,5),'FaceColor','blue','EdgeColor','none');
hcap = patch(isocaps(XR,5),'FaceColor','interp','EdgeColor','none');
colormap(map)
daspect(gca,[1,1,.4])
lightangle(305,30);
fig = gcf;
fig.Renderer = 'zbuffer';
lighting phong
isonormals(Ds,hiso)
hcap.AmbientStrength = .6;
hiso.SpecularColorReflectance = 0;
hiso.SpecularExponent = 50;
ax = gca;
ax.View = [215,30];
ax.Box = 'On';
axis tight
title('Original Data');



focalLength    = [350, 350]; 
principalPoint = [350, 350];
imageSize      = [350, 350];

intrinsics = cameraIntrinsics(focalLength,principalPoint,imageSize);


% load Ribosome_Map.mat;
cd(Particles_images);
D1 = dir('*.png');
fixed  = ((imread(D1(1).name)));
[m,n]=size(fixed);
I1=zeros(m,n,3);
I1(:,:,1)=im2double(fixed);
I1(:,:,2)=im2double(fixed);
I1(:,:,3)=im2double(fixed);
figure;imshow(I1);

% Detect features. Increasing 'NumOctaves' helps detect large-scale
% features in high-resolution images. Use an ROI to eliminate spurious
% features around the edges of the image.
border = 50;
roi = [border, border, size(fixed, 2)- 2*border, size(fixed, 1)- 2*border];
prevPoints   = detectSURFFeatures(fixed, 'NumOctaves', 8, 'ROI', roi);

% Extract features. Using 'Upright' features improves matching, as long as
% the camera motion involves little or no in-plane rotation.
prevFeatures = extractFeatures(fixed, prevPoints, 'Upright', true);

% Create an empty viewSet object to manage the data associated with each
% view.
vSet = viewSet;

% Add the first view. Place the camera associated with the first view
% and the origin, oriented along the Z-axis.
viewId = 1;
vSet = addView(vSet, viewId, 'Points', prevPoints, 'Orientation', ...
    eye(3, 'like', prevPoints.Location), 'Location', ...
    zeros(1, 3, 'like', prevPoints.Location));

% Detect feature points
imagePoints1 = detectMinEigenFeatures(rgb2gray(I2), 'MinQuality', 0.1);

% Visualize detected points
figure
imshow(I1, 'InitialMagnification', 50);
title('150 Strongest Corners from the First Image');
hold on
plot(selectStrongest(imagePoints1, 150));

% Create the point tracker
tracker = vision.PointTracker('MaxBidirectionalError', 1, 'NumPyramidLevels', 5);

% Initialize the point tracker
imagePoints1 = imagePoints1.Location;
initialize(tracker, imagePoints1, I1);         

% Load precomputed camera parameters
load upToScaleReconstructionCameraParameters.mat

% Change the direction to the image folder
   for i = 2:numel(D1)
       tic;
        close all;
        % start to read cancer image
        clc; 
        disp('======================================================================================');        
        disp('  P  A  R  T  I  C  L  E -  D  E  T  E  C  T  I  O  N  - A N D -  P  I  C  K  I  N  G ');
        disp('  ------------- T R A I N I N G - D A T A S E T - P R E P E R A T I O N ------------- ');
        disp('======================================================================================');
        disp(' ');
        fprintf('Cryo-EM Image No. : %d\n',i');  
        
        %% Read the Particles Masks
        particle_image = (imread(D1(i).name));
        I2(:,:,1)=im2double(fixed);
        I2(:,:,2)=im2double(fixed);
        I2(:,:,3)=im2double(fixed);
%         figure; imshowpair(I1,I2,'montage'); title('Original and Preprocessed Particle Image ')
        
        % Track the points
        [imagePoints2, validIdx] = step(tracker, I2);
        matchedPoints1 = imagePoints1(validIdx, :);
        matchedPoints2 = imagePoints2(validIdx, :);

        % Visualize correspondences
        figure
        showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2);
        title('Tracked Features');
         
        % Estimate the fundamental matrix
        [E, epipolarInliers] = estimateEssentialMatrix(...
        matchedPoints1, matchedPoints2, cameraParams, 'Confidence', 99.99);
    
        % Estimate the camera pose of current view relative to the previous view.
        % The pose is computed up to scale, meaning that the distance between
        % the cameras in the previous view and the current view is set to 1.
        % This will be corrected by the bundle adjustment.
%         [relativeOrient, relativeLoc, inlierIdx] = helperEstimateRelativePose(...
%             matchedPoints1, matchedPoints2, intrinsics);
        F = estimateFundamentalMatrix(matchedPoints1,matchedPoints2);
        [relativeOrient, relativeLoc, inlierIdx] =  relativeCameraPose(F,intrinsics, matchedPoints1, matchedPoints2);
        
        % Add the current view to the view set.
        vSet = addView(vSet, i, 'Points', currPoints);
        
        % Store the point matches between the previous and the current views.
        vSet = addConnection(vSet, i-1, i, 'Matches', indexPairs(inlierIdx,:));
        
        
        
        
        
        %% Find Point Correspondences Between The Images
        % Detect feature points
        imagePoints1 = detectMinEigenFeatures(rgb2gray(I1), 'MinQuality', 0.1);

        % Visualize detected points
        figure
        imshow(I1, 'InitialMagnification', 50);
        title('150 Strongest Corners from the First Image');
        hold on
        plot(selectStrongest(imagePoints1, 150));

        % Create the point tracker
        tracker = vision.PointTracker('MaxBidirectionalError', 1, 'NumPyramidLevels', 5);

        % Initialize the point tracker
        imagePoints1 = imagePoints1.Location;
        initialize(tracker, imagePoints1, I1);

        % Track the points
        [imagePoints2, validIdx] = step(tracker, I2);
        matchedPoints1 = imagePoints1(validIdx, :);
        matchedPoints2 = imagePoints2(validIdx, :);

        % Estimate the Essential Matrix
        load upToScaleReconstructionCameraParameters.mat
        % Estimate the fundamental matrix
        [E, epipolarInliers] = estimateEssentialMatrix(...
            matchedPoints1, matchedPoints2, cameraParams, 'Confidence', 99.99);

        % Find epipolar inliers
        inlierPoints1 = matchedPoints1(epipolarInliers, :);
        inlierPoints2 = matchedPoints2(epipolarInliers, :);

        % Display inlier matches
        figure
        showMatchedFeatures(I1, I2, inlierPoints1, inlierPoints2);
        title('Epipolar Inliers');
        
        % Visualize correspondences
        figure
        showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2);
        title('Tracked Features');
        
        % Compute the Camera Pose
        [orient, loc] = relativeCameraPose(E, cameraParams, inlierPoints1, inlierPoints2);

        %% Reconstruct the 3-D Locations of Matched Points
        % Detect dense feature points. Use an ROI to exclude points close to the
        % image edges.
        roi = [30, 30, size(I1, 2) - 30, size(I1, 1) - 30];
        imagePoints1 = detectMinEigenFeatures(rgb2gray(I1), 'ROI', roi, ...
            'MinQuality', 0.001);

        % Create the point tracker
        tracker = vision.PointTracker('MaxBidirectionalError', 1, 'NumPyramidLevels', 5);

        % Initialize the point tracker
        imagePoints1 = imagePoints1.Location;
        initialize(tracker, imagePoints1, I1);

        % Track the points
        [imagePoints2, validIdx] = step(tracker, I2);
        matchedPoints1 = imagePoints1(validIdx, :);
        matchedPoints2 = imagePoints2(validIdx, :);

        % Estimate the camera pose of current view relative to the previous view.
        % The pose is computed up to scale, meaning that the distance between
        % the cameras in the previous view and the current view is set to 1.
        % This will be corrected by the bundle adjustment.
        [relativeOrient, relativeLoc, inlierIdx] = helperEstimateRelativePose(...
        matchedPoints1, matchedPoints2, intrinsics);
        
        
        % Compute the camera matrices for each position of the camera
        % The first camera is at the origin looking along the Z-axis. Thus, its
        % rotation matrix is identity, and its translation vector is 0.
        camMatrix1 = cameraMatrix(cameraParams, eye(3), [0 0 0]);

        % Compute extrinsics of the second camera
        [R, t] = cameraPoseToExtrinsics(orient, loc);
        camMatrix2 = cameraMatrix(cameraParams, R, t);

        % Compute the 3-D points
        points3D = triangulate(matchedPoints1, matchedPoints2, camMatrix1, camMatrix2);

        % Get the color of each reconstructed point
        numPixels = size(I1, 1) * size(I1, 2);
        allColors = reshape(I1, [numPixels, 3]);
        colorIdx = sub2ind([size(I1, 1), size(I1, 2)], round(matchedPoints1(:,2)), ...
            round(matchedPoints1(:, 1)));
        color = allColors(colorIdx, :);

        % Create the point cloud
        ptCloud = pointCloud(points3D, 'Color', color);
        
%         %% Display the 3-D Point Cloud
%         % Visualize the camera locations and orientations
%         cameraSize = 0.3;
%         figure
%         
%         plotCamera('Size', cameraSize, 'Color', 'r', 'Label', '1', 'Opacity', 0);
%         hold on
%         grid on
%         plotCamera('Location', loc, 'Orientation', orient, 'Size', cameraSize, ...
%             'Color', 'b', 'Label', '2', 'Opacity', 0);
% 
%         % Visualize the point cloud
%         pcshow(ptCloud, 'VerticalAxis', 'y', 'VerticalAxisDir', 'down', ...
%             'MarkerSize', 45);
% 
%         % Rotate and zoom the plot
%         camorbit(0, -30);
%         camzoom(1.0);
% 
%         % Label the axes
%         xlabel('x-axis');
%         ylabel('y-axis');
%         zlabel('z-axis')
%         title('Up to Scale Reconstruction of the Scene');
        
   end
