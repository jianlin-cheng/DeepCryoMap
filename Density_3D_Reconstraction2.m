close all;
clear all;

% Binary_mask_images='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Images\Apoferritin\Binary Masks';
Particles_images='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\KLH\Side View\Localized Aligning\Original\Particle Images';
% Particles_images1='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Images\Original\Top_view';

% load Ribosome_Map.mat;
cd(Particles_images);
D1 = dir('*.png');
fixed  = ((imread(D1(1).name)));
% [m,n]=size(fixed);
% I1=zeros(m,n,3);
% I1(:,:,1)=im2double(fixed);
% I1(:,:,2)=im2double(fixed);
% I1(:,:,3)=im2double(fixed);
% figure;imshow(I1);

I=fixed;
% Detect features. Increasing 'NumOctaves' helps detect large-scale
% features in high-resolution images. Use an ROI to eliminate spurious
% features around the edges of the image.
border = 50;
roi = [border, border, size(I, 2)- 2*border, size(I, 1)- 2*border];
prevPoints   = detectSURFFeatures(I, 'NumOctaves', 8, 'ROI', roi);

% Extract features. Using 'Upright' features improves matching, as long as
% the camera motion involves little or no in-plane rotation.
prevFeatures = extractFeatures(I, prevPoints, 'Upright', true);

% Create an empty viewSet object to manage the data associated with each
% view.
vSet = viewSet;

% Add the first view. Place the camera associated with the first view
% and the origin, oriented along the Z-axis.
viewId = 1;
vSet = addView(vSet, viewId, 'Points', prevPoints, 'Orientation', ...
    eye(3, 'like', prevPoints.Location), 'Location', ...
    zeros(1, 3, 'like', prevPoints.Location));

% Change the direction to the image folder
   for i = 1:numel(D1)
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
        I=particle_image;
        
%         I2(:,:,1)=im2double(fixed);
%         I2(:,:,2)=im2double(fixed);
%         I2(:,:,3)=im2double(fixed);
%         figure; imshowpair(I1,I2,'montage'); title('Original and Preprocessed Particle Image ');
 
        % Detect, extract and match features.
        currPoints   = detectSURFFeatures(I, 'NumOctaves', 8, 'ROI', roi);
        currFeatures = extractFeatures(I, currPoints, 'Upright', true);    
        indexPairs = matchFeatures(prevFeatures, currFeatures, ...
            'MaxRatio', .7, 'Unique',  true);
        
        % Select matched points.
        matchedPoints1 = prevPoints(indexPairs(:, 1));
        matchedPoints2 = currPoints(indexPairs(:, 2));
        
        % Estimate the camera pose of current view relative to the previous view.
        % The pose is computed up to scale, meaning that the distance between
        % the cameras in the previous view and the current view is set to 1.
        % This will be corrected by the bundle adjustment.
        [relativeOrient, relativeLoc, inlierIdx] = helperEstimateRelativePose(...
        matchedPoints1, matchedPoints2, intrinsics);
   end
