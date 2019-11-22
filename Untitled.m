close all;
clear all;

cd ('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Masks\Side_View');

% load Ribosome_Map.mat;
cd(cryo_EM_images1);
D1 = dir('*.png');
% map=cell(numel(D1),2);
% Change the direction to the image folder
   for i = 1:numel(D1)
        close all;
        % start to read cancer image
        clc; 
        disp('======================================================================================');        
        disp('  P  A  R  T  I  C  L  E -  D  E  T  E  C  T  I  O  N  - A N D -  P  I  C  K  I  N  G ');
        disp('  ------------- T R A I N I N G - D A T A S E T - P R E P E R A T I O N ------------- ');
        disp('======================================================================================');
        disp(' ');
        fprintf('Cryo-EM Image No. : %d\n',i');  
        
        % read the MRI image from the spesific directory
        Mask_image = (imread(D1(i).name));



fixed  = rgb2gray(imread('ref.png'));
moving   = rgb2gray(imread('ref5.png'));

figure; imshowpair(moving,fixed,'montage'); title('Unregistered')
figure; imshowpair(moving,fixed); title('Unregistered')

[optimizer,metric] = imregconfig('multimodal');

%%
movingRegisteredDefault = imregister(moving,fixed,'affine',optimizer,metric);
% figure; imshowpair(movingRegisteredDefault,fixed); title('A: Default Registration')
%%
optimizer.InitialRadius = optimizer.InitialRadius/3.5;
movingRegisteredAdjustedInitialRadius = imregister(moving,fixed,'affine',optimizer,metric);
% figure; imshowpair(movingRegisteredAdjustedInitialRadius,fixed); title('B: Adjusted InitialRadius')

%%
optimizer.MaximumIterations = 500;
movingRegisteredAdjustedInitialRadius300 = imregister(moving,fixed,'affine',optimizer,metric);
figure; imshowpair(movingRegisteredAdjustedInitialRadius300,fixed); title('C: Adjusted InitialRadius, MaximumIterations = 300')
%%
tformSimilarity = imregtform(moving,fixed,'similarity',optimizer,metric);
Rfixed = imref2d(size(fixed));
movingRegisteredRigid = imwarp(moving,tformSimilarity,'OutputView',Rfixed);
% figure; imshowpair(movingRegisteredRigid, fixed); title('D: Registration Based on Similarity Transformation Model')
%%
tformSimilarity.T
movingRegisteredAffineWithIC = imregister(moving,fixed,'affine',optimizer,metric,...
    'InitialTransformation',tformSimilarity);
% figure; imshowpair(movingRegisteredAffineWithIC,fixed); title('E: Registration from Affine Model Based on Similarity Initial Condition')
%%

%%
movingRegistered=movingRegisteredAdjustedInitialRadius300;

% figure;imshowpair(fixed, movingRegistered,'Scaling','joint')


labeledImage1 = bwlabel(moving);
measurements1 = regionprops(labeledImage1, 'Orientation', 'MajorAxisLength', 'Centroid');
Angle1 =round( [measurements1.Orientation])

labeledImage2 = bwlabel(movingRegistered);
measurements2 = regionprops(labeledImage2, 'Orientation', 'MajorAxisLength', 'Centroid');
Angle2 =round( [measurements2.Orientation])


Angle=round(abs(Angle1-90))

J = imrotate(moving,round(Angle));
% figure;imshow(moving);
% figure; imshow(J);

figure;imshowpair(fixed,J,'montage')

