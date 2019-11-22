filebase = 'C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\KLH\Side View\Localized Aligning\Pre-processed\Particle Images'; 
startFrame = 1; 
endFrame = 26;
% Initialize slice below. If you know the sizes M and N in advance, you can
% skip the following two lines
filename=[filebase, num2str(1,'%2d'),'.png'] ;

% load Ribosome_Map.mat;
cd(Particles_images);
D1 = dir('*.png');

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
    temp = (imread(D1(i).name));
    myslice(:,:,i) = (temp);
%     imagesc(myslice(:,:,i)); 
%     colormap('gray') 
%     drawnow 
        
end
% try isosurface with isovalue = 0.5
figure
isosurface(myslice,0.5)
% % Uncomment to try the sliceomatic function
% figure;
% sliceomatic(myslice)

temp=double(imread(filename)); 
[M,N] = size(temp(:,:,1));
myslice = zeros(M,N,endFrame-startFrame+1);
for i=startFrame:endFrame
    filename=[filebase, num2str(i,'%2d'),'.tif'] ;
    temp=double(imresize(imread(filename), 0.5)); 
    % convert to gray scale
    myslice(:,:,i) = rgb2gray(temp);
    imagesc(myslice(:,:,i)); 
    colormap('gray') 
    drawnow 
end
% try isosurface with isovalue = 0.5
figure
isosurface(myslice,0.5)
% % Uncomment to try the sliceomatic function
% figure;
% sliceomatic(myslice)