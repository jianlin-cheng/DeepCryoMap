close all;
clear all;

Binary_mask_images='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Masks\Side_View';
Particles_images='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Images\Original\Side_View';

fixed  = rgb2gray(imread('Refrences.png'));
% figure;imshow(fixed);

% load Ribosome_Map.mat;
cd(Binary_mask_images);
D1 = dir('*.png');
% map=cell(numel(D1),2);
% Change the direction to the image folder
similarity=zeros(numel(D1),1); 
consuming_time=zeros(numel(D1),1); 
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
        Mask_image = (imread(D1(i).name));
        Logical_particle_mask_image=logical(Mask_image);
        One_particle_mask_image = bwareafilt(Logical_particle_mask_image,1);

        moving=double(One_particle_mask_image);
        
        %% Read the orginal Particles...
        cd(Particles_images);
        D = dir('*.png');
        particle_image = uint16(rgb2gray(imread(D(i).name)));
%         figure; imshowpair(particle_image,Mask_image,'montage'); title('Original Particle Image and Masks')
%         figure; imshowpair(moving,fixed,'montage'); title('Unregistered Particle Masks and the Refrence')
%         figure; imshowpair(moving,fixed); title('Unregistered Masks Projection')

        %% Particles Alignement...
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
%         figure; imshowpair(movingRegisteredAdjustedInitialRadius300,fixed); title('C: Adjusted InitialRadius, MaximumIterations = 300')
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

        %% Extract the Aligned Angles
        movingRegistered=movingRegisteredAdjustedInitialRadius300;
%         figure;imshowpair(fixed, movingRegistered,'Scaling','joint'); title('Registered Particle Masks and the Refrence')
%         figure; imshowpair(moving,movingRegistered,'montage'); title('Aligned Mask Vs. Orginal Mask')

        labeledImage1 = bwlabel(moving);
        measurements1 = regionprops(labeledImage1, 'Orientation', 'MajorAxisLength', 'Centroid');
        Angle1 =round( [measurements1.Orientation])

        labeledImage2 = bwlabel(movingRegistered);
        measurements2 = regionprops(labeledImage2, 'Orientation', 'MajorAxisLength', 'Centroid');
        Angle2 =round( [measurements2.Orientation]);

        Angle=round(abs(Angle1-90));

        % Align the orignal particle
        J1 = imrotate(particle_image,round(Angle));
        consuming_time(i,1)=toc;
        [m n]=size(particle_image);
        K=particle_image;
        p = 350;
        q = 350;
        K_pad = padarray(K, [floor((p-m)/2) floor((q-n)/2)],'post');
        particle_image = padarray(K_pad, [ceil((p-m)/2) ceil((q-n)/2)],'pre');
%         figure;imshow(particle_image,[]);
        
        [m n]=size(J1);
        K=J1;
        p = 350;
        q = 350;
        K_pad = padarray(K, [floor((p-m)/2) floor((q-n)/2)],'post');
        J1 = padarray(K_pad, [ceil((p-m)/2) ceil((q-n)/2)],'pre');
%         figure;imshowpair(particle_image,J1,'montage');  title('Orginal Particle Vs. Aligned Particle')

         % Find the similarity betwwen the orginal the allignement particle
        [ssimval,ssimmap] = ssim(J1,particle_image);
%         figure; imshow(ssimmap,[])
%         title(['Local SSIM Map with Global SSIM Value: ',num2str(ssimval)])
        similarity(i,10)=ssimval;
        
        % Localized the particle Image 
        particle_image = uint16(rgb2gray(imread(D(i).name)));
        moving = uint16(moving);
        
        % Enlarge particle binary mask
        BW3 = imdilate(moving, strel('disk',12));
        BW2 = uint16(BW3);
        particle_image1=particle_image.*BW2;
        J1 = uint8(imrotate(particle_image,round(Angle)));
        J2 = uint8(imrotate(particle_image1,round(Angle)));
        
        figure;imshowpair(J1,J2,'montage');  title('Aligned Particle Vs. Localized Aligned Particle')

        %% Save the alignement results 
          cd('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\Regular Aligning\Original\Particle Images');
          imwrite(J1, ['Particle_Sample_' num2str(i) '.png'],'png'); 
          cd('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\Regular Aligning\Original\Particle Mask');
          imwrite(movingRegistered, ['Particle_Mask_' num2str(i) '.png'],'png'); 
          
cd(Binary_mask_images);
   end



