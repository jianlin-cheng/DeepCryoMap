close all;
clear all;

Binary_mask_images='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Images\Apoferritin\Binary Masks';
Particles_images='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Images\Apoferritin\Preprocessed Particle Images';
Particles_images1='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Images\Apoferritin\Orginal Particle Images';

fixed  = rgb2gray(imread('Refrences.png'));
% figure;imshow(fixed);

% load Ribosome_Map.mat;
cd(Binary_mask_images);
D1 = dir('*.png');
% map=cell(numel(D1),2);
p = 250;
q = 250;
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
        One_particle_mask_image = im2double(bwareafilt(Logical_particle_mask_image,1));
%         figure;imshow(One_particle_mask_image,[]);
        
        %% Read the orginal Particles...
        cd(Particles_images);
        D = dir('*.png');
        particle_image = im2double(rgb2gray((imread(D(i).name))));
        cd(Particles_images1);
        D = dir('*.png');
        orginal_particle_image = im2double(rgb2gray((imread(D(i).name))));
%         figure;imshow(particle_image,[]);
        figure; imshowpair(orginal_particle_image,particle_image,'montage'); title('Original and Preprocessed Particle Image ')

        img=particle_image;
        [row,col,color_mp]= size(img);
        % Convert the image from 2D to 1D image space...
        img_vector = img(:);
        % specify number of clusters
        Clusters=4;    
        Cluster = cell(1,Clusters);
        Cluster(:) = {zeros(size(img_vector,1),1);};
        % Range       
        range = max(img_vector) - min(img_vector);
        % Determine the # of steps
        stepv = range/Clusters;
        % Cluster initialization
        K=stepv:stepv:max(img_vector);
        for ii=1:size(img_vector,1)
            difference = abs(K-img_vector(ii));
            [y,ind]=min(difference);
            Cluster{ind}(ii)=img_vector(ii);
        end
        cluster_1=reshape(Cluster{1,1},[row col]);
        C1=cluster_1;C1(cluster_1~=0)=1;
        cluster_2=reshape(Cluster{1,2},[row col]);
        C2=cluster_2;C2(cluster_2~=0)=2;
        cluster_3=reshape(Cluster{1,3},[row col]);
        C3=cluster_3;C3(cluster_3~=0)=3;
        cluster_4=reshape(Cluster{1,4},[row col]);
        C4=cluster_4;C4(cluster_4~=0)=4;
        %
        cluster1=bwareaopen(cluster_1,5000);
        figure; imshow(cluster_1,[]);title('Clustered (Binary Mask) Particle Image');
        binIM=cluster1;
        SE=strel('disk',1);
        cluster2=imerode(binIM,SE);
        k=imfill(cluster2,'holes');
        BW = imclose(k,SE);
        cluster3=bwareaopen(BW,0);
        cluster3=imdilate(cluster3,SE);
        figure; imshow(cluster3);title('CryoEM-Binary Mask Image');
        
        % Connect the Circle Edges....
        % Convex hull of outer ring points
        BWin = imclearborder(cluster3);
        CH1 = bwconvhull(cluster3);
        % Remove outer ring points
        se = strel('disk',10);
        CH1b = imerode(CH1,se);
        BW1 = BWin & CH1b;
        % Convex hull of inner ring points
        CH2 = bwconvhull(BW1);
        % Extract the area between inner and outer convex hull
        BWout = CH1 - CH2;
%         figure; imshow(BWout);
        % fill the image 
        k=imfill(BWout,'holes');
%         figure; imshow(k);

        % Smoothed The binary Mask
        windowSize = 55;
        kernel = ones(windowSize) / windowSize ^ 2;
        blurryImage = conv2(single(k), kernel, 'same');
        binaryImage = blurryImage > 0.5; % Rethreshold
%         figure; imshow(binaryImage);title('Particles Picking Results');
        
        % Circle Detection using Modified CHT
        [centers, radii, metric] = imfindcircles(binaryImage,[25 200]);
        [m,n]=size(centers);
        if m==1
    %         figure; imshow(binaryImage);title('Particles Picking Results');
    %         hold on;
    %         viscircles(centers, radii,'EdgeColor','b');
    %         plot(centers(:,1), centers(:,2), 'r+')
    %         hold off;
            IMG_final=double(0*cluster3);

            % Draw a Perfect Circle 
            RGB = insertShape(IMG_final,'FilledCircle',[centers(1,1) centers(1,2) radii],'LineWidth',2); 
            final_mask=double(imbinarize(RGB(:,:,2)));
    %         figure;imshow(final_mask,[]);

            % Localized Particle Image  
            BW3 = imdilate(final_mask, strel('disk',15));
            BW2 = im2double(BW3);
%             figure;imshow(BW2,[]);
            localized_particle_image=orginal_particle_image.*BW2;
            figure;imshow(localized_particle_image,[]);

%             figure;imshowpair(particle_image,localized_particle_image,'montage');  title('Whole Particle Vs. Localized Particle Image ')

            % Unify the Particle Images 
            particle_Mask_unit_8 = im2double(mat2gray((BW2)));
            [m n]=size(BW2);
            K=particle_Mask_unit_8;
            K_pad = padarray(K, [floor((p-m)/2) floor((q-n)/2)],'post');
            Unify_particle_Mask = padarray(K_pad, [ceil((p-m)/2) ceil((q-n)/2)],'pre');
            
            localized_particle_image_unit_8=im2double(mat2gray(localized_particle_image));
            [m n]=size(localized_particle_image_unit_8);
            K=localized_particle_image_unit_8;
            K_pad = padarray(K, [floor((p-m)/2) floor((q-n)/2)],'post');
            Unify_Particle_Image = padarray(K_pad, [ceil((p-m)/2) ceil((q-n)/2)],'pre');
%             figure;imshow(particle_2);

            %% Save the alignement results 
              % save the localized particle images 
              cd('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\Top View\Localized Approach\Orginal');
              imwrite(Unify_Particle_Image, ['Particle_Sample_' num2str(i) '.png'],'png'); 
%               cd('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\Top View\Localized Approach\Preprocessed\Binary Mask');
%               imwrite(Unify_particle_Mask, ['Particle_Mask_' num2str(i) '.png'],'png'); 
%               % save the regular particle iomages 
%               cd('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\Top View\Regular Approach\Preprocessed\Particle Image');
%               imwrite(particle_Mask_unit_8, ['Particle_Sample_' num2str(i) '.png'],'png'); 
    %           cd('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\Localized Aligning\Original\Particle Mask');
    %           imwrite(BW2, ['Particle_Mask_' num2str(i) '.png'],'png'); 
    end

cd(Binary_mask_images);
   end



