% close all;
% clear all;

% Binary_mask_images='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Particles Images\Apoferritin\Binary Masks';
Particles_images='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\Ribosome\Irregular Shape\Preprocessed\Particle Image';
Particles_images1='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\Ribosome\Irregular Shape\Orginal\Particle Image';

Code_dire='C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement';


fixed  = rgb2gray(imread('Refrences.png'));
% figure;imshow(fixed);

% load Ribosome_Map.mat;
cd(Particles_images);
D1 = dir('*.png');
% map=cell(numel(D1),2);
p = 500;
q = 500;
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
%         Mask_image = (imread(D1(i).name));
%         Logical_particle_mask_image=logical(Mask_image);
%         One_particle_mask_image = im2double(bwareafilt(Logical_particle_mask_image,1));
%         figure;imshow(One_particle_mask_image,[]);
        
        %% Read the orginal Particles...
        cd(Particles_images);
        D = dir('*.png');
        particle_image = im2double(((imread(D(i).name))));
        preprocessed_particle_image=particle_image;
%         figure;imshow(particle_image);
        cd(Particles_images1);
        D = dir('*.png');
        orginal_particle_image = im2double(rgb2gray((imread(D(i).name))));
        figure; imshowpair(orginal_particle_image,particle_image,'montage'); title('Original and Preprocessed Particle Image ')
        %%
        %%
        A=preprocessed_particle_image;
        [L,N] = superpixels(A,100);
        figure
        BW = boundarymask(L);
        imshow(imoverlay(A,BW,'cyan'));
        
        outputImage = zeros(size(A),'like',A);
        idx = label2idx(L);
        numRows = size(A,1);
        numCols = size(A,2);
        for labelVal = 1:N
            redIdx = idx{labelVal};
            outputImage(redIdx) = mean(A(redIdx));
        end    

        figure;imshow(outputImage)
        %
        figure; imshowpair(particle_image,outputImage,'montage'); title('Original and Preprocessed Particle Image ')
        %%
        %% Particle Clustering...
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
        cluster1=bwareaopen(cluster_1,1000);
        figure;imshow(cluster1,[]);
% % 
%         out = bwskel(cluster1);
%         figure; imshow(labeloverlay(orginal_particle_image,out,'Transparency',0))

% I = imresize(cluster1,[512 512]);
% S = qtdecomp(I,.9);
% blocks = repmat(uint8(0),size(S));
% 
% for dim = [512 256 128 64 32 16 8 4 2 1];    
%   numblocks = length(find(S==dim));    
%   if (numblocks > 0)        
%     values = repmat(uint8(1),[dim dim numblocks]);
%     values(2:dim,2:dim,:) = 0;
%     blocks = qtsetblk(blocks,S,dim,values);
%   end
% end
% blocks(end,1:end) = 1;
% blocks(1:end,end) = 1;
% imshow(I), figure, imshow(blocks,[])
                
                
          
        
%         
%         binIM=cluster1;
%         SE=strel('disk',1);
%         cluster2=imerode(binIM,SE);
%         k=imfill(cluster2,'holes');
%         BW = imclose(k,SE);
%         cluster3=bwareaopen(BW,0);
%         cluster3=imdilate(cluster3,SE);
%         figure; imshowpair(cluster_1,cluster3,'montage'); title('Orginal Clustered (Binary Mask) Vs. preprocessed Particle Image ')
% 
%         
%         
% cd(Code_dire);
%                
% figure;imshow(cluster3);
% [y,x] = find(bwperim(cluster3));
% h = convhull(x,y,'Simplify',true);
% x_hull = x(h);
% y_hull = y(h);
% delete(hull_line);
% 
% figure; imshow(cluster3)
% hold on
% plot(x_hull,y_hull,'r-*','LineWidth',2,'MarkerSize',12)
% hold off
% title('A Blob''s Convex Hull and Its Vertices')
        
        %% Post Processing...
        % Connect the Circle Edges....
        % Convex hull of outer ring points
%         BWin = imclearborder(cluster3);
%         CH1 = bwconvhull(cluster3);
%         figure;imshow(CH1,[]);title('Smoothed Convex hull');
%         
%         % Remove outer ring points
%         se = strel('disk',10);
%         CH1b = imerode(CH1,se);
%         BW1 = BWin & CH1b;
% 
%         % Convex hull of inner ring points
%         CH2 = bwconvhull(BW1);
%         windowSize = 55;
%         kernel = ones(windowSize) / windowSize ^ 2;
%         blurryImage = conv2(single(CH2), kernel, 'same');
%         CH2 = blurryImage > 0.5; % Rethreshold
%         figure;imshow(CH2,[]);title('Convex hull of inner ring points');
% 
%         % Extract the area between inner and outer convex hull
%         BWout = CH1 - CH2;
%         figure;imshow(BWout,[]);title('Area (inner, outer) convex hull');
% 
%         % fill the image 
%         k=imfill(BWout,'holes');
%         figure;imshow(k,[]);title('Full Circular Object');
% 
%         % Smoothed The binary Mask
%         windowSize = 55;
%         kernel = ones(windowSize) / windowSize ^ 2;
%         blurryImage = conv2(single(k), kernel, 'same');
%         binaryImage = blurryImage > 0.5; % Rethreshold
%         figure; imshow(binaryImage);title('Particles Picking Results');
%         
% %%
%         % Circle Detection using Modified CHT
%         [centers, radii, metric] = imfindcircles(binaryImage,[25 200]);
%         [m,n]=size(centers);
%         if m==1
% %             figure; imshow(binaryImage);title('Modified CHT');
% %             hold on;
% %             viscircles(centers, radii,'EdgeColor','b');
% %             plot(centers(:,1), centers(:,2), 'r+')
% %             hold off;
%             IMG_final=double(0*cluster3);
% 
%             % Draw a Perfect Circle 
%             RGB = insertShape(IMG_final,'FilledCircle',[centers(1,1) centers(1,2) radii],'LineWidth',2); 
%             final_mask=double(imbinarize(RGB(:,:,2)));
% %             figure;imshow(final_mask,[]);title('Perfect Circular Picking Mask');
% 
%             % Localized Particle Image  
%             BW3 = imdilate(final_mask, strel('disk',16));
%             BW2 = im2double(BW3);
% %             figure;imshow(BW2,[]);title('Crentral Aligned Mask');
% 
%             Orginal_localized_particle_image=orginal_particle_image.*BW2;
%             Preprocessed_localized_particle_image=preprocessed_particle_image.*BW2;
% %             figure;imshow(Preprocessed_localized_particle_image,[]);title('Localized Top-View Particle');
% %             figure;imshowpair(particle_image,localized_particle_image,'montage');  title('Whole Particle Vs. Localized Particle Image ')
% 
%             % Unify the Particle Mask 
%             [m n]=size(BW2);
%             K=BW2;
%             K_pad = padarray(K, [floor((p-m)/2) floor((q-n)/2)],'post');
%             Unify_particle_Mask = padarray(K_pad, [ceil((p-m)/2) ceil((q-n)/2)],'pre');
% %             figure;imshow(Unify_particle_Mask);
% 
%             % Unify the Orginal Particle Image 
%             [m n]=size(Orginal_localized_particle_image);
%             K=Orginal_localized_particle_image;
%             K_pad = padarray(K, [floor((p-m)/2) floor((q-n)/2)],'post');
%             Unify_Orginal_Particle_Image = padarray(K_pad, [ceil((p-m)/2) ceil((q-n)/2)],'pre');
% %             figure;imshow(Unify_Orginal_Particle_Image);
% 
%             % Unify the Preprocessed Particle Image 
%             [m n]=size(Preprocessed_localized_particle_image);
%             K=Preprocessed_localized_particle_image;
%             K_pad = padarray(K, [floor((p-m)/2) floor((q-n)/2)],'post');
%             Unify_Preprocessed_Particle_Image = padarray(K_pad, [ceil((p-m)/2) ceil((q-n)/2)],'pre');
% %             figure;imshow(Unify_Preprocessed_Particle_Image);
% 
%             % Alighn the Center
%             stat = regionprops(Unify_particle_Mask,'centroid');
%             x1=round(stat.Centroid(1));
%             x2=round(stat.Centroid(2));
%             Unify_Orginal_Particle_Image = imcrop(Unify_Orginal_Particle_Image,[x1-150 x2-150 300 300]);
%             Unify_particle_Mask = imcrop(Unify_particle_Mask,[x1-150 x2-150 300 300]);
%             Unify_Preprocessed_Particle_Image = imcrop(Unify_Preprocessed_Particle_Image,[x1-150 x2-150 300 300]);
%             %% Save the alignement results 
%               % save the localized orginal particle images 
%               cd('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\KLH\Top View\Localized Approach\Orginal\Particle Image');
%               imwrite(Unify_Orginal_Particle_Image, ['Particle_Sample_' num2str(i) '.png'],'png'); 
%               cd('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\KLH\Top View\Localized Approach\Orginal\Binary Mask');
%               imwrite(Unify_particle_Mask, ['Particle_Mask_' num2str(i) '.png'],'png'); 
%                    
%               % save the regular particle iomages 
%               cd('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\KLH\Top View\Localized Approach\Preprocessed\Particle Image');
%               imwrite(Unify_Preprocessed_Particle_Image, ['Particle_Sample_' num2str(i) '.png'],'png'); 
% %               cd('C:\Users\Adil Al-Azzawi\Desktop\Particle Alignement\Alignement Results\Apoferritin\Top View\Localized Approach\Preprocessed\Binary Mask');
% %               imwrite(BW2, ['Particle_Mask_' num2str(i) '.png'],'png'); 
%         end
cd(Particles_images);
pause;
   end



