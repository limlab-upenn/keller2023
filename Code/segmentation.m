clear all
close all

%% load images and set parameters
% 
% [filename, pathname]= uigetfile({'*.tif','Image Files(*.tif)'},'Choose');
pathname = '/Users/bomyilim/Library/CloudStorage/GoogleDrive-bomyilim@seas.upenn.edu/My Drive/Students/Gaochen/Dl_modeling_Brian/Gaochen/Dlven_snaPP7/020823/';
filename = 'MAX_020823-Dl-ven-sna-PP7-1.tif';

fname = [pathname filename];
info = imfinfo(fname);
num_images = numel(info);

export = [pathname 'segmentation/'];
export2 = pathname;


zero = 1; % first frame of nuclear cycle 14
start = 12; % total number of frames in the movie
finish = 72; % frame index where nuclei signal's not accumulated

his = 2; % Histone channel (nuclei)
ms2 = 1; % ms2 channel
pp7 = 3; % pp7 channel

saveseg = 0; % 0 = don't save, 1 = save
codepath = 2; % 1=take original, 2=take corrected seg
numPass = 1; % 1 = first pass/third pass, 2 = second pass

maskON = 0; % 1= mask, 0=no mask
borderON = 0;
if borderON == 1
        load([pathname 'borders.mat']);
end
% The number of components that are allowed to be loss upon resizing and
% watershedding 
% First pass has 2, allowing for a few fusions
% The second pass, as it usually occurs after some hand correction, does
% not allow the loss of any components
if numPass == 1
%     lossAllowed = 2;
    lossAllowed = 5;
else
%     lossAllowed = 0;
    lossAllowed = 0;
end

T = 0.13*(1:finish-start+1)+5; %should be accurate 
numComp = zeros(finish-start+1, 3);

if maskON == 1
    mask = imread([pathname 'mask.tif']);
    h = fspecial('gaussian',3,3);
    mask = imfilter(mask,h);
    mask1 = imadjust(mask);
    bw_m = im2bw(mask1,graythresh(mask1));
    bw_m = bwareaopen(bw_m,50);
    bw_m = imfill(bw_m,'holes');
end
%% Nuclei segmentation
% First step is to process the nuclei channel to enhance the image quality
% for successful segmentation of nuclei
tic
for n=start:finish
    %If numPass = 2, we only want to look at every third frame (this is to
    %speed things up). The second pass must be run after propagating
    %changes to other frames after hand-correction. 
    if numPass == 2
        if mod(((n+1)-start)+2, 5)  ~= 0
            continue
        end
    end
    disp(n);
    I = imread(fname, 3*(n-1)+his); % bfp
    M = imread(fname, 3*(n-1)+ms2); % ms2
    P = imread(fname, 3*(n-1)+pp7); % pp7
    % Define the frame size
    if codepath == 1
        if maskON == 1
            for ii=1:size(I,1)
                for jj=1:size(I,2)
                    if bw_m(ii,jj)==0
                        I(ii,jj)=0;
                    end
                end
            end
        end
        
        h = fspecial('gaussian',2,2); % gaussian filter
        P1 = imfilter(I,h);
        
        se = strel('disk',25);
        P2 = imtophat(P1,se); % top-hat filtering (correct for uneven illumination)
        P2 = imadjust(P2); % further contrast enhancement
        
        % change into a binary image (if gthresh is too big,use gthresh-0.05)
        bw = im2bw(P2, graythresh(P2)); %im2bw(P2,graythresh(P2)-thresh);
        bw = bwareaopen(bw,10);
        
        %early watershedding
        D = -bwdist(~bw);
        mask = imextendedmin(D,2);
        %imshowpair(bw,mask,'blend')
        D2 = imimposemin(D,mask);
        Ld2 = watershed(D2);
        bw5 = bw;
        bw5(Ld2 == 0) = 0;
        if numPass == 2
            F = imread(sprintf('%s%s%03d.tif',export,'',n-start+1));
            F = F(:,:,1);
            bw5 = im2bw(F,graythresh(F));
        end
    elseif codepath == 2
        F = imread(sprintf('%s%s%03d.tif',export,'',n-start+1));
        F = F(:,:,1);
        bw5 = im2bw(F,graythresh(F));
        oldBI = bw5;
    else
        s = 'ERROR!';
        disp(s);
    end
    
%     if maskON == 1 
%         
%         load([pathname 'mask.mat']); 
%         cc_tmp = bwconncomp(bw5);
%         rpo_tmp = regionprops(cc_tmp,'Centroid','PixelList');
%         
%         for k=1:size(rpo_tmp,1)
%             in = inpolygon(x1,y1,rpo_tmp(k).PixelList(:,1),rpo_tmp(k).PixelList(:,2));
%             s(k) = sum(in);
%         end
%         mk = find(s>0);
%         for k=1:length(mk)
%             idx = size(rpo_tmp(mk(k)).PixelList(:,1),1);
%             for kk=1:idx
%                 bw5(rpo_tmp(mk(k)).PixelList(kk,2),rpo_tmp(mk(k)).PixelList(kk,1))=0;
%             end
%         end
%         
%     end
    
    
    if borderON == 1
        cc_tmp = bwconncomp(bw5);
        rpo_tmp = regionprops(cc_tmp,'Centroid','PixelList');
        
        % between top to bottm
        for k=1:size(rpo_tmp,1)
            y11(k) = rpo_tmp(k).Centroid(1)*p1(1)+p1(2);
            y22(k) = rpo_tmp(k).Centroid(1)*p2(1)+p2(2);
            y_tmp(k) = rpo_tmp(k).Centroid(2);
        end

        tmp = find(y22 < y_tmp | y11 > y_tmp); % due to flipped x and y on images
        for k=1:length(tmp)
            idx = size(rpo_tmp(tmp(k)).PixelList(:,1),1);
            for kk=1:idx
                bw5(rpo_tmp(tmp(k)).PixelList(kk,2),rpo_tmp(tmp(k)).PixelList(kk,1))=0;
            end
        end
        clear rpo_tmp y11 y22 y_tmp;
    end
    
    if codepath == 1
        disp('first pass segment');
        %% Define each nucleus (connected components)
        cc = bwconncomp(bw5);
        
        % obtain centroids and pixel lists of each nucleus
        rpo = regionprops(cc, 'Centroid','ConvexImage','PixelList');
        centroids = cat(1,rpo.Centroid);
        
        %% Remove the edge nuclei
        for i=1:size(cc.PixelIdxList,2)
            a = rpo(i,1).PixelList;
            
            % find nuclei at the upper and left boundary
            side = find(a(:,1) < 5 | a(:,2) < 5); %the value of 5 can be changed, but 5 seems pretty good in most cases
            if isempty(side) == 1
                edge(i,1) = 0;
            else
                edge(i,1) = 1; %touching the edge
            end
            
            % find nuclei at the bottom and right boundary
            side1 = find(a(:,1) > size(P,2) - 5 | a(:,2) > size(P,1) - 5); %the value of 5 can be changed, but 5 seems pretty good in most cases
            if isempty(side1) == 0
                edge(i,1) = 1;
            end
        end
        
        edge1 = find(edge == 0); % edge nuclei
        
        %% remove the mask nuclei 
        
        %reindex edge nuclei away
        centroids = centroids(edge1,:);
        rpo = rpo(edge1,:);
        disp('edge nuclei removed');
        
        %%Begin adaptative water-shedding
        justRight = 1;
        resFactor = 1.1;
        ctrl = length(centroids);
        numComp(n-start+1, 1) = ctrl;
        numComp(n-start+1, 2) = ctrl;
        oldBI = bw5;
        while(justRight) %while the image hasn't resized TOO big
            newBI = false(size(I,1), size(I,2));
            for i = 1:length(centroids) % for each nucleus, resize
                c_image = rpo(i).ConvexImage;
                c_image = imresize(c_image, resFactor);
                r1 = regionprops(c_image,'PixelList','centroid');
                dx = centroids(i,1) - r1(1).Centroid(1); dx = round(dx);
                dy = centroids(i,2) - r1(1).Centroid(2); dy = round(dy);
                r1(1).PixelList(:,1) = r1(1).PixelList(:,1) + dx;
                r1(1).PixelList(:,2) = r1(1).PixelList(:,2) + dy;
                % correction for the coordinates (smaller/bigger than the framesize)
                ddum = find(r1(1).PixelList(:,1) < 1);
                r1(1).PixelList(ddum,1) = 1;
                ddum = find(r1(1).PixelList(:,2) < 1);
                r1(1).PixelList(ddum,2) = 1;
                ddum = find(r1(1).PixelList(:,1) > size(I,2)); % size,2 = column (x)
                r1(1).PixelList(ddum,1) = size(I,2);
                ddum = find(r1(1).PixelList(:,2) > size(I,1)); % size,1 = row (y)
                r1(1).PixelList(ddum,2) = size(I,1);
                % x and y coordinates of a give nucleus in the entire frame
                xx = r1(1).PixelList(:,1);
                yy = r1(1).PixelList(:,2);
%                 bound = boundary(xx, yy); 
%                 figure(11), plot(xx(bound), yy(bound))
                for j=1:size(xx,1) % number of pixels per nucleus
                     newBI(yy(j), xx(j)) = 1;
                end
            end
            
            %With newly-resized nuclei, watershed 
            D = -bwdist(~newBI); %distance matrix
            mask = imextendedmin(D,2); %smoothed mask
            D2 = imimposemin(D,mask); %apply smoothed mask
            Ld2 = watershed(D2); %watershed beased on smoothed distance matrix
            bw10 = newBI;
            bw10(Ld2 == 0) = 0;
            [labeledImage, numRegions] = bwlabel(bw10);
            if(numRegions < (ctrl - lossAllowed)) %if we lose too many components, we resized too much
                justRight = 0;
                break
            end
            oldBI = bw10;
            numComp(n-start+1, 2) = numRegions;
            numComp(n-start+1, 3) = resFactor;
            resFactor = resFactor + 0.05;
            if resFactor > 1.4
                justRight = 0;
            break
            end
        end
        disp(resFactor); %this is the factor by which our resizing failed
        cc = bwconncomp(oldBI);
        
        % obtain centroids and pixel lists of each nucleus
        rpo = regionprops(cc, 'Centroid','ConvexImage','PixelList');
        centroids = cat(1,rpo.Centroid);
        
        % Re-define the number of objects (exclude the edge ones)
        for i=1:size(centroids,1)
            coord{n-start+1,1}(i,1:size(rpo(i).PixelList(:,1),1)) = rpo(i).PixelList(:,1);
            coord{n-start+1,2}(i,1:size(rpo(i).PixelList(:,2),1)) = rpo(i).PixelList(:,2);
            convex_image{n-start+1,i} = rpo(i).ConvexImage;
        end
        
        centroids_x(n-start+1,1:size(centroids(:,:),1)) = centroids(:,1);
        centroids_y(n-start+1,1:size(centroids(:,:),1)) = centroids(:,2);
        
        %% Store images
        images_pp7(:,:,n-start+1) = P;
        images_ms2(:,:,n-start+1) = M;
        images_seg(:,:,n-start+1) = oldBI;
        images_nuc(:,:,n-start+1) = I;
        disp('images stored');
        % use to save the segmentation file
        if saveseg == 1 || saveseg==2
            F = zeros(size(images_pp7,1),size(images_pp7,2),3,'uint16');
            F(:,:,1) = images_seg(:,:,n-start+1)*22000;
            F(:,:,3) = images_nuc(:,:,n-start+1);
            F(:,:,2) = images_ms2(:,:,n-start+1)*4;
%             F(:,:,2) = images_pp7(:,:,n-start+1)*1.2 + images_ms2(:,:,n-start+1);
            if saveseg == 1
                imwrite(F,sprintf('%s%s%03d.tif',export,'',n-start+1));
            else
                imwrite(F,sprintf('%s%s%03d.tif',export,'search_',n-start+1));
            end
        end
        clear edge centroids bw5
        toc
    elseif codepath == 2
        %% Define each nucleus (connected components)
        cc = bwconncomp(bw5);
        
        % obtain centroids and pixel lists of each nucleus
        rpo = regionprops(cc, 'Centroid','ConvexImage','PixelList');
        centroids = cat(1,rpo.Centroid);
        % Re-define the number of objects (exclude the edge ones)
        for i=1:size(centroids,1)
            coord{n-start+1,1}(i,1:size(rpo(i).PixelList(:,1),1)) = rpo(i).PixelList(:,1);
            coord{n-start+1,2}(i,1:size(rpo(i).PixelList(:,2),1)) = rpo(i).PixelList(:,2);
            convex_image{n-start+1,i} = rpo(i).ConvexImage;
        end
        
        centroids_x(n-start+1,1:size(centroids(:,:),1)) = centroids(:,1);
        centroids_y(n-start+1,1:size(centroids(:,:),1)) = centroids(:,2);
        
         %% Store images
        images_pp7(:,:,n-start+1) = P;
        images_ms2(:,:,n-start+1) = M;
        images_seg(:,:,n-start+1) = oldBI;
        images_nuc(:,:,n-start+1) = I;
    end
end
    %   show series of images
    seq = [1 round(size(images_seg,3)/4) round(size(images_seg,3)/2) size(images_seg,3)];
    for i=1:4
        figure(1); subplot(1,4,i); imshow(images_seg(:,:,seq(i))); hold on;
    end
if numPass == 2
    return;
end
%% Nuclei tracking 
for i=1:size(centroids_x,1)
    tmp = find(centroids_x(i,:) == 0);
    if isempty(tmp) == 1
        endpt(i) = size(centroids_x,2);
    else
        endpt(i) = tmp(1)-1;
    end
end

% re-save the centroid values
Adj(1,:) = 1:endpt(1);
c_nx(1,:) = centroids_x(1,1:endpt(1));
c_ny(1,:) = centroids_y(1,1:endpt(1));


for i=1:endpt(1) % number of nuclei
    for j=2:size(centroids_x,1)
        % find the nucleus with the minimum distance
        [a,b] = min((c_nx(j-1,i)-centroids_x(j,:)).^2 + (c_ny(j-1,i)-centroids_y(j,:)).^2);
        
        % assign the lineage if the distance is short enough
        if a<200
            Adj(j,i) = b;
            c_nx(j,i) = centroids_x(j,b);
            c_ny(j,i) = centroids_y(j,b);
        else
            Adj(j,i) = 0;
            c_nx(j,i) = c_nx(j-1,i);
            c_ny(j,i) = c_ny(j-1,i);
        end
    end
end

for i=1:size(c_nx,2)
    tmp1 = find(Adj(:,i)==0);
    if size(tmp1,1)<1
        A(i) = 1;
    else
        A(i) = 0;
    end
end
AA = find(A==1);
nuc_lineage = Adj(:,AA); % lineage nuclei index

for i=1:size(nuc_lineage,2)
    tmp = find(nuc_lineage(:,i)==0);
    nuc_lineage(tmp,i)=nuc_lineage(tmp-1,i);
end
% store x,y coordinates of the tracked nuclei
for i=1:size(nuc_lineage,1)
    for j=1:size(nuc_lineage,2)
        lineage_cx(i,j) = centroids_x(i,nuc_lineage(i,j));
        lineage_cy(i,j) = centroids_y(i,nuc_lineage(i,j));
    end
end

% remove overlapping nuclei
if codepath == 2
    for i=1:size(nuc_lineage,2)
        tmp=find(nuc_lineage(size(nuc_lineage,1),:)==nuc_lineage(size(nuc_lineage,1),i));
        if length(tmp)>1
            for j=1:length(tmp)
                dum(j) = sum(diff(lineage_cx(:,tmp(j))).^2+diff(lineage_cy(:,tmp(j))).^2);
            end
            [~,loc]=min(dum);
            if tmp(loc)==i
                a(i)=1;
            else
                a(i)=0;
            end
            clear dum
        else
            a(i)=0;
        end       
    end
end
dum=find(a==1);
AA(:,dum)=[];
nuc_lineage(:,dum)=[]; lineage_cx(:,dum)=[]; lineage_cy(:,dum)=[];

for i=1:size(nuc_lineage,2)   
    aa(1:size(nuc_lineage,1)-1,i) = diff(nuc_lineage(:,i));
    b(i) = max(aa(:,i));
end
dum = find(b>100);
AA(:,dum)=[];
nuc_lineage(:,dum)=[]; lineage_cx(:,dum)=[]; lineage_cy(:,dum)=[];


save(sprintf('%ssegmentation_lineage.mat',export2),'centroids_x','centroids_y','convex_image',...
    'images_pp7','images_ms2','images_seg','images_nuc','nuc_lineage','AA','lineage_cx','lineage_cy','T');
csvwrite(sprintf('%ssegmentation_results.csv', export2), numComp);

figure,imshow(images_seg(:,:,ceil((n-start+1)/2))); hold on;
plot(lineage_cx(ceil((n-start+1)/2),:),lineage_cy(ceil((n-start+1)/2),:),'bo');

