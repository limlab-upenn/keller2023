clear all
% close all
tic
%% load images and set parameters 
%export2 = 'C:\Users\blim\Dropbox\eve_endogenous_bursting\replicate02\stripe_1-3\';  %mac
% pathname = '/Users/bomyilim/Library/CloudStorage/GoogleDrive-bomyilim@seas.upenn.edu/My Drive/Students/Gaochen/Dl_modeling_Brian/Gaochen/Dlven_snaPP7/020823/';  %windows
pathname = '/Users/bomyilim/Library/CloudStorage/GoogleDrive-bomyilim@seas.upenn.edu/Shared drives/Lim_Lab/Sam/Elongation Project/sogdist_sna_MLP/sdsnamlpn2/';
% [filename, pathname]= uigetfile({'*.mat','Data files(*.mat)'},'Choose');
% load([pathname filename]);
% load([export2 'trajectories_ms2.mat']);
load([pathname 'segmentation_lineage14.mat']);


start = 1;
T = start*0.516 + 0.516*(1:size(lineage_cx,1));


%% First, get the threshold for each time point 
for i=1:size(nuc_lineage,2) % nuclei index
    for j=1:size(nuc_lineage,1)
        
%         nuc_lineage(100,281) = 299; 
        c_image = convex_image{j,nuc_lineage(j,i)};
        c_image = imresize(c_image,1.05);
        
        % x, y distance translation
        r1 = regionprops(c_image,'PixelList','centroid');  % get pixellist and centroids for the nucleus
        
        % find the coordinates in the entire image (rather than within a
        % single nucleus)
        dx = lineage_cx(j,i) - r1.Centroid(1); dx = round(dx);
        dy = lineage_cy(j,i) - r1.Centroid(2); dy = round(dy);
        
        r1.PixelList(:,1) = r1.PixelList(:,1) + dx;
        r1.PixelList(:,2) = r1.PixelList(:,2) + dy;
        
        % correction for the coordinates (smaller/bigger than 1, 512)
        ddum = find(r1.PixelList(:,1) < 1); % too small
        r1.PixelList(ddum,1) = 1;
        ddum = find(r1.PixelList(:,2) < 1);
        r1.PixelList(ddum,2) = 1;
        ddum = find(r1.PixelList(:,1) > size(images_pp7,2));
        r1.PixelList(ddum,1) = size(images_pp7,2);
        ddum = find(r1.PixelList(:,2) > size(images_pp7,1));
        r1.PixelList(ddum,2) = size(images_pp7,1);
        
        % x and y coordinates of a given nucleus in the entire frame
        xx = r1.PixelList(:,1);
        yy = r1.PixelList(:,2);
        
        % remove the background, by subtracting the 50 percentile value
        I_pp7 = images_pp7(:,:,j) - prctile(reshape(images_pp7(:,:,j),[1 size(images_ms2,1)*size(images_ms2,2)]),45);
        I_ms2 = images_ms2(:,:,j) - prctile(reshape(images_ms2(:,:,j),[1 size(images_ms2,1)*size(images_ms2,2)]),45);
%         I_ms2 = images_ms2(:,:,j);
        
        PP = zeros(size(images_pp7,1),size(images_pp7,2),'uint16');
        O = zeros(size(images_ms2,1),size(images_ms2,2),'uint16');
        O1 = zeros(size(images_ms2,1),size(images_ms2,2),'uint16');
        P1 = zeros(size(images_pp7,1),size(images_pp7,2),'uint16');
        A = zeros(size(images_pp7,1),size(images_pp7,2),'uint16');
        
        % extract pp7 and ms2 signals from all the pixel indeces
        for k=1:size(xx,1)
            Ns(k) = images_nuc(yy(k),xx(k),j);
            Ps(k) = images_pp7(yy(k),xx(k),j);
            Ms(k) = images_ms2(yy(k),xx(k),j); % ms2 signal
%             Ps(k) = I_pp7(yy(k),xx(k));
%             Ms(k) = I_ms2(yy(k),xx(k));
            O(yy(k),xx(k)) = I_ms2(yy(k),xx(k));
            O1(yy(k),xx(k)) = images_ms2(yy(k),xx(k),j);
            PP(yy(k),xx(k)) = I_pp7(yy(k),xx(k));
            P1(yy(k),xx(k)) = images_pp7(yy(k),xx(k),j);
            A(yy(k),xx(k),3) = 20000;
            
        end
        
%         Ps = double(Ps); Ms = double(Ms);
        
%         P_ratio(j,i) = max(Ps)/median(Ps); M_ratio(j,i) = max(Ms)/median(Ms);
        
        O = imcrop(O,[min(xx) min(yy) max(xx)-min(xx)+1 max(yy)-min(yy)+1]);
        O = imresize(O,2);
        O1 = imcrop(O1,[min(xx) min(yy) max(xx)-min(xx)+1 max(yy)-min(yy)+1]);
        O1 = imresize(O1,2);
        PP = imcrop(PP,[min(xx) min(yy) max(xx)-min(xx)+1 max(yy)-min(yy)+1]);
        PP = imresize(PP,2);
        P1 = imcrop(P1,[min(xx) min(yy) max(xx)-min(xx)+1 max(yy)-min(yy)+1]);
        P1 = imresize(P1,2);
        
        O_bw = im2bw(O,graythresh(O));
        rpo = regionprops(O_bw,'centroid');
        s = size(rpo,1);
        
        if s>2 || graythresh(O)<0.01% too many objects
            O_bw = im2bw(O,0.09);
            rpo = regionprops(O_bw,'centroid');
            if size(rpo,1)>2
                O_bw = im2bw(O,0.1);
                rpo = regionprops(O_bw,'centroid');
                s = size(rpo,1);
            else
                s = size(rpo,1);
            end
        end
            
        
        A(:,:,2) = I_ms2*5;
        A(:,:,1) = images_pp7(:,:,j)*5;
        
        
        Mm(j,i) = median(Ms);
        Mmax(j,i) = max(Ms);
%         M(j,i) = max(Ms);
        Pm(j,i) = mean(Ps);
        Pmax(j,i) = max(Ps);
        Nm(j,i) = median(Ns);
        Nmax(j,i) = max(Ns);
        
        Msort = sort(Ms,'descend');
        M(j,i) = mean(Msort(1:2));
        Psort = sort(Ps,'descend');
        P(j,i) = mean(Psort(1:3));
        Nsort = sort(Ns,'descend');
        N(j,i) = mean(Nsort(1:2));

        
        
%             save(sprintf('%strajectories.mat',pathname),'nuc_lineage','lineage_cx','lineage_cy',...
%         'M','T','Mm','Pm','P','Mmax','Pmax','Nm','Nmax','N');



%         figure(10); plot(T,M(:,i),'bo-'); hold on;
%         figure(11); plot(T,P(:,i),'ro-'); hold on;
%         
        clear xx yy Ms Ps Ns
        
%         w = waitforbuttonpress;
%         if w==1
%             continue
%         end
%         clf(1); clf(2); clf(3);
    end
    
end
save(sprintf('%strajectories.mat',pathname),'nuc_lineage','lineage_cx','lineage_cy',...
        'M','T','Mm','Pm','P','Mmax','Pmax','Nm','Nmax','N');         
toc            
%             
            
        
        
       
