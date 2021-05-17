%% Stabilize frames

%Camera is moving at known rate avg_x_disp (set to 0 if camera isn't moving)
avg_x_disp =  -0.1; %pix/time step

%Size of images (y,x)
imsize = [700,1200];
[X,Y]=meshgrid([1:imsize(2)],[1:imsize(1)]);

%Frame rate
fps = 120;
deltaT = 1/fps;

%Define region with a constant feature(s) in absolute frame of reference
map_ROI = [201,1000, 551,670];
crop_ROI = [2, 100]; % y and x additional cropping to "fit" second frame within first and avoid boundary problems

%If processing multiple videos (directories of frames), name sets and specify which to 'for' over
namesets = {'sample1','sample2'};
%  ^should be names of data sets that have images in directory read_dir (below)

%Create directory to save displacement corrections, directory to save transformed images
save_dir = 'Tformcorrection/';
write_dir = 'Tformimages/';
mkdir(save_dir)

%Read directory
read_dir = 'Readimages/';

%If recorded, import timestamps for each frame in each set
time_dir = 'Timeimages/';
importtime = true;

for index = 1:numel(namesets)
    
    name_dir = namesets{index};
    
    fileID = fopen([save_dir name_dir '_Tform.txt'],'w');
    
    mkdir([write_dir name_dir '/'])
    
    direc_im = dir([read_dir name_dir '*.bmp']); %adjust file name and type as necessary
    
    if importtime
        time_file = dir([time_dir name_dir '*.txt']);
        time_im = importdata([time_dir name_dir time_file.name]);
    else
        time_im = deltaT * [0:numel(direc_im)-1]; %Assumes stable frame rate and no frame dropping
    end
    
    startim = true;
    for imindex = 1:numel(direc_im)
        
        time_im2 = time_im(imindex);
        im2 = imread([direc_im(1).folder '/' direc_im(imindex).name]);
        
        %******************************************************
        %******************************************************
        %***Any pre-processing of images should happen here,***
        %***such as correction for rotations                ***
        %******************************************************
        %******************************************************
        
        if startim
            %No translation corrections for first image
            yx = [0,0];
            time_im1 = time_im2;
            startim = false;
        else
            %Crop both images in a region that with unchanging features in absolute frame
            %But crop previous frame with larger area in which the new image should fit given fps
            im1_crop = single( im1(map_ROI(3)-crop_ROI(1):map_ROI(4)+crop_ROI(1),...
                map_ROI(1)-crop_ROI(2):map_ROI(2)+crop_ROI(2)) ) ;
            im2_crop = single( im2(map_ROI(3):map_ROI(4),map_ROI(1):map_ROI(2)) ) ;
            
            %%*** Translation through normalized correlation
            %Casting as single makes more efficient computation here
            transcorr = normxcorr2( im_rot_crop, im_prev_crop );
            %Find peak of correlation
            cen = find(transcorr == max(transcorr(:)));
            %Pick neighborhood around peak
            [ceny,cenx] = ind2sub(size(transcorr),cen);
            BB = [cenx-3,ceny-3,7,7]; %7x7 neighborhood for fit
            com_reg = transcorr( ceil(BB(2)):ceil(BB(2))+BB(4)-1, ceil(BB(1)):ceil(BB(1))+BB(3)-1 );
            M = sum(com_reg(:));
            C = NaN(1,2);
            ind = zeros([size(com_reg),2]);
            %Centroid of region
            for dims = 1:2
                shp = ones(1,3);
                shp(dims) = size(com_reg,dims);
                rep = size(com_reg);
                rep(dims) = 1;
                ind(:,:,dims) = repmat(reshape(1:size(com_reg,dims),shp),rep);
                C(dims) = sum(sum(ind(:,:,dims).*com_reg))./M;
            end
            ind = cat(3, ind(:,:,2), ind(:,:,1));
            
            %Estimate widths of 2D gaussian near xcorr peak
            BB2 = [cenx-10,ceny-10,21,21]; %9x9 neighborhood
            regs = regionprops( imbinarize(transcorr(ceil(BB2(2)):ceil(BB2(2))+BB2(4)-1, ...
                ceil(BB2(1)):ceil(BB2(1))+BB2(3)-1), max(transcorr(:))/sqrt(2.71828)), ...
                'BoundingBox','Area'); %estimate 1/e height for est width
            if numel(regs)>1
                allarea = NaN(numel(regs),1);
                for regar = 1:numel(regs)
                    allarea(regar) = regs(regar).Area;
                end
                [~, theregind] = max(allarea);
                regs = regs(theregind);
            end
            width_est_y = regs.BoundingBox(4)/2;
            width_est_x = regs.BoundingBox(3)/2;
            
            %2D Gaussian fit with best guesses from above (adjustment of parameters necessary)
            g = @(A,X) A(1)*exp( -((X(:,:,1)-A(2)).^2/(2*A(3)^2) + (X(:,:,2)-A(4)).^2/(2*A(5)^2)) );
            A0 = [max(com_reg(:)),C(2),width_est_x,C(1),width_est_y];   % Inital (guess) parameters
            lb = [max(com_reg(:))-0.05, C(2)-1, width_est_x-1.5, C(1)-1, width_est_y-1.5];
            ub = [max(com_reg(:))+0.05, C(2)+1, width_est_x+1.5, C(1)+1, width_est_y+1.5];
            A = lsqcurvefit(g,A0,ind,com_reg,lb,ub,optimset('Display','off')); %[A,resnorm,res,flag,output] =
            
            yx = [A(4)+ceil(BB(2))-1, A(2)+ceil(BB(1))-1]- size(transcorr)/2 - [0.5,0.5] ;
            %[-0.5,-0.5] may not be necessary. Check if transcorr has even or odd number of elements
            
        end
        
        %Compute imposed average displacement needed to keep up with camera motion
        xdisp = avg_x_disp * (time_im2-time_im1)/deltaT;
        
        %Translate image im2 to match with image im1, then advance by xdisp (which is 0 if cams stationary)
        im2 = imtranslate( im2,flip(yx)-[xdisp,0],'linear'); %advance by fixed amt
        
        %Save image
        imwrite(im2, [write_dir name_dir '/' direc_im(imindex).name]);
        
        if writedata
            fprintf(fileID, ['%.7f %.7f %.7f %.7f\n'], time_im2, yx(2)-xdisp, yx(1), (time_im2-time_im1)/deltaT);
        end
        
        %Store image and time for next frame to be correlated with
        im1 = im2;
        time_im1 = time_im2;
        
    end
    
end

if writedata
    fclose(fileID);
end

%clean up
clearvars
close all
clc

