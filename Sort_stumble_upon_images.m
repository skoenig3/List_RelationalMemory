%% Automatically Process images from Flickr Pics
% modified from SortListSQImages.m by Seth Konig 3/4/15
% mofified to sort though 10,000's of stumbleupon images simulatenously.
% Algorithms the same, just sorting into subfoldrs and removing images that
% are exact duplicates of others. Code does steps 1-3 simultaneously to
% improve processing timesby importing images once instead of 3x
% since I have sooo many images.
% [0] Remove duplicate images
% [1] Remove images that are too small or are grayscale
% [2] Resize images
% [3] Select for "interesting" Pictures
%% [0] Remove dubplicated images in each directory
% check for duplicate images because mutlipe urls for the same image could
% be present on the same webpage. Each folder can contain images pulled
% from multiple webpages, but each webpage should be contained to a single
% folder.

main_directory = 'C:\Users\seth.koenig\Documents\SU images 3\';%where all the images are stored
files = dir(main_directory);  % assume starting from current directory
filenames = {files.name};
subdirs = filenames([files.isdir]);

subdirs(1:2) = [];%these are '.' and '..';

for sb = 1:length(subdirs);
    image_list = ls([main_directory subdirs{sb},'\*.bmp']); %get all the .bmp files in this directory
    total_images = size(image_list,1); %total number of bmp files in directory
    last_image = zeros(600,800,3); %start with a blank image so it' can't match anything
    imgnum = 1;
    while imgnum < total_images
        next_img = imread([main_directory subdirs{sb} '\' image_list(imgnum,:)]);
        if all(size(last_image) == size(next_img)); %if they are the same size could be the same image
            if all(all(all(next_img == last_image)));%if this image is the same as the last image
                delete([main_directory subdirs{sb} '\' image_list(imgnum,:)]); %delete this file
                imgnum = imgnum+1;
            else%if it's not a duplicate
                last_image = next_img;
                imgnum = imgnum+1;
            end
        else %if not the same size probably not a duplicate
            last_image = next_img;
            imgnum = imgnum+1;
        end
    end
end


%% [1-3] Going to try to do all steps at once
% step 1 checking for images that are too small or grayscalle shouldn't be an
% issues since they shouldn't have been pulled with the original code
% going to resize and then check for entropy frequency stuff simultaneously
% within 1 for loop

imageX = 800;% desired horizontal image size
imageY = 600;% desired vertical image size

% sobel filters detect edges!
sobelx = [1     2   1;
    0     0   0;
    -1    -2  -1;];

sobely = [1     2   1;
    0     0   0;
    -1    -2  -1;];

main_directory = 'C:\Users\seth.koenig\Documents\SU images 3\';%where all the images are stored

% smalldir and graydir shouldn't overflow with too many images so there can
% only be 1 directory for each. The good and the bad may overflow so going
% to create directories for every 5000 images in each type
smalldir = [main_directory 'Images that are too small\'];
mkdir(smalldir);
graydir = [main_directory 'Gray scale images\'];
mkdir(graydir);

gooddir_base = 'Good';
baddir_base = 'Bad';

imgs_in_gooddir = 5000;%no images starting off in good directory but need to trigger 1st if statement
imgs_in_baddir = 5000; %no images starting off in bad directory but need to trigger 2bd if statement
last_dir = ones(1,2); %number of directories for good and bad directories

files = dir(main_directory);  % assume starting from current directory
filenames = {files.name};
subdirs = filenames([files.isdir]);

for sb = 1:length(subdirs);
    if imgs_in_gooddir >= 5000
        gooddir = [main_directory gooddir_base num2str(last_dir(1)) '\'];
        mkdir(gooddir)
        last_dir(1) = last_dir(1)+1;
        imgs_in_gooddir = 0;
    end
    if imgs_in_baddir >= 5000
        baddir = [main_directory baddir_base num2str(last_dir(2)) '\'];
        mkdir(baddir)
        last_dir(2) = last_dir(2)+1;
        imgs_in_baddir = 0;
    end
    image_list = ls([main_directory subdirs{sb},'\*.bmp']); %get all the .bmp files in this directory
    
    for im = 1:size(image_list,1)%total number of bmp files in directory
        imgname = [main_directory subdirs{sb} '\' image_list(im,:)];%name and path to image
        img = imread(imgname);%load image
        
        %check if image is grayscale, if it is move to graydir
        if size(img,3) == 1
            movefile(imgname,graydir);
            continue
        elseif all(all(img(:,:,1) == img(:,:,2))) ...
                || all(all(img(:,:,2) == img(:,:,3))) || all(all(img(:,:,1) == img(:,:,3)))
            movefile(imgname,graydir);
            continue
        end
        
        %check if image is too small
        if size(img,1) < imageY || size(img,2) < imageX
            movefile(imgname,smalldir);
            continue
        end
        
        %resize the image, changing size may alter entropy/frequency
        %content so resize first
        if ~ all(size(img) == [imageY,imageX,3])%if the image isn't already the desired size
            %then resize it, otherwize don't waste computer power
            img = imresize(img,[imageY,imageX]);
        end
        
        % then  select for intereseting pictures/removed bad ones
        % use variability in pixel intensities (image entropy) and percent edginess
        % (high freqency content) to weed out images lacking in dynamic content.
        % This is an objective process but isn't totally perfect. Essentially code
        % removes pictures with a lot of background and low range of image
        % intensities (i.e. images that are too bright or too dark)
        
        img = rgb2gray(img);
        entropyvalues = entropy(img);%pixel intesnity entropy
        xedges = imfilter(img,sobelx);
        yedges = imfilter(img,sobely);
        edgevalues = mean2(xedges+yedges); %edgineess
        
        if entropyvalues > 7 && (edgevalues > 15) &&  (edgevalues <85) %if good
            movefile(imgname,gooddir);
            imgs_in_gooddir = imgs_in_gooddir +1; %1 more image added to the gooddir
        else %bad
            movefile(imgname,baddir);
            imgs_in_badddir = imgs_in_baddir +1; %1 more image added to the baddir
        end
    end
end
%% Resize all images
clear, clc

imageX = 800;% desired horizontal image size
imageY = 600;% desired vertical image size

imgdir = 'C:\Users\seth.koenig\Documents\SU images 3\';
base_dir = 'Good';

for g = 1:3;
    current_dir = [imgdir base_dir num2str(g) '\'];
    image_list = ls([current_dir,'\*.bmp']); %get all the .bmp files in this directory
    for im = 1:size(image_list,1)%total number of bmp files in directory
        imgname = [current_dir '\' image_list(im,:)];%name and path to image
        img = imread(imgname);%load image
        if any(size(img) ~= [imageY,imageX,3]);%if image isn't 600x800 already, then resize it
            img = imresize(img,[imageY,imageX]);
            imwrite(img,imgname,'bmp');
        end
    end
end