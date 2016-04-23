%code written to easily double check if there are duplicate images in a
%folder or across folder. Code briefly calculates mean and standar
%deviation of image intensity as well as entropy to calculate 3 simple
%values per image. 

%for images within 1 directory
current_dir = 'C:\Users\seth.koenig\Documents\Use these\';
image_list = ls([current_dir,'\*.bmp']); %get all the .bmp files in this directory
values = NaN(size(image_list,1),3);
for im = 1:size(image_list,1)%total number of bmp files in directory
    imgname = [current_dir '\' image_list(im,:)];%name and path to image
    img = imread(imgname);%load image
    img = rgb2gray(img);
    entropyvalues = entropy(img);%pixel intesnity entropy
    values(im,1) = entropyvalues;
    values(im,2) = mean2(img);%mean pixel intensity
    values(im,3) = std2(img);%standard deviation of pixel intensity
end

repeat_dir = 'Repeats\';%move repeated images here
x = pdist2(values,values);
x(find(tril(x))) = NaN;
x(find(eye(im))) = NaN;
[i,j] = find(x == 0 ); 
%%
while ~isempty(j)
    repeatimg = image_list(j(1),:);
    movefile([current_dir repeatimg],[current_dir repeat_dir]);
    all_js = find(j == j(1) | i == j(1)); %should include j(1)
    i(all_js) = [];
    j(all_js) = [];
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for images in 2 directories 
current_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Sequence Task_RelationalMemory\Used\';
image_list = ls([current_dir,'\*.bmp']); %get all the .bmp files in this directory
values = NaN(size(image_list,1),3);
for im = 1:size(image_list,1)%total number of bmp files in directory
    imgname = [current_dir '\' image_list(im,:)];%name and path to image
    img = imread(imgname);%load image
    img = rgb2gray(img);
    entropyvalues = entropy(img);%pixel intesnity entropy
    values(im,1) = entropyvalues;
    values(im,2) = mean2(img);%mean pixel intensity
    values(im,3) = std2(img);%standard deviation of pixel intensity
end

%2nd directory should be the directory you are removeing images from 
current_dir = 'C:\Users\seth.koenig\Documents\Use these\';
repeat_dir = 'Repeats\';%move repeated images here
image_list = ls([current_dir,'\*.bmp']); %get all the .bmp files in this directory
values2 = NaN(size(image_list,1),3);
for im = 1:size(image_list,1)%total number of bmp files in directory
    imgname = [current_dir '\' image_list(im,:)];%name and path to image
    img = imread(imgname);%load image
    img = rgb2gray(img);
    entropyvalues = entropy(img);%pixel intesnity entropy
    values2(im,1) = entropyvalues;
    values2(im,2) = mean2(img);%mean pixel intensity
    values2(im,3) = std2(img);%standard deviation of pixel intensity
end
%%
x = pdist2(values,values2);
[i,j] = find(x == 0); 
for jj =1:length(j)
    repeatimg = image_list(j(jj),:);
    movefile([current_dir repeatimg],[current_dir repeat_dir]);
end
%% Double Check For ones that are really close since they may be compressed/decrompress verions of each other
[i,j] = find(x < 0.2); %should only need 0.1 but may want to doule check higher values just in case. 0.1 estimated empircally. 
for ii = 1:length(i);
    img1 = imread([current_dir image_list(j(ii),:)]);
    
    img2 = imread([current_dir image_list(i(ii),:)]);
    figure
    subplot(1,2,1)
    imshow(img1) 
    title(image_list(j(ii),:))
    box off 
    axis off
    subplot(1,2,2)
    imshow(img2)
        title(image_list(i(ii),:))
    box off
    axis off
    pause
    close 
end
