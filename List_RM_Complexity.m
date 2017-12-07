%check image complexity

% sobel filters detect edges!
sobelx = [1     2   1;
    0     0   0;
    -1    -2  -1;];

sobely = [1     2   1;
    0     0   0;
    -1    -2  -1;];

entropyvalues = NaN(1,3840);
salience_entropyvalues = NaN(1,3840);
edgevalues =NaN(1,3840);
img_count = 1;

for set = 1:40
    disp(['Set #' num2str(set)])
    if set < 10
        img_dir = ['C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Image Sets\LRM0' num2str(set) '\'];
    else
        img_dir = ['C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Image Sets\LRM' num2str(set) '\'];
    end
      for img = 1:96
        if set < 10
            if img < 10
                imgr = imread([img_dir 'S0' num2str(set) 'I0' num2str(img) '.bmp']);
                load([img_dir 'S0' num2str(set) 'I0' num2str(img) '-saliencemap.mat'],'fullmap')
            else
                imgr = imread([img_dir 'S0' num2str(set) 'I' num2str(img) '.bmp']);
                load([img_dir 'S0' num2str(set) 'I' num2str(img) '-saliencemap.mat'],'fullmap')
            end
        else
            if img < 10
                imgr = imread([img_dir 'S' num2str(set) 'I0' num2str(img) '.bmp']);
                load([img_dir 'S' num2str(set) 'I0' num2str(img) '-saliencemap.mat'],'fullmap')
            else
                imgr = imread([img_dir 'S' num2str(set) 'I' num2str(img) '.bmp']);
                load([img_dir 'S' num2str(set) 'I' num2str(img) '-saliencemap.mat'],'fullmap')
            end 
        end
        

        entropyvalues(img_count) = entropy(imgr);%pixel intesnity entropy
        salience_entropyvalues(img_count) = entropy(fullmap);
        xedges = imfilter(imgr,sobelx);
        yedges = imfilter(imgr,sobely);
        edgevalues(img_count) = mean2(xedges+yedges); %edgineess
        
        img_count = img_count+1;
    end
end

%%