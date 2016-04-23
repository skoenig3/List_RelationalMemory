function [itmlist,clrchng_locations,first_img_item,imgs] = read_ListRM_itm_and_cnd_files(itmfile)
% writteen by Seth Konig May, 2015. Modified from read_ListSQ_itm_and_cnd_files.mat
% Function imports item file and and grabs condition file to determine which items
% are associated with which condition (itmlist) since conditions are randomly
% organized. 

ITMFile = ['C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Item and CND files\' itmfile];
CNDFile = ['C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Item and CND files\' itmfile(1:end-4) '.cnd'];

itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

first_img_item = [];
for i = 6:size(itmfil,1); %first 5 are hearder, background, etc
    str = textscan(itmfil(i,:),'%s');
    if ~isempty(strfind(str{1}{end},'.bmp')) && isempty(first_img_item)
        first_img_item  = str2num((str{1}{1}));
        break
    end
end

itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

%probably a more efficient way of doing this since there's only 25 color
%change locations, but it may be easier for processing later 
% also get which images were displayed. Should be the same across all items
% and condition sets but what the hey. 
imgs = zeros(2,96);
clrchng_locations = cell(1,length(itmlist));
for cnd = 1:length(itmlist)
    if itmlist(cnd) < first_img_item %then it is a clrchng trial
        str = textscan(itmfil(itmlist(cnd)+6,:),'%d');
        clrchng_locations{cnd} =  double(str{1}(4:5)); %location in dva
    else
        str = textscan(itmfil(itmlist(cnd)+6,:),'%s');
        imgnum = str2double(str{1}{end}(14:15)); 
        if imgs(1,imgnum) == 0; %novel image
            imgs(1,imgnum) = cnd;
        else
            imgs(2,imgnum) = cnd;
        end
    end
end