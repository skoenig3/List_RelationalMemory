%% [1] Generate Item files
% Modified from SMT_generate_setnum_CND Files Seth Konig 4.14.15
numimages = 96;
for setnum = 1:40
    if setnum < 10
        set = ['ListRM0' num2str(setnum)];
    else
        set = ['ListRM' num2str(setnum)];
    end
    
    fid = fopen([set '.itm'],'w+');
    
    spacingx = {'-12.0','-6.00',' 0.00',' 6.00',' 12.0'};
    spacingy = {'-8.00','-4.00',' 0.00',' 4.00',' 8.00'};
    ngrid = length(spacingx);
    
    line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------';
    line2 =' -4   14      1    0.00    0.00      0   0.50  0.50  0.00              255 255 255 x';
    line3 =' -3   14      1    0.00    0.00      0   0.50  0.50  0.00               75  75  75 x';
    line4 =' -2    1      1    0.00    0.00      0   0.20  0.20  0.00  0.50         10  10  10 x';
    line5 =' -1    1      0    0.00    0.00      0   1.00  1.00  0.00              200 200 200 x';
    line6 ='  0    1      1    0.00    0.00      0   0.15  0.15  0.00              150 150 150 x';
    
    
    for line = 1:6
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    for i = 1:length(spacingx)
        for ii = 1:length(spacingy)
            num1 = num2str(0+3*ngrid*(i-1)+2*(ii-1)+ii);
            num3 = num2str(0+3*ngrid*(i-1)+2*(ii-1)+ii+1);
            num2 = num2str(0+3*ngrid*(i-1)+2*(ii-1)+ii+2);
            if str2double(num1) < 10
                header1 = ['  ' num1];
            else
                header1 = [' ' num1];
            end
            if str2double(num2) < 10
                header2 = ['  ' num2];
            else
                header2 = [' ' num2];
            end
            if str2double(num3) < 10
                header3 = ['  ' num3];
            else
                header3 = [' ' num3];
            end
            str1 = [header1 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50        150 150 150 x' '\r\n'];
            str3 = [header3 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50        150 150 150 x' '\r\n'];
            str2 = [header2 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50        175 175 130 x' '\r\n'];
            fprintf(fid,str1);
            fprintf(fid,str3);
            fprintf(fid,str2);
        end
    end
    
    if setnum < 10
        setzero = '0';
    else
        setzero = '';
    end
    for i = 1:numimages;
        if i < 10
            imgzero = '0';
        else
            imgzero = '';
        end
        if 75+i < 100
            str = [' ' num2str(75+i) '    8           0.00    0.00      0                                  75  75  75 x   C:\\'...
                'LRM' setzero num2str(setnum)  '\\S' setzero num2str(setnum) 'I' imgzero num2str(i) '.bmp\r\n'];
        else
            str = ['' num2str(75+i) '    8           0.00    0.00      0                                  75  75  75 x   C:\\' ...
                'LRM' setzero num2str(setnum)  '\\S' setzero num2str(setnum) 'I' imgzero num2str(i) '.bmp\r\n'];
        end
        fprintf(fid,str);
    end
    fclose(fid);
end
%%
%%---[3] Sort images into Image Sets---%%
% secition sorts images from unused library into image sets and keeps track
% of the original name of files. Original files are put into used libary.
root_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\';
unused_dir = 'Use these\';
used_dir ='Used\';

cd([root_dir unused_dir]);
num_images_per_set = 96;
for set = 1:40
    d=dir([root_dir unused_dir '*.bmp']);
    if set < 10
        set_dir = [root_dir 'LRM0' num2str(set) '\']; %shortened from ListSQ or directory name too long for cortex
    else
        set_dir = [root_dir 'LRM' num2str(set) '\']; %shortened from ListSQ or directory name too long for cortex
    end
    if exist(set_dir,'dir')
        disp('Image Set already exists') %do not let it rewrite original folders
        continue;
    else
        mkdir(set_dir)
    end
    rr = randperm(length(d));
    for img = 1:num_images_per_set;
        original = d(rr(img)).name;
        if set < 10;
            if img < 10
                new = ['S0' num2str(set) 'I0' num2str(img) '.bmp'];
            else
                new = ['S0' num2str(set) 'I' num2str(img) '.bmp'];
            end
        else
            if img < 10
                new = ['S' num2str(set) 'I0' num2str(img) '.bmp'];
            else
                new = ['S' num2str(set) 'I' num2str(img) '.bmp'];
            end
        end
        copyfile([root_dir unused_dir original],[set_dir new])
        movefile([root_dir unused_dir original],[root_dir used_dir original],'f')
    end
end
%%
%----------------------------------------------------------------%
%write unique random conditions files for each day to psuedorandomize image
%clrchng trials
line1 = 'COND# BCKGND TIMING FIX_ID ---COLOR-PALETTE--- TEST0 TEST1 TEST2 TEST3 TEST4 TEST5 TEST6 TEST7'; %header
line2 = '  1     -3    1      37                         38    39'; %color change line
rand('seed',150415);%do not change seed so that can reproduce each file exactly
for setnum = 1:40
    if setnum < 10
        set = ['ListRM0' num2str(setnum)];
    else
        set = ['ListRM' num2str(setnum)];
    end
    fid = fopen([set '.cnd'],'w+');
    
    for line = 1:2
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    all_clrblocks = [];
    num_clrblocks = 24;
    for blk = 1:num_clrblocks
        all_clrblocks = [all_clrblocks; randperm(25)'];
    end
    
    all_clr_items = [3*(all_clrblocks-1)+1 3*(all_clrblocks-1)+2 3*(all_clrblocks-1)+3];
    
    number_of_images = 96;
    image_spacing = 16;
    number_clrchng_trials_btwn_images = 3;
    
    cndline = 2;
    base_image_itm = 76;
    clrchng_cnd = 1;
    %write images with number_clrchng_trials_btwn_images trials in between
    for block = 1:number_of_images/image_spacing
        for nov_rep = 1:2; %write 2x for novel then for repeat presentation
            for imgpair = 1:image_spacing/2;
                imgnums = [imgpair*2-1 imgpair*2]+(block-1)*image_spacing;
                if cndline < 10
                    cndspace = '  ';
                elseif cndspace < 100
                    cndspace = ' ';
                else
                    cndspace = '';
                end
                str = [cndspace num2str(cndline) '     -2    2      -4                         ' ...
                    num2str(imgnums(1)+base_image_itm-1) '\r\n'];
                fprintf(fid,str);
                cndline = cndline+1;
                for clrchng = 1:number_clrchng_trials_btwn_images
                    if cndline < 10
                        cndspace = '  ';
                    elseif cndline < 100
                        cndspace = ' ';
                    else
                        cndspace = '';
                    end
                    btfc = '     -3    1  '; %background timing, fixid, color palate
                    teststr = [];
                    for t = 1:3;
                        if t == 2
                            if all_clr_items(clrchng_cnd,t) < 10;
                                testspace = '                          ';
                            else
                                testspace = '                         ';
                            end
                        else
                            if all_clr_items(clrchng_cnd,t) < 10;
                                testspace = '     ';
                            else
                                testspace = '    ';
                            end
                        end
                        teststr = [teststr testspace num2str(all_clr_items(clrchng_cnd,t))];
                    end
                    fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                    cndline = cndline+1;
                    clrchng_cnd = clrchng_cnd+1;
                end
                str = [cndspace num2str(cndline) '     -2    2      -4                         ' ...
                    num2str(imgnums(2)+base_image_itm-1) '\r\n'];
                fprintf(fid,str);
                cndline = cndline+1;
                for clrchng = 1:number_clrchng_trials_btwn_images
                    if cndline < 10
                        cndspace = '  ';
                    elseif cndline < 100
                        cndspace = ' ';
                    else
                        cndspace = '';
                    end
                    btfc = '     -3    1  '; %background timing, fixid, color palate
                    teststr = [];
                    for t = 1:3;
                        if t == 2
                            if all_clr_items(clrchng_cnd,t) < 10;
                                testspace = '                          ';
                            else
                                testspace = '                         ';
                            end
                        else
                            if all_clr_items(clrchng_cnd,t) < 10;
                                testspace = '     ';
                            else
                                testspace = '    ';
                            end
                        end
                        teststr = [teststr testspace num2str(all_clr_items(clrchng_cnd,t))];
                    end
                    fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                    cndline = cndline+1;
                    clrchng_cnd = clrchng_cnd+1;
                end
            end
        end
    end
    fclose(fid);
end