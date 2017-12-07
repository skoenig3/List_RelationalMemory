% Written by Seth König 7/13/15
% set(0,'DefaultFigureVisible','OFF');
screen_size = get(0, 'ScreenSize');
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Cortex Data\';
fixwin = 3.5; %size of the fixation window/2. Was a width of 7 dva
imageX = 800;
imageY = 600;

pause(10)
% cortex_files = {'PW150416.4','PW150417.2','PW150420.2','PW150421.2',...
%     'PW150422.2','PW150423.2','PW150424.2','PW150427.2',...
%     'PW150428.2','PW150429.2','PW150430.2','PW150501.2',...
%     'PW150504.2','PW150505.2','PW150506.2','PW150507.2',...
%     'PW150511.2'};
% 
% item_files = {'ListRM21.itm','ListRM22.itm','ListRM23.itm','ListRM24.itm',...
%     'ListRM25.itm','ListRM26.itm','ListRM27.itm','ListRM28.itm',...
%     'ListRM29.itm','ListRM30.itm','ListRM31.itm','ListRM32.itm',...
%     'ListRM33.itm','ListRM34.itm','ListRM35.itm','ListRM36.itm',....
%     'ListRM37.itm'};

% cortex_files = {'RR150423.2','RR150424.2','RR150427.2','RR150428.2','RR150429.2',...
%     'RR150430.2','RR150501.2','RR150504.2','RR150505.2','RR150506.2',...
%     'RR150507.2','RR150508.2','RR150511.2','RR150512.2','RR150513.2',...
%     'RR150515.2','RR150518.2','RR150519.2','RR150520.2','RR150521.2'};
% item_files = {'ListRM01.itm','ListRM02.itm','ListRM03.itm','ListRM04.itm','ListRM05.itm',...
%     'ListRM06.itm','ListRM07.itm','ListRM08.itm','ListRM09.itm','ListRM10.itm',...
%     'ListRM11.itm','ListRM12.itm','ListRM13.itm','ListRM14.itm','ListRM15.itm',...
%     'ListRM16.itm','ListRM17.itm','ListRM18.itm','ListRM19.itm','ListRM20.itm'};

% cortex_files = {'TO150513.2','TO150514.2','TO150515.2','TO150518.2','TO150519.2',...
%     'TO150520.2','TO150521.2','TO150522.2','TO150526.2','TO150527.2',...
%     'TO150528.2','TO150529.2','TO150601.2','TO150602.2','TO150603.2',...
%     'TO150604.2','TO150605.2'};
% item_files = {'ListRM01.itm','ListRM02.itm','ListRM03.itm','ListRM04.itm','ListRM05.itm',...
%     'ListRM06.itm','ListRM07.itm','ListRM08.itm','ListRM09.itm','ListRM10.itm',...
%     'ListRM11.itm','ListRM12.itm','ListRM13.itm','ListRM14.itm','ListRM15.itm',...
%     'ListRM16.itm','ListRM17.itm'};

% cortex_files = {'TT150427.2','TT150428.2','TT150429.2','TT150430.2','TT150501.2',...
%     'TT150504.2','TT150505.2','TT150506.2','TT150507.2','TT150508.2',...
%     'TT150511.2','TT150512.2','TT150513.2','TT150514.2','TT150515.2',...
%     'TT150518.2'};
% item_files = {'ListRM01.itm','ListRM02.itm','ListRM03.itm','ListRM04.itm','ListRM05.itm',...
%     'ListRM06.itm','ListRM07.itm','ListRM08.itm','ListRM09.itm','ListRM10.itm',...
%     'ListRM11.itm','ListRM12.itm','ListRM13.itm','ListRM14.itm','ListRM15.itm',...
%     'ListRM17.itm'};
% 
% cortex_files = {'TO170404.2','TO170405.2','TO170406.2','TO170407.2',...
%                 'TO170410.2','TO170411.2','TO170412.2','TO170414.2',...
%                 'TO170419.2','TO170421.2','TO170425.2','TO170426.2',...
%                 'TO170427.2','TO170428.2','TO170501.2'};
% item_files = {'ListRM22.itm','ListRM23.itm','ListRM24.itm','ListRM25.itm',...
%               'ListRM26.itm','ListRM27.itm','ListRM28.itm','ListRM30.itm',...
%               'ListRM33.itm','ListRM35.itm','ListRM36.itm','ListRM37.itm',...
%               'ListRM38.itm','ListRM39.itm','ListRM40.itm'};

% cortex_files = {'TT150515.2'};%,...
% item_files = {'ListRM15.itm'};%,...

cortex_files = {'TO170404.2','TO170405.2','TO170406.2','TO170407.2',...
                'TO170410.2','TO170411.2','TO170412.2','TO170414.2',...
                'TO170419.2','TO170421.2','TO170425.2','TO170426.2',...
                'TO170427.2','TO170428.2','TO170501.2'};
item_files = {'ListRM22.itm','ListRM23.itm','ListRM24.itm','ListRM25.itm',...
              'ListRM26.itm','ListRM27.itm','ListRM28.itm','ListRM30.itm',...
              'ListRM33.itm','ListRM35.itm','ListRM36.itm','ListRM37.itm',...
              'ListRM38.itm','ListRM39.itm','ListRM40.itm'};


for file = 1:length(cortex_files)
    load([data_dir cortex_files{file}(1:8) '_' cortex_files{file}(end) '-fixation.mat'])
    
    viewed = zeros(2,96);
    img_dur = zeros(2,96);
    for trial = 1:length(per);
        [r,c] = find(imgs == per(trial).cnd-1000);
        if viewed(r,c) ~= 0
            error('This should be the first viewing of this image, but it is not!')
        else
            viewed(r,c) = trial;
            img_dur(r,c) = per(trial).alltim(per(trial).allval == 24)-...
                per(trial).alltim(per(trial).allval == 23);
        end
    end
    
    %1st image trial was displayed with the wrong timing file, then when fixed
    %started the task over again. The first image therefore is not novel and thus
    %we will not be analyzing data for this image on this session.
    if strcmpi('PW150416_2',cortex_files{file}(1:10))
        viewed(:,1) = 0;
        imgdur(:,1) = 0;
    end
    
    %Let's not even look at images in which the monkey looked away more than the
    %image was displayed for. Also should take care of EOG overflow trials
    viewed(viewed & img_dur > 14000) = 0; %1000 ms buffer for accidental lookaways and cortex timing
    %     viewed(:,viewed(1,:) == 0 | viewed(2,:) == 0) = [];
    
    setnum = str2double(item_set(7:8));
    if setnum < 10
        imgdir = ['C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Image Sets\LRM0' num2str(setnum) '\'];
    else
        imgdir = ['C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Image Sets\LRM' num2str(setnum) '\'];
    end
    
    figure_dir = ['C:\Users\seth.koenig\Documents\MATLAB\List_RelationalMemory\Figures\' num2str(setnum) '\'];
    mkdir(figure_dir)
    
    
    for img = 1:96
        if viewed(1,img) ~= 0 && viewed(2,img) ~= 0
            
            nov_allval = per(viewed(1,img)).allval;
            nov_alltim = per(viewed(1,img)).alltim;
            nov_img_on = nov_alltim(nov_allval == 23)-nov_alltim(nov_allval == 100);%image on relative to eye data start
            nov_img_off = nov_alltim(nov_allval == 24)-nov_alltim(nov_allval == 100);%image on relative to eye data start
            
            nov_x = fixationstats{viewed(1,img)}.XY(1,nov_img_on:nov_img_off);
            nov_y = fixationstats{viewed(1,img)}.XY(2,nov_img_on:nov_img_off);
            nov_pupil = pupildata{viewed(1,img)}(round(nov_img_on/5):round(nov_img_off/5));
            
            nov_fix = fixationstats{viewed(1,img)}.fixations;
            nov_fixtimes = fixationstats{viewed(1,img)}.fixationtimes;
            nov_sactimes = fixationstats{viewed(1,img)}.saccadetimes;
            
            pre_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
            post_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
            pre_img_saccades = find(nov_sactimes(1,:) <= nov_img_on);
            post_img_saccades = find(nov_sactimes(2,:) > nov_img_off);
            
            nov_fix(:,post_img_fix) = [];
            nov_fix(:,pre_img_fix) = [];
            nov_fixtimes(:,post_img_fix) = [];
            nov_fixtimes(:,pre_img_fix ) = [];
            nov_sactimes(:,post_img_saccades) = [];
            nov_sactimes(:,pre_img_saccades) = [];
            
            
            %for repeat images
            rep_allval = per(viewed(2,img)).allval;
            rep_alltim = per(viewed(2,img)).alltim;
            rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
            rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
            
            rep_x = fixationstats{viewed(2,img)}.XY(1,rep_img_on:rep_img_off);
            rep_y = fixationstats{viewed(2,img)}.XY(2,rep_img_on:rep_img_off);
            rep_pupil = pupildata{viewed(2,img)}(round(rep_img_on/5):round(rep_img_off/5));
            
            rep_fix = fixationstats{viewed(2,img)}.fixations;
            rep_fixtimes = fixationstats{viewed(2,img)}.fixationtimes;
            rep_sactimes = fixationstats{viewed(2,img)}.saccadetimes;
            
            pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
            post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
            pre_img_saccades = find(rep_sactimes(1,:) <= rep_img_on);
            post_img_saccades = find(rep_sactimes(2,:) > rep_img_off);
            
            rep_fix(:,post_img_fix) = [];
            rep_fix(:,pre_img_fix) = [];
            rep_fixtimes(:,post_img_fix) = [];
            rep_fixtimes(:,pre_img_fix) = [];
            rep_sactimes(:,post_img_saccades) = [];
            rep_sactimes(:,pre_img_saccades) = [];
            
            %---Remove Central Fixations that continued to Occur once image was on---%
%             first_out = find((nov_fix(1,1) > imageX/2+fixwin*24 | nov_fix(1,1) <  imageX/2-24*fixwin) | ...
%                 (nov_fix(2,1) > imageY/2+fixwin*24 | nov_fix(2,1) < imageY/2-fixwin*24));
%             nov_fix(:,1) = [];
%             nov_fixtimes(:,1) = [];
%             
%             first_out = find((rep_fix(1,1) > imageX/2+fixwin*24 | rep_fix(1,1) <  imageX/2-24*fixwin) | ...
%                 (rep_fix(2,1) > imageY/2+fixwin*24 | rep_fix(2,1) < imageY/2-fixwin*24));
%             rep_fix(:,1) = [];
%             rep_fixtimes(:,1) = [];
            
            
            if setnum < 10
                setzero = '0';
            else
                setzero = '';
            end
            if img < 10
                imgzero = '0';
            else
                imgzero = '';
            end
            imag = imread([imgdir 'S' setzero num2str(setnum) 'I' imgzero num2str(img) '.bmp']);
            
            figure
%             subplot(4,4,[1,2,5,6])
            subplot(2,2,1)
            imshow(imag);
            hold on
            plot(nov_x,600-nov_y,'b');
            nov_x = fixationstats{viewed(1,img)}.XY(1,:);
            nov_y = fixationstats{viewed(1,img)}.XY(2,:);
            for s = 1:size(nov_sactimes,2)
                plot(nov_x(nov_sactimes(1,s):nov_sactimes(2,s)),600-nov_y(nov_sactimes(1,s):nov_sactimes(2,s)),'g')
            end
            for f = 1:size(nov_fixtimes,2)
                plot(nov_x(nov_fixtimes(1,f):nov_fixtimes(2,f)),600-nov_y(nov_fixtimes(1,f):nov_fixtimes(2,f)),'r')
            end
            hold off
            axis off
            title('Novel')
            
%             subplot(4,4,[3,4,7,8])
            subplot(2,2,2)

            imshow(imag);
            hold on
            plot(rep_x,600-rep_y,'b');
            rep_x = fixationstats{viewed(2,img)}.XY(1,:);
            rep_y = fixationstats{viewed(2,img)}.XY(2,:);
            for s = 1:size(rep_sactimes,2)
                plot(rep_x(rep_sactimes(1,s):rep_sactimes(2,s)),600-rep_y(rep_sactimes(1,s):rep_sactimes(2,s)),'g')
            end
            for f = 1:size(rep_fixtimes,2)
                plot(rep_x(rep_fixtimes(1,f):rep_fixtimes(2,f)),600-rep_y(rep_fixtimes(1,f):rep_fixtimes(2,f)),'r')
            end
            hold off
            axis off
            title('Repeat')
            
%             subplot(4,4,[9 13])
            subplot(2,2,3)

            plot(nov_x,'b');
            hold on
            for s = 1:size(nov_sactimes,2)
                plot(nov_sactimes(1,s):nov_sactimes(2,s),nov_x(nov_sactimes(1,s):nov_sactimes(2,s)),'g')
            end
            for f = 1:size(nov_fixtimes,2)
                plot(nov_fixtimes(1,f):nov_fixtimes(2,f),nov_x(nov_fixtimes(1,f):nov_fixtimes(2,f)),'r')
            end
            hold off
            axis off
            title('Novel')
            
            %subplot(4,4,[12 16])
                        subplot(2,2,4)

            plot(rep_x,'b');
                        hold on
            for s = 1:size(rep_sactimes,2)
                plot(rep_sactimes(1,s):rep_sactimes(2,s),rep_x(rep_sactimes(1,s):rep_sactimes(2,s)),'g')
            end
            for f = 1:size(rep_fixtimes,2)
                plot(rep_fixtimes(1,f):rep_fixtimes(2,f),rep_x(rep_fixtimes(1,f):rep_fixtimes(2,f)),'r')
            end
            hold off
            axis off
            title('repeat')
            
%             subplot(4,4,[14:15])
%             plot(nov_fixtimes(2,:)-nov_fixtimes(1,:)+1,'b')
%             hold on
%             plot(rep_fixtimes(2,:)-rep_fixtimes(1,:)+1,'r')
%             hold off
%             ylim([50 500])
            
            subtitle(['ListRM image# ' num2str(img) ' ' cortex_files{1}(1:2)]);
%             set(gcf, 'Position', [0 0 screen_size(3) screen_size(4)]);
            %close
            save_and_close_fig(figure_dir,['ListRM_' cortex_files{1}(1:2) '_' num2str(setnum) '_' num2str(img)])
        end
    end
end