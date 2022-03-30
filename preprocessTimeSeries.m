function [save_fn] = preprocessTimeSeries(filelist,template,register_flag,nonrigid_flag,)
% Modified version of A_ProcessTimeSeries, see original description below:
% Written 13Sep2021 KS
%
% First step in analysis pipeline:
% Extracts calcium imaging data from multi-page TIFs, performs X-Y
% registration, and calculates activity map (using local corr or modified
% DF/F) for further segmentation into ROIs
%
% Input:
% filelist: 1xN Cell array containing filenames as strings. If no input is
% given, user will be prompted to load file from current folder
%
% Multiple files can be batch processed, but files *must* belong to same
% imaging plane. If you wish to process multiple files from different
% imaging planes, run function within a loop.
%
% register_flag: (0,1) Whether to perform image registration (recommended,
% but increases processing time)
%
% nonrigid_flag: Whether to warp seperate files together using nonrigid
% registration. Use for aligning corresponding recordings from different
% sessions. Note: Computer Vision Systems toolbox required for this step.
%
% movie_flag: (0,1) Whether to save an avi movie file of time series
%
% Output: None - Data will be saved to a MAT file.
%
% Code written by Michael Goard - updated: Oct 2016
% Nonrigid registration added by James Roney - June 2017

%% Load files
if nargin==0
    filelist = uigetfile('.tif','MultiSelect','on');
    nonrigid_flag = 'No';
    register_flag = questdlg('Perform X-Y registration (movement correction)?','Dialog box','Yes','No','Yes');
elseif nargin ==1
    register_flag = questdlg('Perform X-Y registration (movement correction)?','Dialog box','Yes','No','Yes');
    if iscell(filelist)
        nonrigid_flag = questdlg('Perform Nonrigid Registration (for aligning seperate recordings)?','Dialog box','Yes','No','Yes');
    end
    movie_flag = questdlg('Create movie?','Dialog box','Yes','No','Yes');
end

if(strcmp(nonrigid_flag, 'Yes'))
    i = listdlg('PromptString', 'Select master recording:','SelectionMode','single','ListString',filelist);
    master = filelist{i};
    filelist{i} = filelist{1};
    filelist{1} = master; %shift master file to be in 1st index
end
%% Calculate length of file list
if iscell(filelist)==0
    lengthList = 1;
else
    lengthList = length(filelist);
end

%% Extract time series parameters from TIF file and save to data structure
for i = 1:lengthList
    
    %% Determine filename
    if iscell(filelist)==0
        data.filename = filelist;
    else
        data.filename = filelist{i};
    end
    
    %% Attempt to extract number of frames and image size from TIF header
    % If TIF header cannot be read, prompt user for values
    try
        header = imfinfo(data.filename);
        data.numFrames = length(header);
        data.xPixels = header(1).Width;
        data.yPixels = header(1).Height;
    catch
        data.numFrames = input('Enter number of frames: ');
        data.xPixels = input('Enter number of pixels (X dimension): ');
        data.yPixels = input('Enter number of pixels (Y dimension): ');
    end
    data.map_type = 'kurtosis';
    
    %% Attempt to extract frame rate from cfg file (PrairieView format)
    % If cfg file cannot be read, prompt user for frame rate
    directory = dir;
    for j = 1:length(directory)
        if ~isempty(strfind(directory(j).name,'env')) && sum(cfg_idx)<1
            cfg_idx(j) = 1;
        else
            cfg_idx(j) = 0;
        end
    end
    if sum(cfg_idx)==1
        cfg_filename = directory(cfg_idx==1).name;
    elseif sum(cfg_idx)>1
        cfg_filename = directory(cfg_idx==1).name;
    else
        cfg_filename = [];
    end
    if ~isempty(cfg_filename)
        cfg_file = importdata(cfg_filename);
        for j = 1:length(cfg_file)
            if strfind(cfg_file{j},'repetitionPeriod') > 0
                cfg_line = cfg_file{j};
                index = strfind(cfg_line,'repetitionPeriod');
                data.frameRate = 1/sscanf(cfg_line(index:end),'repetitionPeriod="%f"');
                if isinf(data.frameRate)
                    data.frameRate = input('Enter frame rate (Hz): ');
                end
            end
        end
    else
        data.frameRate = input('Enter frame rate (Hz): ');
    end
    
    %% Perform X-Y Registration (movement correction)
    if strcmp(register_flag,'Yes')
        disp('Performing X-Y registration...')
        % if i==1 || strcmp(nonrigid_flag,'Yes') 
        %     reference = [];
        % end

	reference = template; 
        overwrite = 1;
        maxOffset = 100;
        [new_filename] = registerStack(data,reference,0,maxOffset,1);
        data.filename = new_filename;
        if strcmp(nonrigid_flag,'Yes')
            data.filename = [data.filename(1:end-4) '_warped.tif'];
        end
        if lengthList==1
            filelist = data.filename;
        else
            filelist{i} = data.filename;
        end
    end
    
    %% Calculate kurtosis map
    
    % mean frame flourescence (for measuring photobleaching)
    frame_F = zeros(1,data.numFrames);
    
    % Params
    win = 7;             % in pixels (e.g., default: 7x7 neighborhood)
    block_size = 1000;    % frames in each block (to avoid memory crash)
    
    % Initalize
    num_blocks = ceil(data.numFrames/block_size);
    k_xy = zeros(data.yPixels,data.xPixels); % activity map
    avg_projection = zeros(data.yPixels,data.xPixels);
    frame_F = [];
    disp('Calculating activity map...')
    
    for j = 1:num_blocks
        idx_vec = (j-1)*block_size+1:min(j*block_size,data.numFrames);
        curr_block_size = length(idx_vec);
        tc = zeros(data.yPixels,data.xPixels,curr_block_size, 'single');
        
        for n = 1:curr_block_size
            tc(:,:,n) = single(imread(data.filename,idx_vec(n)));
        end 
        
        % average projection
        avg_projection = avg_projection+sum(tc,3);
        
        % mean fluorescence across frame
        frame_F = cat(2,frame_F,squeeze(mean(mean(tc,1),2))');
        
        % filter
        kernel = ones(win,win); % rectangular kernal
        for n = 1:curr_block_size
            tc(:,:,n) = imfilter(tc(:,:,n),kernel,'same')/sum(sum(kernel));
        end
        
        % Kurtosis
        k_image = kurtosis(tc,[],3);
        
        % Add to running average
        k_xy = k_xy+k_image*curr_block_size;
    end
    
    data.avg_projection = avg_projection/data.numFrames;
    data.frame_F = frame_F;
    reference = data.avg_projection;
    
    %% plot fluorescence over time (To look for photobleaching)
    if lengthList == 1
        plot(frame_F,'linewidth',2)
        F_fit = polyfit(1:data.numFrames,frame_F,1);
        PhotoBl = round((F_fit(1)*data.numFrames)/F_fit(2)*100);
        ylim([0 max(frame_F)*1.25])
        title(['Fluorescence timecourse, Photobleaching = ' num2str(PhotoBl) '%'])
        xlabel('Frame #')
        ylabel('Raw fluorescence')
        set(gcf,'color',[1 1 1])
        saveas(gcf,'Fluorescence_timecourse')
        close
    end
    
    % Divide activity map by frame number and save to data structure
    k_xy = k_xy/data.numFrames;
    k_xy(isnan(k_xy)) = 0;
    data.activity_map = k_xy;
                                                                              
    save_fn = strcat(data.filename(1:end-4), '_data.mat')
    save(save_fn, 'data');
    clear data
end

function [new_filename] = registerStack(data,reference,overwrite,maxOffset,use_fft)

%% subroutine_registerStack.m
%
% Register a stack to the average projection using Matlab imaging processing
% toolbox.
%
% Note: If offset is greater than value of 'maxOffset' value, then the
%       previous frame will be used (to prevent large fluorescence
%       artifacts), and an error message will be posted.
%       If there are many missed frames, the data may be of poor
%       quality or have excessive movement
%
% Called by A_ProcessTimeSeries.m
%
% Code written by Michael Goard - updated: Oct 2016

%% Default parameters
if nargin==4
    use_fft = 1;
end
if nargin==3
maxOffset = 10; % maximum pixel offset (Default = 12 for 16x zoom)
end

if nargin==1
    reference = [];
    overwrite = 0;
elseif nargin==2
    overwrite = 0;
end

%% Load file and image parameters
filename = data.filename;
numFrames = data.numFrames;
xPixels = data.xPixels;
yPixels = data.yPixels;
image_matrix = zeros([yPixels xPixels numFrames],'uint16');
data.offsets = zeros(data.numFrames, 2);

%% Make reference frame
if isempty(reference)
    disp('generating target frame...');
    if numFrames < 1000
        idx_vec = 1:numFrames;
    else
        idx_vec = [1:50];%100:numFrames];
    end
    sum_proj = zeros(yPixels, xPixels, length(idx_vec), 'single');
    for i = 1:length(idx_vec)
        sum_proj(:,:,i) = single(imread(data.filename,idx_vec(i)));
    end
    template = mean(sum_proj,3);
    close all
else
    template = reference;
end
ref_frame = zeros(size(template,1)+2*maxOffset,size(template,2)+2*maxOffset);
ref_frame(1+maxOffset:yPixels+maxOffset,1+maxOffset:xPixels+maxOffset) = template;

is_window = bwmorph(logical(template), 'close');
r = LargestRectangle(is_window);
x_edges = r(2:3, 1);
y_edges = r(3:4, 2);

%% Register
jump_ct = 0;
disp('aligning frames...');
new_filename = [filename(1:end-4) '_registered.tif'];
for i = 1:numFrames
    curr_frame = single(imread(filename,i));
    % Measure 2D xCorr

    if(use_fft)
        shifts = subroutine_dftregistration(fft2(template(y_edges(1):y_edges(2), x_edges(1):x_edges(2))), fft2(curr_frame(y_edges(1):y_edges(2), x_edges(1):x_edges(2))));
        corr_offset = -shifts(3:4);
    else
        cc = normxcorr2(template(y_edges(1):y_edges(2), x_edges(1):x_edges(2)), curr_frame(y_edges(1):y_edges(2), x_edges(1):x_edges(2)));
        [~,imax] = max(abs(cc(:)));
        [ypeak,xpeak] = ind2sub(size(cc),imax(1));
        corr_offset = [(ypeak-yPixels) (xpeak-xPixels)];
    end
    data.offsets(i,:) = corr_offset;
    % Determine offset and register
    if sum(abs(corr_offset)<maxOffset==[1 1])==2
        reg_frame = ref_frame;
        y_vec = 1+maxOffset-corr_offset(1):yPixels+maxOffset-corr_offset(1);
        x_vec = 1+maxOffset-corr_offset(2):xPixels+maxOffset-corr_offset(2);
        reg_frame(y_vec,x_vec) = curr_frame;
        reg_frame = reg_frame(1+maxOffset:yPixels+maxOffset,1+maxOffset:xPixels+maxOffset);
        prev_frame = reg_frame;
    else % If offset is greater than 'maxOffset', use previous frame
        disp(['Error: Frame ' num2str(i) ' offset greater than ' num2str(maxOffset) ' pixels, no registration performed'])
        if i>1
            reg_frame = prev_frame;
        else
            reg_frame = curr_frame;
            prev_frame = reg_frame;
        end
        jump_ct = jump_ct + 1;
    end
    image_matrix(:,:,i) = uint16(reg_frame);
    
    if jump_ct > floor(numFrames / 10)
	    keyboard
        error('Too many jumps, terminated registration.')
    end
end
close all

% write to Tif file (tif files >4GB supported)
disp('Writing to multi-page Tif file...')
options.big = true;
subroutine_saveastiff(image_matrix,new_filename,options);  

% delete old file
if overwrite==1
    disp('Deleting unregistered Tif file...')
    try imread(new_filename,numFrames);
        delete(filename)
        disp('Conversion complete.')
    catch disp('Error: new file appears not to have saved properly, unregistered Tif file preserved')
    end
else
    disp('Conversion complete.')
end
