%% Step 0: Set up
addpath(genpath('E:\_Code\2P-Derotation')); % Make sure you add the entire folder and subfolders to have all the necessary helper functions

%% Step 1: Choose your recordings for head and platform rotations
[head_fn, head_pn] = uigetfile('*.tif');
[plat_fn, plat_pn] = uigetfile('*.tif');

if contains(head_fn, '000001.ome.tif') % not yet run through subroutine_tifConvert
    cd(head_pn)
    heading_recording = subroutine_tifConvert(head_fn);
else
    heading_recording = strcat(head_pn, head_fn);
end

if contains(plat_fn, '000001.ome.tif')
    cd(plat_pn)
    platform_recording = subroutine_tifConvert(plat_fn);
else
    platform_recording = strcat(plat_pn, plat_fn);
end

% correct weird slashes because of stupid windows
heading_recording(strfind(heading_recording, '\')) = '/';
platform_recording(strfind(platform_recording, '\')) = '/';

%% Step 2: Perform standard stuff on the platform recording to generate a
% template.

cd(plat_pn);
plat_matfile_fn = processTimeSeries(platform_recording);
plat_data = importdata(plat_matfile_fn);

platform_recording = plat_data.filename;

%% Step 3: Perform rotation correction (derotation) on the head rotation (not necessary for platform rotation)
dr = Derotater(heading_recording, plat_data.avg_projection);

supplementary_files = dir(strcat(head_pn, '*_supplementary_files.mat'));
if isempty(supplementary_files)
    prev = cd(head_pn);
    dr.dirtyRegister(); % first pass at registration
    dr.cleanRegister(); % adding some smoothing and fine grained adjustments
    dr.cleanUp();
    cd(prev)
else
    dr.load(strcat(supplementary_files.folder, '\\', supplementary_files.name))
end
dr.save();

%% Step 4: Coregister the platform rotation recording to the template from head rotation
temp = imadjust(rescale(plat_data.avg_projection));

% Pad it
autoregister_flag = false; % feel free to try autoregistration, but becaus ethere are few field to register, i usually manually do it
if autoregister_flag
    [optimizer, metric] = imregconfig('multimodal');
    tform = imregtform(temp, dr.template, 'rigid', optimizer, metric);
    output_view = affineOutputView(size(dr.template), tform, 'BoundsStyle', 'SameAsInput');
    
else
    % Recommended starting points:
    % Change registration type to "Intensity: Multimodal"
    % Check box for "Normalize"
    % Change Align Centers to "Center of Mass"
    % Set Transformation to "Rigid"
    % Try tweaking initial radius to improve fit
    
    % When finished, press "Export", and allow it to export your results
    % into a variable caleld "movingReg". Then press any key in the command
    % window to continue.
    [optimizer, metric] = imregconfig('multimodal');
    tform = imregtform(temp, dr.template, 'rigid', optimizer, metric);
    output_view = affineOutputView(size(dr.template), tform, 'BoundsStyle', 'SameAsInput');
    
    registrationEstimator(temp, dr.template);
    disp('Press any key to continue: ')
    pause
    output_view = movingReg.SpatialRefObj;
    tform = movingReg.Transformation;
end

out = imwarp(plat_data.avg_projection, tform, 'OutputView', output_view);
imshowpair(out, dr.template)

[folder, noxt] = fileparts(platform_recording);
new_fn = strcat(folder, '\', noxt, '_coregistered.tif');

writer = Tiff(new_fn, 'w8');
disp('Writing coregistered images...')
for ii = 1:numel(plat_data.frame_F)
    img = imread(platform_recording, ii);
    warp_img = imwarp(img, tform, 'OutputView', output_view);
    
    writer.setTag('ImageWidth', size(warp_img, 2));
    writer.setTag('ImageLength', size(warp_img, 1));
    writer.setTag('Photometric', Tiff.Photometric.MinIsBlack);
    writer.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    writer.setTag('Compression', Tiff.Compression.None);
    writer.setTag('BitsPerSample', 16);
    writer.setTag('SamplesPerPixel', 1);
    writer.write(warp_img); % write the image
    writer.writeDirectory(); % move on for multipage
end
writer.close();
clear('writer')

%% Step 5: Crop the platform rotation and head direction to the same FOV
cd(fileparts(heading_recording))
head_cropped_fn = cropToCenter(dr.output_filename, dr.occupancy, strcat(fileparts(heading_recording))); 

cd(fileparts(platform_recording))
plat_cropped_fn = cropToCenter(new_fn, dr.occupancy, strcat(fileparts(platform_recording)));

%% Step 6: Process recordings with A_ProcessTimeSeries
cd(fileparts(platform_recording))
plat_matfile_fn = processTimeSeries(plat_cropped_fn);
plat_data = importdata(plat_matfile_fn);

% use the platform registered recording to register to head one
cd(fileparts(heading_recording))
head_matfile_fn = processTimeSeries(head_cropped_fn, plat_data.avg_projection);

% check here for good alignment
head_data = importdata(head_matfile_fn);

imshowpair(rescale(plat_data.avg_projection), rescale(head_data.avg_projection), 'montage');

% Adjusting activity maps
% This is because the movement of the imaging field within the circle
% causes a "ring" of really high activity, which we remove here.
data = importdata(plat_matfile_fn);
window_mask = logical(data.activity_map);
data.activity_map(~window_mask) = NaN; % set outside to NaNs;
window_mask = bwmorph(window_mask, 'remove');
window_mask = bwmorph(window_mask, 'dilate', 5);
data.activity_map(window_mask) = NaN;
data.activity_map(isnan(data.activity_map)) = min(data.activity_map(:));
save(plat_matfile_fn, 'data');

data = importdata(head_matfile_fn);
median_val = median(data.activity_map(:));
std_val = std(data.activity_map(:));
data.activity_map(data.activity_map > median_val + 5 * std_val) = 0;
save(head_matfile_fn, 'data');

% Define ROIs based on the platform recording
B_DefineROI(plat_matfile_fn)

% transfer ROIs to the head rotation
donor = importdata(plat_matfile_fn);
cd(fileparts(heading_recording))
load(head_matfile_fn)
data.cellMasks = donor.cellMasks;
save(head_matfile_fn, 'data');

% Extract DFFs
cd(fileparts(heading_recording))
C_ExtractDFF(head_matfile_fn, 'Local Neuropil', 'Yes');

cd(fileparts(platform_recording))
C_ExtractDFF(plat_matfile_fn, 'Local Neuropil', 'Yes');

% Done.
disp('Finished.')
