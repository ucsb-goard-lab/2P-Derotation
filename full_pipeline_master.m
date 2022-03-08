% Step 0: Set up
addpath(genpath(pwd)); % Make sure you add the entire folder and subfolders to have all the necessary helper functions

% Step 1: Choose your recordings for head and platform rotations
[head_fn, head_pn] = uigetfile('*.tif');
[plat_fn, plat_pn] = uigetfile('*.tif');

heading_recording = strcat(head_pn, head_fn);
platform_recording = strcat(plat_pn, plat_fn);

% Step 2: Perform standard stuff on the platform recording to generate a
% template.

cd(plat_pn);
plat_matfile_fn = A_ProcessTimeSeries_SpinnyMod(platform_recording, [], 'Yes', 'No', 'No');
plat_data = importdata(plat_matfile_fn);

platform_recording = plat_data.filename;

% Step 2: Perform rotation correction (derotation) on the head rotation (not necessary for platform rotation)
sr = SpinnyRegisterer(heading_recording, plat_data.avg_projection);

supplementary_files = dir(strcat(head_pn, '*_supplementary_files.mat'));
if isempty(supplementary_files)
    prev = cd(head_pn);
    sr.dirtyRegister2();
    sr.cleanRegister();
    sr.cleanUp();
    cd(prev)
else
    sr.load(strcat(supplementary_files.folder, '\\', supplementary_files.name))
end
sr.save();

% Step 3: Coregister the platform rotation recording to the template from head rotation
temp = imadjust(rescale(plat_data.avg_projection));

% Pad it
autoregister_flag = false;
if autoregister_flag
    [optimizer, metric] = imregconfig('multimodal');
    tform = imregtform(temp, sr.template, 'rigid', optimizer, metric);
    output_view = affineOutputView(size(sr.template), tform, 'BoundsStyle', 'SameAsInput');
    
else
    [optimizer, metric] = imregconfig('multimodal');
    tform = imregtform(temp, sr.template, 'rigid', optimizer, metric);
    output_view = affineOutputView(size(sr.template), tform, 'BoundsStyle', 'SameAsInput');
    
        registrationEstimator(temp, sr.template);
    disp('Press any key to continue: ')
    pause
    output_view = movingReg.SpatialRefObj;
    tform = movingReg.Transformation;
end

out = imwarp(plat_data.avg_projection, tform, 'OutputView', output_view);
imshowpair(out, sr.template)

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

% Step 4: Crop the platform rotation and head direction to the same FOV
cd(fileparts(heading_recording))
head_cropped_fn = cropToCenter(sr.output_filename, sr.occupancy, strcat(fileparts(heading_recording))); 

cd(fileparts(platform_recording))
plat_cropped_fn = cropToCenter(new_fn, sr.occupancy, strcat(fileparts(platform_recording)));

% Step 5: Pass the recordings through suite2p_coregister.py to get the same 
% probably going to have to open up anaconda prompt to take care of this manually...
cd(fileparts(platform_recording))
plat_matfile_fn = A_ProcessTimeSeries_SpinnyMod(plat_cropped_fn, [], 'Yes', 'No', 'No');
plat_data = importdata(plat_matfile_fn);


% use the platform registered recording to register to head one
cd(fileparts(heading_recording))
head_matfile_fn = A_ProcessTimeSeries_SpinnyMod(head_cropped_fn, plat_data.avg_projection, 'Yes', 'No', 'No');

% check here for good alignment
head_data = importdata(head_matfile_fn);

imshowpair(rescale(plat_data.avg_projection), rescale(head_data.avg_projection), 'montage');

% adjusting the activity maps to get rid of the ring of high activity?
data = importdata(plat_matfile_fn);
window_mask = logical(data.activity_map);
data.activity_map(~window_mask) = NaN; % set outside to NaNs;
window_mask = bwmorph(window_mask, 'remove');
window_mask = bwmorph(window_mask, 'dilate', 5);
data.activity_map(window_mask) = NaN;
data.activity_map(isnan(data.activity_map)) = min(data.activity_map(:));
save(plat_matfile_fn, 'data');

% head?
data = importdata(head_matfile_fn);
median_val = median(data.activity_map(:));
std_val = std(data.activity_map(:));
data.activity_map(data.activity_map > median_val + 5 * std_val) = 0;
save(head_matfile_fn, 'data');

% Define ROIs based on the platform recording (clearer?)
B_DefineROI(plat_matfile_fn)

% transfer ROIs
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
