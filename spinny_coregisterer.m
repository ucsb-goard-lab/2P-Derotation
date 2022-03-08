platform_fn = ('./SKKS091-PlatformRotation-002/SKKS091-PlatformRotation-002.tif');
head_fn = ('./SKKS091-HeadRotation-001/cropped/SKKS091-HeadRotation-001_derotated_cropped.tif');

% read imgs
platform_avg = 0;
for ii = 1:2550
    platform_avg = platform_avg + double(imread(platform_fn, ii));
end
    

load('./SKKS091-HeadRotation-001/SKKS091-HeadRotation-001_supplementary_files.mat');
[optimizer, metric] = imregconfig('multimodal');
% tform =imregtform(head.avg_projection, platform.avg_projection, 'rigid', optimizer, metric);
% imwarp(head.avg_projection, tform)


% tform_rough = imregcorr(platform_avg, template, 'rigid');
tform = imregtform(platform_avg, template, 'rigid', optimizer, metric);

output_view = affineOutputView(size(template), tform, 'BoundsStyle', 'SameAsInput');
out = imwarp(platform_avg, tform, 'OutputView', output_view);
imshowpair(out, template)


[~, noxt] = fileparts(platform_fn);
writer = Tiff(strcat('./cropped', noxt, '_coregistered.tif'), 'w8');
for ii = 1:2550
    disp(ii)
    img = imread(platform_fn, ii);
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