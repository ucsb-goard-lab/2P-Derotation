function save_fn = cropToCenter(fn, occupancy, save_path)

if nargin < 3 || isempty(save_path)
    save_path = pwd;
end
if ~exist(save_path)
    mkdir(save_path)
end

image_info = imfinfo(fn);
occupancy(occupancy < max(occupancy(:))) = 0; % get only the major overlap areas
occupancy = bwmorph(occupancy, 'clean'); % clean up lone pixels

bbox_struct = regionprops(occupancy,'BoundingBox');

bbox = floor(bbox_struct.BoundingBox);
crop_mask = logical(occupancy(bbox(2) : bbox(2) + bbox(4), bbox(1) : bbox(1) + bbox(3)));
[~, fn_only] = fileparts(fn);save_fn = strcat(save_path, '/', fn_only, '_cropped.tif');

writer = Tiff(save_fn, 'w8'); % big tiff
fprintf('Writing cropped image...\n')
for im = 1:numel(image_info)
    img = imread(fn, im);
    img = img(bbox(2) : bbox(2) + bbox(4), bbox(1) : bbox(1) + bbox(3));
    img(~crop_mask) = NaN;
    writer = setDefaultTags(writer, size(img));
    writer.write(img); % write the image
    writer.writeDirectory(); % move on for multipage
end
writer.close();
end

function writer = setDefaultTags(writer, img_size)
writer.setTag('ImageWidth', img_size(2));
writer.setTag('ImageLength', img_size(1));
writer.setTag('Photometric', Tiff.Photometric.MinIsBlack);
writer.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
writer.setTag('Compression', Tiff.Compression.None);
writer.setTag('BitsPerSample', 16);
writer.setTag('SamplesPerPixel', 1);
end
