addpath 'D:\_NeurotarRecordings\_HeadRotation\_code'
% run subroutine_tifConvert first...
% subroutine_tifConvert();

sr = SpinnyRegisterer('SKKS091-HeadRotation-001.tif');
% figure; % for display
% sr.dirtyRegister();
tmp = importdata('SKKS091-HeadRotation-001_supplementary_files.mat');
sr.cleanRegister(tmp);
% sr.cleanUp();

% crop time
sr.crop(tmp.occupancy, './cropped/', 'out');

img = 0;
for ii = 1:2550
    temp = double(imread('SKKS091-HeadRotation-001_derotated_cropped.tif',ii));
    imagesc(temp)
    img = img + temp;
    pause(1/60)
end

sum_proj = img;

heading = atan2d(squeeze(tform_mat(2, 1, :)), squeeze(tform_mat(2, 2, :)));

% Changed this to be bin_width of 3, and a smoothing filter, a la Giocomo et al 2014, Curr Bio

dat = data.s2p_spks;
bin_width = 3;
bin_edges = -180:bin_width:180;
groups = discretize(heading, bin_edges);

u_groups = 1:length(bin_edges) - 1; % Get all the possible groups
out = zeros(size(dat, 1), length(u_groups));

for g = 1:length(u_groups)
    for c = 1:size(dat, 1)
        temp = dat(c, groups == u_groups(g));
        out(c, g) = mean(temp);
    end
end

out = movmean(out, 30/3, 2); % 15 degree smoothing filter (as in Giocomo et al 2014)

for ii = 1:size(out, 1)
    plot(out(ii, :))
    pause
end