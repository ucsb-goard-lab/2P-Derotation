classdef SpinnyRegisterer < handle
	properties
		avg_proj
		filename
		output_filename
		width
		height
		n_images
		recorded_angle
		registration_type

		tform_mat
		template
		occupancy
		idx
	end

	methods
		function obj = SpinnyRegisterer(filename, avg_projection)
			if nargin < 2 || isempty(avg_projection)
				if obj.n_images > 2550
					avg_projection = obj.calculateAvgProjection2();
				else
					avg_projeciton = obj.calculateAvgProjection();
				end
			end

			fprintf('Initializing SpinnyRegisterer...\n')
			obj.filename = filename;
			obj.avg_proj = avg_projection;
			obj.getBasicImageInfo();
			obj.checkForAffine2DRelaxed();
			obj.checkForRecordedAngle();
			obj.output_filename = strcat(obj.removeExt(obj.filename), '_derotated.tif');
		end

		function dirtyRegister(obj, debug_stop_frame)
			if nargin < 2 || isempty(debug_stop_frame)
				debug_stop_frame = Inf;
			end
			% Main registration function

			% initialize matrices
			if sum(obj.tform_mat(:)) == 0;
				obj.tform_mat = zeros(3, 3, obj.n_images);     % matrix of transform coorindates
			end
			% pad so the template doesn't try to grow...
			obj.template = single(padarray(obj.avg_proj, round([size(obj.avg_proj)]/2))); % big boi
			obj.occupancy = ones(size(obj.template));

			% set up defaults for registering
			[optimizer, metric] = imregconfig('multimodal');

			% prepare for writing
			writer = Tiff(strcat(obj.removeExt(obj.filename), '_temporary.tif'), 'w8'); % big tiff

			% Start the main loop
			tic;
			is_jump = zeros(1, obj.n_images);
			im = 1;
			ct = 1;
			while im < obj.n_images
				obj.idx = im;
				tic;
				fprintf('Frame #%d\n', im);

				% read the next image
				raw_img = imread(obj.filename, im);
				processed_img = obj.preProcessImage(raw_img);
				if im >= debug_stop_frame
					keyboard
				end
				if im == 1 % first image might need extra tweaking
					tform = obj.registerFirstImage(processed_img);
				else

					best_guess = affine2d_relaxed(obj.getBestGuess());
					% initial attempt at registration
					tform = imregtform(processed_img, obj.template, obj.registration_type, optimizer, metric,'InitialTransformation', best_guess);

					% continue to try if it is not ready
					try_ct = 1;
					tform_store = []; % reset
					while ~obj.registrationReady(tform)
						[~, warn_id] = lastwarn();
						is_jump(im) = true;

						if ~isempty(warn_id)
							fprintf('Stepping down initial radius, current value: %0.02g new value: %0.02g\n', optimizer.InitialRadius, optimizer.InitialRadius/10);
							optimizer.InitialRadius = optimizer.InitialRadius/2;
						else
							% store previous transforms and try again
							tform_store(:, :, try_ct) = tform.T;
							try_ct = try_ct + 1;
						end

						if optimizer.InitialRadius < 1e-03 || try_ct == 5
							if obj.checkTriesForConsistency(tform_store) && ~obj.checkForJump(affine2d_relaxed(mean(tform_store, 3)))
								tform = affine2d_relaxed(mean(tform_store, 3));
								break
							else
								% give up
								fprintf('Unable to reach a good registration, applying previous transform\n')
								tform = affine2d_relaxed(obj.getBestGuess());
								break
							end
						end

						% test new transform parameters
						tform = imregtform(processed_img, obj.template, obj.registration_type, optimizer, metric, 'InitialTransformation', best_guess);
					end

					% reset optimizer
					optimizer.InitialRadius = 6.25e-03; % default
					tform = affine2d_relaxed(tform.T); % update to an affine2d_relaxed to override the issues with isRigid
				end
				% save tform values and registered image
				obj.tform_mat(:, :, im) = tform.T;
				% prepare output view
				output_view = affineOutputView(size(obj.template), tform, 'BoundsStyle', 'SameAsInput');

				% Apply transform
				warped_img = imwarp(raw_img, tform, 'OutputView', output_view);

				% for normalization
				current_warped_occupancy = imwarp(ones(size(raw_img)), tform, 'OutputView', output_view);
				obj.template = obj.template + single(warped_img); %re-added, 25Aug2021
				obj.occupancy = obj.occupancy + current_warped_occupancy;
				writer = obj.setDefaultTags(writer, [], true);
				writer.write(warped_img); % write the image
				writer.writeDirectory(); % move on for multipage

				time_vec(im) = toc;
				% visualization
				if im > 1
					best_guess_img = imwarp(raw_img, best_guess, 'OutputView', output_view);
					obj.visualize(raw_img, warped_img, best_guess_img, is_jump, time_vec)

				end

				if sum(is_jump(max(1, im - 10) : im)) > 5
					% if too many consecutive jumps, go back
					start_jump = find(is_jump(max(1, im - 10) : im), 1); % find the first jump
					is_jump(max(1, im - 10) : im) = false;

					im = im - (10 - start_jump); % go back to there
					ct = ct - (10 - start_jump);
				end
				im = im + 1;
				ct = ct + 1;
			end
			writer.close();
			obj.save()
			fprintf('Successfully registered in: %0.02f\n', toc)
		end

		function is_consistent = checkTriesForConsistency(obj, tform_store)
			angle_thresh = 10;
			x_thresh = 50;
			y_thresh = 50;

			guess_std = std(tform_store, [], 3);
			angle_std = std(acosd(squeeze(tform_store(1, 1, :))));
			x_std = guess_std(end, 1);
			y_std = guess_std(end, 2);

			if all([angle_std < angle_thresh, x_std < x_thresh, y_std < y_thresh])
				is_consistent = true;
			else
				is_consistent = false;
			end
		end

		function cleanRegister(obj, supplementary_info)
			if nargin < 2 || isempty(supplementary_info)
				tform_matrix = obj.tform_mat;
				template = obj.template;
			else
				tform_matrix = supplementary_info.tform_mat;
				template = supplementary_info.template;
			end

			if isempty(tform_matrix) || isempty(template)
				error('Need to either supply a transform matrix or have run the registerer beforehand')
			end

			smooth_angles = movmean(tform_matrix(1:2, 1:end, :), 5, 3);
			smooth_tform_mat = cat(1, smooth_angles, tform_matrix(end, 1:end, :));
			% warning('Line 176 clean registration, some weird issue with the last value of tform_matrix is all 0s, causing a bug, trace this back later')
			smooth_tform_mat(:, :, end) = tform_matrix(:, :, end-1); % not sure what happened here, fix later
			writer = Tiff(obj.output_filename, 'w8');

			for im = 1:obj.n_images
				raw_img = imread(strcat(obj.filename), im);
				tform = affine2d(smooth_tform_mat(:, :, im));
				output_view = affineOutputView(size(template), tform, 'BoundsStyle', 'SameAsInput');
				warp_img = imwarp(raw_img, tform, 'OutputView', output_view);
				writer = obj.setDefaultTags(writer, size(template), true);
				writer.write(warp_img);
				writer.writeDirectory();

				if rem(im,10)==0
					disp(['Frame ' num2str(im) '/' num2str(obj.n_images) ' complete.'])
				end
			end
			writer.close();
			disp('Finished cleaning rotation.')
		end

		function cleanUp(obj)
			disp('Cleaning up temporary files...')
			% check for properly read
			try
				imread(obj.output_filename, obj.n_images);
			catch ME
				if strcmp(ME.identifier, 'MATLAB:imagesci:rtifc:invalidDirIndex')
					error('Something went wrong with writing the tif file...')
				else
					error(sprintf(ME.message))
				end
			end
			delete(strcat(obj.removeExt(obj.filename), '_temporary.tif')); % get rid of teh temporary tif file
			disp('Finished cleaning up.')
		end

		function heading = recoverHeading(obj)
			heading = atan2d(squeeze(obj.tform_mat(2, 1, :)), squeeze(obj.tform_mat(2, 2, :)));
		end

		function load(obj, supplementary_info_fn)
			s_info = importdata(supplementary_info_fn);
			obj.template = s_info.template;
			obj.occupancy = s_info.occupancy;
			obj.tform_mat = s_info.tform_mat;
		end

		function out = removeExt(obj, filename)
			[~, out] = fileparts(obj.filename);
		end

		function crop(obj, occupancy, save_path, crop_style)
			if nargin < 2 || isempty(occupancy)
				occupancy = obj.occupancy;
			end

			if nargin < 3 || isempty(save_path)
				save_path = pwd;
			end

			if nargin < 4 || isempty(crop_style)
				crop_style = 'in';
			end

			if isempty(occupancy)
				error('Derotation not performed or loaded, cannot continue.')
			end
			if ~exist(save_path)
				mkdir(save_path)
			end

			occupancy(occupancy < max(occupancy(:))) = 0; % get only the major overlap areas
			occupancy = bwmorph(occupancy, 'clean'); % clean up lone pixels
			switch crop_style
				case 'in'
					r = obj.getLargestRectangle(occupancy);
					r_mask = poly2mask(r(2:end, 1), r(2:end, 2), size(occupancy, 1), size(occupancy, 2)); % convert to mask
					rotated_r_mask = imrotate(r_mask, r(1, 2), 'crop');
					bbox_struct = regionprops(rotated_r_mask, 'BoundingBox');
					occupancy = imrotate(occupancy, r(1, 2), 'crop');
				case 'out'
					bbox_struct = regionprops(occupancy,'BoundingBox');
			end
			bbox = floor(bbox_struct.BoundingBox);
			crop_mask = logical(occupancy(bbox(2) : bbox(2) + bbox(4), bbox(1) : bbox(1) + bbox(3)));
			[~, fn_only] = fileparts(obj.output_filename);
			writer = Tiff(strcat(save_path, fn_only, '_cropped.tif'), 'w8'); % big tiff
			fprintf('Writing cropped image...\n')
			for im = 1:obj.n_images
				img = imread(obj.output_filename, im);
				if strcmp(crop_style, 'in')
					img = imrotate(img, r(1, 2), 'crop');
				end
				img = img(bbox(2) : bbox(2) + bbox(4), bbox(1) : bbox(1) + bbox(3));
				img(~crop_mask) = 0;
				writer = obj.setDefaultTags(writer, size(img), false);
				writer.write(img); % write the image
				writer.writeDirectory(); % move on for multipage
			end
			writer.close();
		end

		function save(obj)
			% Saving some supplemetnal files to go along with the thing
			tform_mat = obj.tform_mat;
			template = obj.template;
			occupancy = obj.occupancy;
			heading = obj.recoverHeading();
			[pathname, no_ext] = fileparts(obj.filename);
			save_fn = strcat(pathname, '\', no_ext, '_supplementary_files.mat');
			save(save_fn, 'tform_mat', 'template', 'occupancy', 'heading')
		end

		function getBasicImageInfo(obj)
			image_info = imfinfo(obj.filename);
			obj.n_images = numel(image_info);
			obj.width = image_info(1).Width;
			obj.height = image_info(1).Height;
		end

		function out = calculateAvgProjection2(obj)
			fprintf('Calculating avg projection (front)...\n')
			img_mat = zeros(obj.height, obj.width, 50);
			for im = 1:50
				img_mat(:, :, im) = imread(obj.filename, im);
			end
			out = mean(img_mat, 3);
		end
		function out = calculateAvgProjection(obj)
			fprintf('Calculating average projection...\n')
			% usually the last 50 are not moving so..
			checkback = 49;
			img_mat = zeros(obj.height, obj.width, checkback+1);
			for im = obj.n_images-checkback:obj.n_images
				img_mat(:, :, im) = imread(obj.filename, im);
			end
			out = mean(img_mat, 3);
		end

		function checkForAffine2DRelaxed(obj)
			% Because the sequential rigid transformations has issues, this resolves it by forcing all transformations to
			% always be "isRigid" == true
			if exist('affine2d_relaxed.m', 'file') ~= 2
				error('This code will not work without the ''affine2d_relaxed'' class, make sure it''s on your path')
			end
		end

		function tform = registerFirstImage(obj, first_img)
			% Special registration for first image, because it needs to registered to the avg projection, which might require
			% additional tweaking here at the expense of speed...
			[optimizer, metric] = imregconfig('multimodal');
			optimizer.InitialRadius = optimizer.InitialRadius/100; % tiny tiny
			tform = imregtform(wiener2(first_img, [50, 50]), obj.template,'translation',optimizer, metric);
		end

		function processed_img = preProcessImage(obj, raw_img)
			processed_img = medfilt2(raw_img);
		end

		function writer = setDefaultTags(obj, writer, img_size, compression_flag)
			if nargin < 3 || isempty(img_size)
				img_size = size(obj.template);
			end

			if nargin < 3 || isempty(compression_flag)
				compression_flag = false;
			end

			if compression_flag
				compress = Tiff.Compression.LZW;
			else
				compress = Tiff.Compression.None;
			end

			% Setting default tags for the tiff writer
			writer.setTag('ImageWidth', img_size(2));
			writer.setTag('ImageLength', img_size(1));
			writer.setTag('Photometric', Tiff.Photometric.MinIsBlack);
			writer.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
			writer.setTag('Compression', compress);
			writer.setTag('BitsPerSample', 16);
			writer.setTag('SamplesPerPixel', 1);
		end

		function out = getCount(obj)
			% Determine what the current counter is by removing zeros off the end of tform_mat, careful with this..
			out = obj.idx;
		end

		function out = registrationReady(obj, tform)
			% Check if the registration is good to go
			[~, warn_id] = lastwarn();
			if ~isempty(warn_id) || obj.checkForJump(tform) || obj.driftingAway(tform) % if there's a warning (nonconvergence) or the jump is big, it fails
				out = false;
			else
				out = true;
			end

			lastwarn('',''); % reset warning
		end

		function out = driftingAway(obj, tform)
			idx = obj.getCount();
			if idx < 10
				out = false;
				return % no use in fitting if < 10 samples
			end
			lookback = 10;
			prev_tform = obj.tform_mat(:, :, max(1, idx - lookback) : max(1, idx - 1));
			for m = 1:size(prev_tform, 1)
				for n = 1:size(prev_tform, 2)
					f = polyfit(1:size(prev_tform, 3), squeeze(prev_tform(m, n, :)), 1);
					predicted_tform(m, n) = polyval(f, size(prev_tform, 3) + 1);
				end
			end

			max_ang = 20;
			max_pos = 100;

			if abs(acosd(tform.T(1, 1)) - acosd(predicted_tform(1, 1))) > max_ang || abs(predicted_tform(end, 1) - tform.T(end, 1)) > max_pos || abs(predicted_tform(end, 2) - tform.T(end, 2)) > max_pos

				out = true;
				disp("DRIFT!!")
			else
				out = false;
			end
		end

		function out = calculateMaxAng(obj)
			previous_angles = squeeze(acosd(obj.tform_mat(1, 1, max(1, obj.idx - 10) : obj.idx)));
			delta_angle = abs(diff(previous_angles));
			out = mean(delta_angle) + 5 * std(delta_angle);
		end

		function out = checkForJump(obj, tform, max_ang, max_pos)
			% Check for major jumps between frames, one of the criteria for having a good registration
			if nargin < 3 || isempty(max_angle)
				max_ang = obj.calculateMaxAng();
			end
			if nargin < 4 || isempty(max_pos)
				max_pos = 50;
			end
			[next_ang, next_pos_x, next_pos_y] = obj.getPreviousInfo();

			if abs(acosd(tform.T(1, 1)) - next_ang) > max(5, max_ang) || abs(next_pos_x - tform.T(end, 1)) > max_pos || abs(next_pos_y - tform.T(end, 2)) > max_pos
				out = true;
				disp('JUMP!!')
			else
				out = false;
			end
		end

		function [next_ang, next_pos_x, next_pos_y] = getPreviousInfo(obj)
			% Predict values for the next transformation matrix based on the previous, right now the predicted is just the
			% previous one...
			prev_tform = obj.tform_mat(:, :, max(1, obj.idx - 1));
			next_ang = acosd(prev_tform(1, 1));

			next_pos_x = prev_tform(end, 1);
			next_pos_y = prev_tform(end, 2);
		end

		function out = getBestGuess(obj)
			% Keep it simple, just based on the previous transform
			idx = obj.getCount();
			out = obj.tform_mat(:, :, idx-1);
		end

		function out = getAvgImage(obj, img)
			% Average image scaled by distance from the current image
			idx = obj.getCount() - 2 : obj.getCount() + 2;
			idx = max(min(idx, size(img, 3)), 1);
			img_mat = img(:, :, idx);
			% scale it
			scale_vector = [0.1, 0.25, 1, 0.25, 0.1];
			for j = 1:length(scale_vector)
				img_mat(:, :, j) = img_mat(:, :, j) * scale_vector(j);
			end
			out = mean(img_mat, 3);
		end

		function checkForRecordedAngle(obj)
			gt_answer = 'No';
			% 		gt_answer = questdlg('Did you record ground truth angle?', 'Recorded Angle', 'Yes', 'No', 'Yes');
			switch gt_answer
				case 'Yes'
					fprintf('Choose your matfile containing the encoder counts: \n')
					[fn, pn] = uigetfile('.mat');
					obj.recorded_angle = obj.processRecordedAngle(strcat([pn, '/', fn]));
					obj.registration_type = 'rigid';
				case 'No'
					obj.recorded_angle = [];
					obj.registration_type = 'rigid';
			end
		end

		function out = processRecordedAngle(obj, recording_path)
			raw_counts = importdata(recording_path);
			big_diameter = 64.8; %mm
			small_diameter = 11.83; %mm
			counts_per_rev = 128 * 2; % i think this is correct...
			resampled_counts = obj.equalizeFs(raw_counts);
			out = mod(resampled_counts, counts_per_rev * (big_diameter/small_diameter));
			% final angle for avg projection set to "0"
			out = rescale(out) * 2*pi;
		end

		function resampled_counts = equalizeFs(obj, raw_counts)
			raw_counts = raw_counts(find(diff(raw_counts ~= 0), 1, 'first'):end) ; %sync to start
			bins = discretize(1:length(raw_counts), obj.n_images-50); % accounting for the first 50 not moving, will add later
			resampled_counts = zeros(1, length(unique(bins)));
			for u = unique(bins)
				resampled_counts(u) = median(raw_counts(bins == u));
			end
			resampled_counts = cat(2, zeros(1, 50), resampled_counts);
		end

		function visualize(obj, raw_img, derotated_img, best_guess_img, is_jump, time_vec)
			% Basic visualization for seeing how its doing in realtime
			ct = obj.getCount();

			subplot(3, 3, 1)
			imagesc(best_guess_img)
			title(sprintf('frame #%d, unregistered', ct))
			axis image
			axis off

			subplot(3, 3, 2)
			imagesc(derotated_img)
			title('registered')
			axis image
			axis off

			subplot(3, 3, 4)
			imagesc(obj.template)
			title('running template')
			axis image
			axis off

			subplot(3, 3, 5)
			imagesc(obj.occupancy)
			title('running occupancy')
			axis image
			axis off

			% performance metrics
			subplot(3, 3, 3)
			plot(1:obj.idx, squeeze(obj.tform_mat(end, 1, 1:obj.idx)));
			ylabel('x translation')

			subplot(3, 3, 6)
			plot(1:obj.idx, squeeze(obj.tform_mat(end, 2, 1:obj.idx)));
			ylabel('y translation')

			subplot(3, 3, 9)
			plot(1:obj.idx, squeeze(acosd(obj.tform_mat(1, 1, 1:obj.idx))));
			ylabel('rotation')

			subplot(3, 3, 7)
			scatter(1:obj.idx, is_jump(1:obj.idx), 'filled');
			hold on
			plot(1:obj.idx, is_jump(1:obj.idx));
			hold off
			ylim([0, 1.2])
			ylabel('jump')

			subplot(3, 3, 8)
			plot(1:obj.idx, time_vec(1:obj.idx));
			ylabel('time to register')

			pause(0.01) % necessary to get it to show upc
		end


		function LRout = getLargestRectangle(obj, image)

			%% from https://www.mathworks.com/matlabcentral/fileexchange/71491-largest-inscribed-rectangle-square-or-circle

			%Function to find the largest inscribed rectangle in an arbitrary shape with multiple holes.
			%INPUT:
			% image:        Image, RGB, grey or BW. By preference BW.
			% RotationStep: Default: 5°
			%               Rotation Step in degrees. In order to find tilted rectangles, the image is rotated.
			%               Range: 0 < RotationStep <= 90.
			%               If RotationStep>(LastAngle-FirstAngle) then the image is rotated once.
			% iterate:      Default: 1 (with iteration)
			%               If 0 then no iteration, if 1 then iteration.
			%               No iteration if (FirstAngle == LastAngle)
			%               It might not find the largest rectangle in the range determined by the RotationStep,
			%               it might lock on any local maximum within a rotation step.
			%FirstAngle:    Default: 0°
			%               First angle:  0 <= FirstAngle < 90
			%LastAngle:     Default: 89.9999°
			%               Last angle: FirstAngle <= LastAngle <= 90
			%Graphic:       Default: 1 (Plot graphic), 0: no graphic
			%OUTPUT LRout:
			% 1st row: Area of the largest rectangle in px, Rotation angle in degrees counterclockwise
			% 2nd row: x,y of top corner of largest rectangle
			% 3rd row: x,y of right corner of largest rectangle
			% 4th row: x,y of bottom corner of largest rectangle
			% 5th row:x,y of left corner of largest rectangle
			%EXAMPLES:
			% LRout=LargestRectangle(myImage)% Run with default values
			% LRout=LargestRectangle(myImage,0.1,0)% Run with small rotation steps, no iteration
			% LRout=LargestRectangle(myImage,0,0,0,0)% For axis parallel rectangle
			% Example for rotation step=5, no iteration, starting at 10°, ending at 20°, no graphic:
			%    LRout=LargestRectangle(myImage,5,0,10,20,0)
			%REMARK:
			% Any hole (black pixel), even only one pixel wide,
			%  will not be inside the largest rectangle
			%Revisions:
			%19.05.20: removed bug when left or top border of input image has a white pixel
			inArg=[5,1,0,89.9999,1];%Default values, you may want to change it
			RotationStep=inArg(1);
			iterate=inArg(2);
			StartAngle=inArg(3);
			EndAngle=inArg(4);
			numRotSteps=floor((EndAngle-StartAngle)/RotationStep)+1;
			if EndAngle==StartAngle
				iterate=0;% NO iteration
				RotationStep=1;% dummy value >0
				numRotSteps=1;
			end
			if iterate && floor((EndAngle-StartAngle)/RotationStep)+1<10
				% make rotatationstep smaller and keep user selected values
				RotationStep=RotationStep/ceil(10*RotationStep/(EndAngle-StartAngle));
				numRotSteps=floor((EndAngle-StartAngle)/RotationStep)+1;
			end
			if ~islogical(image)
				image=im2bw(image);
			end
			ImBWc=image;
			%add 1px black border arround ImBWc:
			%You gain approx. 1% area if you remove 'logical', but 3 times slower
			ImBWc=logical([zeros(1,size(ImBWc,2)+2);...%add top border
				zeros(size(ImBWc,1),1),ImBWc,zeros(size(ImBWc,1),1);...%add left+right border  + image
				zeros(1,size(ImBWc,2)+2)]);%add bottom border
		[rows, columns] = find(ImBWc);
		if isempty(rows)
			%totaly black image
			disp ('Totaly black image! No rectangle found.');
			LRout=[0,0;1,1;1,1;1,1;1,1];
			return;
		end
		%determine smallest black border arround image
		leftCrop = min(columns)-1;
		topCrop = min(rows)-1;
		rightCrop = max(columns)+1;
		bottomCrop = max(rows)+1;
		%crop image:
		ImBWc=ImBWc(topCrop:bottomCrop,leftCrop:rightCrop);
		%correct origin for later use (1px was added as border):
		leftCrop = leftCrop-1;
		topCrop = topCrop-1;
		LR=[0,0,0,0,0,0,0,0,0];%area,xmin,xmax,ymin,ymax, angle, sizeImR(1), sizeImR(2),AreaNow
		%Main loop, rotate image get border, get largest rectangle
		AN=zeros(numRotSteps,2);% Angle, LR, used for iterate only
		ANi=0;
		for angle=StartAngle:RotationStep:EndAngle
			[border,ColLims,RowLims,sizeImR]=GetContour(ImBWc,angle);
			if border(1)==0
				%Only 1 px areas
				[border,ColLims,RowLims,sizeImR]=GetContour(ImBWc,0);
				LR=getLR(LR,border,ColLims,RowLims);
				LR(6:8)=[0,sizeImR];
				iterate=0;
				break;%Exit for loop
			end
			areaOld=LR(1);
			LR=getLR(LR,border,ColLims,RowLims);
			ANi=ANi+1;
			AN(ANi,:)=[angle,LR(9)];
			if LR(1)>areaOld
				LR(6:8)=[angle,sizeImR];
			end
		end
		if iterate
			%iterate
			MinRotationStep=min(0.01,RotationStep/10);
			AN=[AN;AN(end,1)+RotationStep,0];
			ANi=ANi+1;
			LimitDiff=0.05*0.85^-log2(RotationStep/15);
			angles=[AN(1,1),AN(end,1)];
			Transit=true;%for first iteration at <0.5°
			while RotationStep>MinRotationStep
				RotationStep2=RotationStep/2;
				numelAN=size(AN,1);
				LimitDiff=LimitDiff*0.85;%decreasing last number results in narrower search range
				Limit=LR(1)*(1-LimitDiff);
				iLold=find(AN(:,1)>angles(1)-RotationStep2,1,'first');
				iRold=find(AN(:,1)>angles(end)-RotationStep2,1,'first');
				iLRmaxM=find(AN(iLold:iRold,2)==LR(1))+iLold-1;%multiple iLRmax
				%check if the LRmax are at around 0° AND around 90° (border problem!)
				if AN(iLRmaxM(1,1))<15 && AN(iLRmaxM(end,1))>85
					%border problem!
					iLRmaxML=iLRmaxM(AN(iLRmaxM(:),1)<45);
					iLRmaxMR=iLRmaxM(iLRmaxML+1:end);
					%take the side with the most max LR
					if numel(iLRmaxML)>=numel(iLRmaxMR)
						iLRmaxM=iLRmaxML;
					else
						iLRmaxM=iLRmaxMR;
					end
				end
				iLRmax=iLRmaxM(floor((numel(iLRmaxM)+1)/2));%take position of midle LRmax
				if Transit && RotationStep<0.5
					%correct iLold and iRold
					if inArg(1)<0.5
						%no previous iteration
						iLstart=max(1,iLRmaxM(1)-20);
						iRend=min(numelAN,iLRmaxM(end)+20);
						iLold=find(AN(iLstart:iLRmaxM(1),2)<Limit,1,'last')+iLstart-1;%Left index
						iRold=find(AN(iLRmaxM(end):iRend,2)<Limit,1,'first')+iLRmaxM(end)-1;%Right index
						if isempty(iLold);iLold=iLstart;end
						if iLold<1; iLold=1;end
						if isempty(iRold);iRold=iRend;end
						if iRold>iRend;iRold=iRend;end
					else
						iLold=find(AN(:,1)==angles(find(angles>AN(iLRmaxM(1),1)-5,1,'first')));
						iRold=find(AN(:,1)==angles(find(angles<AN(iLRmaxM(end),1)+5,1,'last')));
					end
					Transit=false;
				end
				iLstart=max(iLRmax-25,iLold);
				iRend=min(iLRmax+25,iRold);
				if RotationStep<0.5 %last number by experiment
					%include limited range arround max peak and second largest peak
					ANs=sort(AN(iLstart:iRend,2),'descend');
					ANs=ANs(ANs<LR(1));%exclude (multiple) LR(1)
					if isempty(ANs)
						LimitE=LR(1)-1;
					else
						LimitE=ANs(1)-1;%second largest-1
					end
					iLE=find(AN(iLstart:iLRmax-1,2)>LimitE,1,'first')+iLstart-2;
					if iLE<1;iLE=iLstart;end
					iRE=find(AN(iLRmax+1:iRend,2)>LimitE,1,'last')+iLRmax+1;%Right index
					if iRE>numelAN;iRE=numelAN;end
					iL=find(AN(iLstart:iLRmax,2)<Limit,1,'last')+iLstart-1;%Left index
					if isempty(iL) || iL<1;iL=iLstart;end
					iR=find(AN(iLRmax:iRend,2)<Limit,1,'first')+iLRmax-1;%Right index
					if isempty(iR) ||iR>iRend;iR=iRend;end
					if iLE<iL
						iL=iLE;
					end
					if iRE>iR
						iR=iRE;
					end
					maxRotIt=min(22,round(6.5+1.7./sqrt(RotationStep)));%9 at 0.4°, 22 at 0.1°
					maxRot=maxRotIt;
					while iR-iL>maxRot
						%delete smaller LR
						if AN(iL,2)<= AN(iR,2) && iL<iLRmaxM(1)
							iL=iL+1;
						elseif iR>iLRmaxM(end)
							iR=iR-1;
						else
							iL=max(iL,iLRmax-12);
							iR=min(iR,iLRmax+12);
							break;
						end
						if iL==iLRmaxM(1) ||  iR==iLRmaxM(end)
							%propably unsymetric peak position
							maxRot=24;
						end
					end
					if iR-iL<4 || (iL>=iLRmaxM(1) && iR<=iLRmaxM(end))
						%only a few rot steps, increase by 2
						angles=AN(iL,1)-RotationStep2:RotationStep:AN(iR,1)+1.1*RotationStep2;
					elseif iL>=iLRmaxM(1)
						angles=AN(iL,1)-RotationStep2:RotationStep:AN(iR,1);
					elseif iR<=iLRmaxM(end)
						angles=AN(iL,1)+RotationStep2:RotationStep:AN(iR,1)+1.1*RotationStep2;
					else
						%regular case
						angles=AN(iL,1)+RotationStep2:RotationStep:AN(iR,1);
					end
				else
					%RotationStep>=0.5, select largest peaks
					iL=find(AN(iLstart:iLRmax,2)>Limit,1,'first')+iLstart-2;%Left index
					iR=find(AN(iLRmax:iRend,2)>Limit,1,'last')+iLRmax;%Right index
					if isempty(iL);iL=iLRmax-2;end
					if iL<1; iL=1;end
					if isempty(iR);iR=iLRmax+2;end
					if iR>iRend;iR=iRend;end
					LimitLow=2*Limit-LR(1);
					angles=AN(iL:iR,:);
					angles=angles(angles(:,2)>LimitLow,:);%exclude small LRs
					maxRot=7;
					if length(angles)>maxRot
						angles=sortrows(angles,-2);
						angles=angles(1:maxRot,1);
					else
						angles=angles(:,1);
					end
					angles=sort([angles-RotationStep2;angles+RotationStep2]);
					%exclude duplicate angles:
					angles=[angles(abs(angles(1:end-1)-angles(2:end))>RotationStep2);angles(end)];
				end% if RotationStep<0.5
				for j=1:numel(angles)
					angle=angles(j);
					[border,ColLims,RowLims,sizeImR]=GetContour(ImBWc,angle);
					areaOld=LR(1);
					LR=getLR(LR,border,ColLims,RowLims);
					ANi=ANi+1;
					AN(ANi,:)=[angle,LR(9)];
					if LR(1)>areaOld
						LR(6:8)=[angle,sizeImR];
					end
				end
				AN=sortrows(AN);
				RotationStep=RotationStep2;
			end
		end%End iterate

		%Prepaire LRout
		ImRcenter=([LR(7),LR(8)]+1)/2;
		ImBWcCenter=(size(ImBWc)+1)/2;
		xy=[[LR(2),LR(3),LR(3),LR(2)]-ImRcenter(2);[LR(4),LR(4),LR(5),LR(5)]-ImRcenter(1)];
		phi=LR(6)/180*pi;
		cosPhi=cos(phi);
		sinPhi=sin(phi);
		RM= [cosPhi -sinPhi; sinPhi cosPhi];
		xy=RM*xy;
		x=xy(1,:)+ImBWcCenter(2)+leftCrop-1;
		y=xy(2,:)+ImBWcCenter(1)+topCrop-1;
		LRout=[LR(1),LR(6);x(1),y(1);x(2),y(2);x(3),y(3);x(4),y(4)];

			function [borderD2,ColLims,RowLims,sizeImR]=GetContour(ImBWc,angle)
				ImBWc=imrotate(ImBWc,angle,'bilinear');
				sizeImR=size(ImBWc);
				sizeImR1=sizeImR(1);
				boundary = bwperim(ImBWc);
				[row,col] = find(boundary);
				borderLinIdxL=(col-2)*sizeImR1+row;% indices for left side of border point
				borderLinIdxR=borderLinIdxL+2*sizeImR1;% indices for right of border points
				borderLinIdxT=borderLinIdxL+sizeImR1-1;% indices for above of border points
				borderLinIdxB=borderLinIdxT+2;% indices for below of border points
				%top-bottom limits:
				borderLinIdx=[borderLinIdxT;borderLinIdxB];
				indx=unique(borderLinIdx(ImBWc(borderLinIdx)==0));
				c=floor((indx-1)/sizeImR1)+1;
				borderD=[c,indx-(c-1)*sizeImR1];%x,y
				[~,mm]=mode(borderD(:,1));
				numRows=borderD(end,1);
				indx=borderD(:,1);
				indx3=indx;
				ColLims=NaN(numRows,mm);
				for i=2:numel(indx)
					if indx(i)==indx(i-1)
						indx3(i)=indx3(i-1)+numRows;
					else
						indx3(i)=indx(i);
					end
				end
				ColLims(indx3)=borderD(:,2);
				%Left-right limits
				borderLinIdx=[borderLinIdxL;borderLinIdxR];
				indx=unique(borderLinIdx(ImBWc(borderLinIdx)==0));
				c=floor((indx-1)/sizeImR1)+1;
				borderD=sortrows([c,indx-(c-1)*sizeImR1],2);%x,y
				[~,mm]=mode(borderD(:,2));
				numRows=borderD(end,2);
				indx=borderD(:,2);
				indx3=indx;
				RowLims=NaN(numRows,mm);
				for i=2:numel(indx)
					if indx(i)==indx(i-1)
						indx3(i)=indx3(i-1)+numRows;
					else
						indx3(i)=indx(i);
					end
				end
				RowLims(indx3)=borderD(:,1);
				%find top border px only
				indxT=borderLinIdxT(ImBWc(borderLinIdxT)==0);
				%find left border px only
				indxL=borderLinIdxL(ImBWc(borderLinIdxL)==0);
				if numel(indxT)<numel(indxL)
					%Use top border only for evaluation
					c=floor((indxT-1)/sizeImR1)+1;
					borderD2=[c,indxT-(c-1)*sizeImR1+1];
				else
					c=floor((indxL-1)/sizeImR1)+1;
					borderD2=[c+1,indxL-(c-1)*sizeImR1];
				end
			end

			function LR=getLR(LR,border,ColLims,RowLims)
				%get largest rectangle
				AreaNow=0;%only for iteration
				for i=1:size(border,1)
					%loop border points
					xNow=border(i,1);
					yNow=border(i,2);
					Indx=1;while ColLims(xNow,Indx)<yNow;Indx=Indx+1;end%much faster than find.m
					yLimit1=ColLims(xNow,Indx-1)+1;yLimit2=ColLims(xNow,Indx)-1;%ylimits @xNow
					Indx=1;while RowLims(yNow,Indx)<xNow;Indx=Indx+1;end
					xLimit1=RowLims(yNow,Indx-1)+1;xLimit2=RowLims(yNow,Indx)-1;%xlimits @yNow
					heightMax=yLimit2-yLimit1+1;
					widthMax=xLimit2-xLimit1+1;
					if widthMax*heightMax>LR(1)
						if (yNow-yLimit1+1)*(yLimit2-yNow+1)<(xNow-xLimit1+1)*(xLimit2-xNow+1)
							%Max Number of Operation y direction smaller than for x direction
							if yNow==yLimit1
								%Border pixel at top  ===========================================
								xLimitTempT1=xLimit1;
								xLimitTempT2=xLimit2;
								xLimitTempB1=xLimit1;
								xLimitTempB2=xLimit2;
								Indx=1;while RowLims(yNow,Indx)<xNow;Indx=Indx+1;end
								xT2=RowLims(yNow,Indx)-1;xT1=RowLims(yNow,Indx-1)+1;%faster in one line
								if xT1>xLimitTempT1;xLimitTempT1=xT1;end
								if xT2<xLimitTempT2;xLimitTempT2=xT2;end
								if (xLimitTempT2-xLimitTempT1+1)*heightMax<LR(1)
									break;
								end
								for yB=yNow:yLimit2
									Indx=1;while RowLims(yB,Indx)<xNow;Indx=Indx+1;end
									xB2=RowLims(yB,Indx)-1;xB1=RowLims(yB,Indx-1)+1;%faster in one line
									if xB1>xLimitTempB1;xLimitTempB1=xB1;end
									if xB2<xLimitTempB2;xLimitTempB2=xB2;end
									if xLimitTempB1>xLimitTempT1;xL=xLimitTempB1;else xL=xLimitTempT1;end
									if xLimitTempB2>xLimitTempT2;xR=xLimitTempT2;else xR=xLimitTempB2;end
									if (xR-xL+1)*heightMax<LR(1)
										break;
									end
									area=(xR-xL+1)*(yB-yNow+1);
									if area>AreaNow
										AreaNow=area;
										if area>LR(1)
											LR(1:5)=[area,xL,xR,yNow,yB];
										end
									end
								end%for yB=yNow:yLimit2
							elseif yNow==yLimit2
								%Border pixel at bottom  ========================================
								xLimitTempB1=xLimit1;
								xLimitTempB2=xLimit2;
								xLimitTempT1=xLimit1;
								xLimitTempT2=xLimit2;
								Indx=1;while RowLims(yNow,Indx)<xNow;Indx=Indx+1;end
								xB2=RowLims(yNow,Indx)-1;
								xB1=RowLims(yNow,Indx-1)+1;
								if xB1>xLimitTempB1;xLimitTempB1=xB1;end
								if xB2<xLimitTempB2;xLimitTempB2=xB2;end
								if (xLimitTempB2-xLimitTempB1+1)*heightMax>LR(1)
									for yT=yNow:-1:yLimit1
										Indx=1;while RowLims(yT,Indx)<xNow;Indx=Indx+1;end
										xT2=RowLims(yT,Indx)-1;
										xT1=RowLims(yT,Indx-1)+1;
										if xT1>xLimitTempT1;xLimitTempT1=xT1;end
										if xT2<xLimitTempT2;xLimitTempT2=xT2;end
										if xLimitTempB1>xLimitTempT1;xL=xLimitTempB1;else xL=xLimitTempT1;end
										if xLimitTempB2<xLimitTempT2;xR=xLimitTempB2;else xR=xLimitTempT2;end
										if (xR-xL+1)*heightMax<LR(1)
											break;
										end
										area=(xR-xL+1)*(yNow-yT+1);
										if area>AreaNow
											AreaNow=area;
											if area>LR(1)
												LR(1:5)=[area,xL,xR,yT,yNow];
											end
										end
									end%
								end%if (xLimitTempB2-xLimitTempB1+1)*heightMax>LR(1)
							else
								%Border pixel between top and bottom ============================
								xLimitTempT1=xLimit1;
								xLimitTempT2=xLimit2;
								%prepare data for inner loop
								xB1a=zeros(yLimit2-yNow+1,1);
								xB2a=xB1a;
								for yB=yNow:yLimit2
									Indx=1;while RowLims(yB,Indx)<xNow;Indx=Indx+1;end
									xB2a(yB-yNow+1)=RowLims(yB,Indx)-1;xB1a(yB-yNow+1)=RowLims(yB,Indx-1)+1;%faster in one line
								end
								for yT=yNow:-1:yLimit1
									xLimitTempB1=xLimit1;
									xLimitTempB2=xLimit2;
									Indx=1;while RowLims(yT,Indx)<xNow;Indx=Indx+1;end
									xT2=RowLims(yT,Indx)-1;xT1=RowLims(yT,Indx-1)+1;%faster in one line
									if xT1>xLimitTempT1;xLimitTempT1=xT1;end
									if xT2<xLimitTempT2;xLimitTempT2=xT2;end
									if (xLimitTempT2-xLimitTempT1+1)*heightMax<LR(1)
										break;
									end
									for yB=yNow:yLimit2
										if xB1a(yB-yNow+1)>xLimitTempB1;xLimitTempB1=xB1a(yB-yNow+1);end
										if xB2a(yB-yNow+1)<xLimitTempB2;xLimitTempB2=xB2a(yB-yNow+1);end
										if xLimitTempB1>xLimitTempT1;xL=xLimitTempB1;else xL=xLimitTempT1;end
										if xLimitTempB2>xLimitTempT2;xR=xLimitTempT2;else xR=xLimitTempB2;end
										if (xR-xL+1)*(yLimit2-yT+1)<LR(1)
											break;
										end
										area=(xR-xL+1)*(yB-yT+1);
										if area>AreaNow
											AreaNow=area;
											if area>LR(1)
												LR(1:5)=[area,xL,xR,yT,yB];
											end
										end
									end%for yB=yNow:yLimit2
								end%for yT=yNow:-1:yLimit1
							end%
						else%
							%Max Number of Operation x direction smaller than for y direction ==================
							if xNow==xLimit1
								%Border pixel on left side ======================================
								yLimitTempL1=yLimit1;%ylimit Left
								yLimitTempL2=yLimit2;%ylimit Left
								yLimitTempR1=yLimit1;
								yLimitTempR2=yLimit2;
								Indx=1;while ColLims(xNow,Indx)<yNow;Indx=Indx+1;end
								yL2=ColLims(xNow,Indx)-1;yL1=ColLims(xNow,Indx-1)+1;%faster in one line
								if yL1>yLimitTempL1;yLimitTempL1=yL1;end
								if yL2<yLimitTempL2;yLimitTempL2=yL2;end
								if widthMax*(yLimitTempL2-yLimitTempL1+1)>LR(1)
									for xR=xNow:xLimit2
										Indx=1;while ColLims(xR,Indx)<yNow;Indx=Indx+1;end
										yR2=ColLims(xR,Indx)-1;yR1=ColLims(xR,Indx-1)+1;%faster in one line
										if yR1>yLimitTempR1;yLimitTempR1=yR1;end%much faster than yR1=max([yLimitTempR1,yR1]);
										if yR2<yLimitTempR2;yLimitTempR2=yR2;end
										if yLimitTempR1>yLimitTempL1;yT=yLimitTempR1;else yT=yLimitTempL1;end
										if yLimitTempR2>yLimitTempL2;yB=yLimitTempL2;else yB=yLimitTempR2;end
										if widthMax*(yB-yT+1)<LR(1)
											break;
										end
										area=(xR-xNow+1)*(yB-yT+1);
										if area>AreaNow
											AreaNow=area;
											if area>LR(1)
												LR(1:5)=[area,xNow,xR,yT,yB];
											end
										end
									end%
								end%if widthMax*(yLimitTempL2-yLimitTempL1+1)>LR(1)
							elseif xNow==xLimit2%if xNow==xLimit1
								%Border pixel on right side =====================================
								yLimitTempR1=yLimit1;
								yLimitTempR2=yLimit2;
								Indx=1;while ColLims(xNow,Indx)<yNow;Indx=Indx+1;end
								yR2=ColLims(xNow,Indx)-1;yR1=ColLims(xNow,Indx-1)+1;%faster in one line
								if yR1>yLimitTempR1;yLimitTempR1=yR1;end
								if yR2<yLimitTempR2;yLimitTempR2=yR2;end
								if widthMax*(yLimitTempR2-yLimitTempR1+1)<LR(1)
									break;
								end
								yLimitTempL1=yLimit1;
								yLimitTempL2=yLimit2;
								for xL=xNow:-1:xLimit1
									Indx=1;while ColLims(xL,Indx)<yNow;Indx=Indx+1;end
									yL2=ColLims(xL,Indx)-1;yL1=ColLims(xL,Indx-1)+1;%faster in one line
									if yL1>yLimitTempL1;yLimitTempL1=yL1;end
									if yL2<yLimitTempL2;yLimitTempL2=yL2;end
									if yLimitTempR1>yLimitTempL1;yT=yLimitTempR1;else yT=yLimitTempL1;end
									if yLimitTempR2>yLimitTempL2; yB=yLimitTempL2;else yB=yLimitTempR2;end
									if widthMax*(yB-yT+1)<LR(1)
										break;
									end
									area=(xNow-xL+1)*(yB-yT+1);
									if area>AreaNow
										AreaNow=area;
										if area>LR(1)
											LR(1:5)=[area,xL,xNow,yT,yB];
										end
									end
								end%for xL=xNow:-1:xLimit1
							else
								%Border pixel between left and right side =======================
								yLimitTempR1=yLimit1;
								yLimitTempR2=yLimit2;
								%prepare data for inner loop
								yL1a=zeros(xNow-xLimit1+1,1);
								yL2a=yL1a;
								for xL=xNow:-1:xLimit1
									Indx=1;while ColLims(xL,Indx)<yNow;Indx=Indx+1;end
									yL2a(xL-xLimit1+1)=ColLims(xL,Indx)-1;yL1a(xL-xLimit1+1)=ColLims(xL,Indx-1)+1;
								end
								for xR=xNow:xLimit2
									Indx=1;while ColLims(xR,Indx)<yNow;Indx=Indx+1;end
									yR2=ColLims(xR,Indx)-1;yR1=ColLims(xR,Indx-1)+1;%faster in one line

									if yR1>yLimitTempR1;yLimitTempR1=yR1;end
									if yR2<yLimitTempR2;yLimitTempR2=yR2;end
									if widthMax*(yLimitTempR2-yLimitTempR1+1)<LR(1)
										break;
									end
									yLimitTempL1=yLimit1;
									yLimitTempL2=yLimit2;
									for xL=xNow:-1:xLimit1
										if yL1a(xL-xLimit1+1)>yLimitTempL1;yLimitTempL1=yL1a(xL-xLimit1+1);end
										if yL2a(xL-xLimit1+1)<yLimitTempL2;yLimitTempL2=yL2a(xL-xLimit1+1);end
										if yLimitTempR1>yLimitTempL1;yT=yLimitTempR1;else yT=yLimitTempL1;end
										if yLimitTempR2>yLimitTempL2; yB=yLimitTempL2;else yB=yLimitTempR2;end
										if (xR-xLimit1+1)*(yB-yT+1)<LR(1)
											break;
										end
										area=(xR-xL+1)*(yB-yT+1);
										if area>AreaNow
											AreaNow=area;
											if area>LR(1)
												LR(1:5)=[area,xL,xR,yT,yB];
											end
										end
									end%for xL=xNow:-1:xLimit1
								end%for xR=xNow:xLimit2
							end%if xNow==xLimit1
						end%if xLimit2-xLimit1>yLimit2-yLimit1
					end%if widthMax*heightMax>LR(1)
				end%for i=1:size(border,1)
				LR(9)=AreaNow;
			end
		end
	end

