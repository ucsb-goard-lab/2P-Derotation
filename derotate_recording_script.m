filename_to_tiff_stack = '~/your_tiff_stack' % fill this in
reference_img = fill_this_in % should be a matrix of x by y pixels that's a reference
% image you want everything to register against, you can input the first frame, but it's
% nicer to have an average of more frames to have a cleaner reference
dr = Derotater(filename_to_tiff_stack, reference_img); % instantiate rotator object object

dr.dirtyRegister(); % first pass at registration
dr.cleanRegister(); % adding some smoothing and fine grained adjustments
dr.cleanUp(); 
dr.save();