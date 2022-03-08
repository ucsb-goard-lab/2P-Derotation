%affine2d 2-D Affine Geometric Transformation
%
%   An affine2d object encapsulates a 2-D affine geometric transformation. 
%
%   affine2d properties:
%      T - 3x3 matrix representing forward affine transformation
%      Dimensionality - Dimensionality of geometric transformation
%
%   affine2d methods:
%      affine2d - Construct affine2d object
%      invert - Invert geometric transformation
%      isTranslation - Determine if transformation is pure translation special case
%      isRigid - Determine if transformation is rigid transformation special case
%      isSimilarity - Determine if transformation is similarity transformation special case
%      outputLimits - Find output spatial limits given input spatial limits
%      transformPointsForward - Apply forward 2-D geometric transformation to points
%      transformPointsInverse - Apply inverse 2-D geometric transformation to points
%
%   Example 1
%   ---------
%   % Construct an affine2d object that defines a rotation of 10 degrees
%   % counter-clockwise.
%   theta = 10;
%   tform = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1]);
%
%   % Apply forward geometric transformation to an input (U,V) point (5,10)
%   [X,Y] = transformPointsForward(tform,5,10)
%
%   % Apply inverse geometric transformation to output (X,Y) point from
%   % previous step. We recover the point we started with from
%   % the inverse transformation.
%   [U,V] = transformPointsInverse(tform,X,Y)
%
%   Example 2
%   ---------
%   % Apply 10 degree counter-clockwise rotation to an image using the function imwarp
%   A = imread('pout.tif');
%   theta = 10;
%   tform = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1]);
%   outputImage = imwarp(A,tform);
%   figure, imshow(outputImage);
%
%   See also AFFINE3D, PROJECTIVE2D, RIGID2D, IMWARP, GEOMETRICTRANSFORM2D, GEOMETRICTRANSFORM3D

% Copyright 2012-2020 The MathWorks, Inc.

%#ok<*EMCA>

classdef affine2d_relaxed < affine2d
        
    properties
    end
    
    
    methods
        
        function self = affine2d_relaxed(A)
            %affine2d Construct affine2d object
            %
            %   tform = affine2d() constructs an affine2d object with default
            %   property settings that correspond to the identity
            %   transformation.
            %
            %   tform = affine2d(A) constructs an affine2d object given an
            %   input 3x3 matrix A that specifies a valid 3x3 affine
            %   transformation matrix. A must be of the form:
            %
            %    A = [a b 0;
            %         c d 0;
            %         e f 1];
            self = self@affine2d(A);           
        end
        

        function TF = isRigid(self)
            % overloads the superclass isRigid

            %isRigid Determine if transformation is rigid transformation
            %
            %   TF = isRigid(tform) determines whether or not affine
            %   transformation is a rigid transformation. TF is a scalar
            %   boolean that is true when tform is a rigid transformation. The
            %   tform is a rigid transformation when tform.T defines only
            %   rotation and translation.

            coder.inline('always');
            coder.internal.prefer_const(self);
            
            % TF = isSimilarity(self) && abs(det(self.T)-1) < 10.^10 * eps(class(self.T));
            TF = true;

        end
    end

end
