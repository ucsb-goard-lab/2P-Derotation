# 2P Derotation Pipeline
This code is to de-rotate calcium imaging timeseries from slowly rotating recordings. 
There is no provision for correcting for the nonrigid warping of the image due to the movement of the scan lines.
As a result, this pipeline is best suited for recordings in which the maximum instantaneous rotation of the image is ~18&deg/s.

## Installation
Clone this repository and add it to your path.
Everything is pretty self contained, so not much else is needed. 

## Usage
`full_pipeline_master.m` is a script that contains all the relevant information to process your own recordings.
It starts from .TIF timeseries data (either multi-page or many single page tifs), and goes through the entire procedure for 1) derotating rotating timeseries ("head_rotation") and 2) aligning derotated timeseries to a second, non-rotated recording ("platform_rotation").

## Dependencies
