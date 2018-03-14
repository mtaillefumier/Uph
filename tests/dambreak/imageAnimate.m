% make_avi_movie_example1.m
%
%  Discussion:
%
%    This program loads a series of sequentially indexed png files to 
%    create an Audio Video Interleaved (AVI) movie that can be played
%    independently of Matlab, using, for example, the XINE player
%    for Linux. We capture NUMFRAMES frames for the x-y plot of breaking
%    waves with phase changing from 0 to 2 pi.  
%
%    Before running this program close all figures and don't interfere 
%    with the generated figures until the recording process is done.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 January 2015
%
%  Author:
%
%    Chola Kalale
%
% Create video with image sequence

clear all
clc

% set current folder
cd('/home/kchola/project/results');

% obtain all the png files in the current folder
output = dir('anime*.png');
name = {output.name};
% sort the files sequentially

[~,index] = sort_nat(name);
output = output(index);

% find the total number of png files in the current folder 
numfiles = size(output,1); % alternatively; numfiles = length(output);


%  Set the total number of frames that we will generate.
  num_frames_per_second = 10;


% create avi object
  aviobj = VideoWriter( 'breakingWave.avi');  % to write video file
  aviobj.FrameRate = num_frames_per_second;
  % define video quality [0 to 100]
  aviobj.Quality = 80;
  
  % open the file sphMovie.avi
  open(aviobj);
  for k = 1:numfiles
      % read images in the current folder one by one
      f{k} = imread(output(k).name);
      % the sizes of the images are made the same
      ResizeImg = imresize(f{k}, [200 3000]);
      
      % convert image to mvoie frame
      frame = im2frame(ResizeImg);
      
      % each image is copied num_frames_per_second times
      for j = 1 : num_frames_per_second
          writeVideo(aviobj, frame);
      end
  end
 
 % close the file 'sphMovie.avi'
 close(aviobj);
      
% the fps is 5. since each image is replicated 5 times
% one image will be seen per second. Therefore, our video
% will run for 179 seconds with 179*5 = 895 images
  

%  Tell MATLAB we have completed the movie.
%
  hold off
  
  fprintf ( 1, '\n' );
  fprintf ( 1, 'MAKE_AVI_MOVIE_EXAMPLE1\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '  The movie file "sphMovie.avi" has been created.\n' );
