% BIOENG1586, Computer Vision Homework
% Adam Smoulder

%% 3a)
load('hgrid.mat');
figure(8);
%subplot(2,1,1);
imshow(hgrid);

% 3b)
%subplot(2,1,2);
%recField = retina(.31, .15, hgrid);
figure(9)
recField = retina(.95, .01, hgrid);
