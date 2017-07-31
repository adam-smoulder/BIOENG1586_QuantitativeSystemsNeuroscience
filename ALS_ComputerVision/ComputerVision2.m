% BIOENG1586, Computer Vision Homework
% Adam Smoulder

% use "load(rose)" manually by double clicking before running script!

%% 2a) 
OR = 0;
SF = 0.01;

% 2b) 
[X,Y] = meshgrid(-20:1:20);

% 2c)
std_x = 7;
std_y = 17;

X = X*cos(OR)+Y*sin(OR);
Y = Y*cos(OR)+X*sin(OR);

A = sin(2*pi*SF*X)/(2*pi*std_x*std_y);
gabor = A.*exp(-((X.^2/(2*std_x^2))+(Y.^2/(2*std_y^2))));

% 2d)
figure(6);
imagesc(gabor)
colormap gray
colorbar

% 2e)
figure(7);
subplot(2,2,1);
imagesc(rose)
colormap gray
%xlabel(strcat('Gabor, OR = ',num2str(OR),' SF = ',num2str(SF)))
subplot(2,2,4);
imagesc(conv2(double(rose), gabor, 'valid'));
colormap gray

% 2f)
OR = 0;
SF = 0.1;

X = X*cos(OR)+Y*sin(OR);
Y = Y*cos(OR)+X*sin(OR);
A = sin(2*pi*SF*X)/(2*pi*std_x*std_y);
gabor = A.*exp(-((X.^2/(2*std_x^2))+(Y.^2/(2*std_y^2))));

figure(7);
colormap gray
subplot(2,2,2);
%xlabel(strcat('Gabor, OR = ',num2str(OR),' SF = ',num2str(SF)))
imagesc(conv2(double(rose), gabor, 'valid'));
colormap gray

OR = pi/2;
SF = 0.01;

X = X*cos(OR)+Y*sin(OR);
Y = Y*cos(OR)+X*sin(OR);
A = sin(2*pi*SF*X)/(2*pi*std_x*std_y);
gabor = A.*exp(-((X.^2/(2*std_x^2))+(Y.^2/(2*std_y^2))));

colormap gray
subplot(2,2,3);
%xlabel(strcat('Gabor, OR = ',num2str(OR),' SF = ',num2str(SF)))
imagesc(conv2(double(rose), gabor, 'valid'));
colormap gray