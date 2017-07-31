% BIOENG1586, Computer Vision Homework
% Adam Smoulder

%% 1a) Create the visual input for the Mach Band Illusion.
machBand = zeros(64, 128);
currentBrightness = 10;
machBand(:,1:32) = currentBrightness*ones(64,32);

for i=33:1:128
    if currentBrightness < 75
        currentBrightness = currentBrightness+1;
    end
    machBand(:,i) = currentBrightness*ones(64,1);
end

figure(1);
hold on;
title('Mach Band')
axis([1,128,1,64]);
imagesc(machBand)
colorbar
colormap gray
hold off;

% 1b) Plot brightness as a function of horizontal position.
midRow = machBand(32, :);
figure(2);
hold on;
ylabel('Brightness')
xlabel('Horiz. position on Mach Band Diagram')
plot(midRow, 'b', 'LineWidth', 3);
axis([-10, 140, -5, 85]);
hold off;

% 1c) Create the receptive field of a retinal ganglion cell. 
thinVar = 2;
broadVar = 6;
S = 500;
IE = 1;

[X,Y] = meshgrid(-2:1:2);
R = sqrt(X.^2+Y.^2);
eGaus = normpdf(R,0,sqrt(thinVar));
iGaus = normpdf(R,0,sqrt(broadVar));
%eGaus = (1/(2*pi*thinVar))*exp(-R.^2/(2*thinVar));
%iGaus = (1/(2*pi*broadVar))*exp(-R.^2/(2*broadVar));

recField = S*(eGaus-IE*iGaus);

figure(3);
imagesc(recField);
hold on
colorbar
colormap gray


% 1d) Convolve your receptive field and input image to observe how the first stages of your
% visual system perceives the image.

convProd = conv2(machBand, recField, 'valid');
sliceMachBand = machBand(1,:)./max(machBand(1,:));
sliceConvProd = convProd(1,:)./max(convProd(1,:));
figure(4);
imagesc(convProd);
colorbar
colormap gray;

figure(5);
hold on;
ylabel('Brightness')
xlabel('Horiz. position on Mach Band Diagram')
plot(sliceMachBand, 'b', 'LineWidth', 4);
plot(sliceConvProd, 'k.-', 'LineWidth', 3);
%axis([-10, 140, -20, 110]);
legend('Mach Band Slice','Convoluted Slice');