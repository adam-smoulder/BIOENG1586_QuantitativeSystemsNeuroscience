% data structure: MTdata(trial#, FR, coherence)
% trial# = 1-100
% FR: 1 = left neuron FR, 2 = right neuron FR, 3 = result (1 = R, 0 = L)
% Coherence = 0, 0.04, 0.08, 0.16, 0.32, 0.64
% -Actual motion is always right (so 1 is correct prediction)
% -psychometric curve eq is: p = 1 - 0.5*exp(-(c/a)^b);
%   where a = threshold, b = sensitivity to coherence, c = coherence

% 1)

% finding accuracy at each coherence value
load('MTdata.mat');
ntrials = size(MTdata, 1);
decisions = squeeze(MTdata(:,3,:))'; % each row is one coherence value
accuracyPsycho = (1/ntrials).*(decisions*ones(ntrials,1)); % sums decisions and divides by length

% ran once to get coefficients:
% coherence=coherence+[0.0001,0,0,0,0,0];
% cftool(coherence, accuracyPsycho); % x's must be nonzero
%
% result - a = 0.11855, b = 3.2663


a = 0.119;
b = 3.266;
c = 0:0.001:1;
curve = 1-0.5*exp(-(c./a).^b);

figure(1);
axis([0,.70, 0.35, 1.10]);
set(gca,'FontSize',14)
hold on
plot(c, curve, 'LineWidth', 3);
plot(coherence, accuracyPsycho, 'k*', 'LineWidth', 10);
plot(a, 0.82, 'r*', 'LineWidth', 10);
title('Psychometric Curve')
xlabel('Coherence')
ylabel('Accuracy')
legend('Est. Curve', 'Data', 'Threshold')

% alpha estimate = 0.11855, beta estimate = 3.2663

%% trying stuff for problem 2

load('MTdata.mat');
ntrials = size(MTdata, 1);

% first let's just plot the distrubtions, just for 16% coh (4th index)
x=0:0.1:100;
neuron = 1;

rightNeuron = sort(MTdata(:,2,neuron));
rightNeuronmean = mean(rightNeuron);
rightNeuronvar = var(rightNeuron);
rightNeurongaussian=1/sqrt(2*pi)/sqrt(rightNeuronvar)*exp(-(x-rightNeuronmean).^2/2/sqrt(rightNeuronvar)/sqrt(rightNeuronvar));

leftNeuron = sort(MTdata(:,1,neuron));
leftNeuronmean = mean(leftNeuron);
leftNeuronvar = var(leftNeuron);
leftNeurongaussian=1/sqrt(2*pi)/sqrt(leftNeuronvar)*exp(-(x-leftNeuronmean).^2/2/sqrt(leftNeuronvar)/sqrt(leftNeuronvar));

figure(2);
hold on;
plot(x, rightNeurongaussian);
plot(x, leftNeurongaussian);

%% 2) Neurometric curve time
% For the curve:
% y-axis = odds null dist. is right of FR
% x-axis = odds pref dist. is right of FR

load('MTdata.mat');
ntrials = size(MTdata, 1);
FRmax = 100;

% for each coherence
% 1. iterate FR 0-100
% 2. if FR is in distribution, we find the index of the first value above
% it, then percentAbove = (100-(index-1))/100;

percentAbove = zeros(100,2,6); % (testing FR, L/R, Coh);
leftNeuron = squeeze(MTdata(:,1,:));
rightNeuron = squeeze(MTdata(:,2,:));

for i=1:length(coherence)
    sortedLeft = sort(squeeze(leftNeuron(:,i)));
    sortedRight = sort(squeeze(rightNeuron(:,i)));
    for j=1:FRmax
        LC_RN_Index = find(j < sortedLeft, 1, 'first');
        RC_RN_Index = find(j < sortedRight, 1, 'first');
        if ~isempty(LC_RN_Index)
            percentAbove(j,1,i) = 1-(LC_RN_Index-1)/length(leftNeuron);
        end
        
        if ~isempty(RC_RN_Index)
            percentAbove(j,2,i) = 1-(RC_RN_Index-1)/length(rightNeuron);
        end
    end
end

% ROC Curve

figure(3)
hold on
title('ROC Curve')
accuracyNeuro = zeros(6,1);
for i=1:length(coherence)
    plot(squeeze(percentAbove(:,1,i)), squeeze(percentAbove(:,2,i)), 'LineWidth', 3);
    accuracyNeuro(i) = trapz(squeeze(percentAbove(:,1,i)), squeeze(percentAbove(:,2,i)));
end
legend('0','4','8','16','32','64')

%values are negative because integration was backwards:
accuracyNeuro = accuracyNeuro.*(-1);


% Weibull curve fitting:
% ran once to get coefficients:
%
% coherence=coherence+[0.0001,0,0,0,0,0];
% cftool(coherence, accuracyNeuro); % x's must be nonzero
% 
% % result - a = 0.11855, b = 3.2663

% estimated a = 0.1745, b = 1.086
a = 0.1745;
b = 1.086;
c = 0:0.001:1;
curve = 1-0.5*exp(-(c./a).^b);

figure(4)
set(gca,'FontSize',14)
hold on
title('Neurometric Curve')
axis([0,1, 0.35, 1.10]);
hold on
plot(c, curve, 'LineWidth', 3)
plot(coherence, accuracyNeuro, 'k*', 'LineWidth', 10)
plot(a, 0.82, 'r*', 'LineWidth', 10)
xlabel('Coherence')
ylabel('Accuracy')
legend('Est. Curve', 'Data', 'Threshold')


%% 3)
load('MTdata.mat');
ntrials = size(MTdata, 1);
FRmax = 100;

RN_RC = zeros(size(MTdata,1),length(coherence)); % dims: trial, coherence
RN_LC = zeros(size(MTdata,1),length(coherence)); % dims: trial, coherence
LN_RC = zeros(size(MTdata,1),length(coherence)); % dims: trial, coherence
LN_LC = zeros(size(MTdata,1),length(coherence)); % dims: trial, coherence

percentAbove_RN_RC = zeros(size(MTdata,1),length(coherence)); % dims: trial, coherence
percentAbove_RN_LC = zeros(size(MTdata,1),length(coherence)); % dims: trial, coherence
percentAbove_LN_RC = zeros(size(MTdata,1),length(coherence)); % dims: trial, coherence
percentAbove_LN_LC = zeros(size(MTdata,1),length(coherence)); % dims: trial, coherence

nRight = zeros(1,6);
nLeft = zeros(1,6);

leftNeuron = MTdata(:,1,:);
rightNeuron = MTdata(:,2,:);

for i=1:length(coherence)
    
    % split data into left (0) and right (1) decisions 
    nRight(i) = sum(MTdata(:,3,i));
    nLeft(i) = 100-nRight(i);
    
    rightcount = 1;
    leftcount = 1;  
    
    for j=1:FRmax
        if MTdata(j,3,i) == 1
            RN_RC(rightcount,i) = MTdata(j,2,i);
            LN_RC(rightcount,i) = MTdata(j,1,i);
            rightcount = rightcount+1;
        else
            RN_LC(leftcount,i) = MTdata(j,2,i);
            LN_LC(leftcount,i) = MTdata(j,1,i);
            leftcount = leftcount+1;
        end
    end
    
    % ROC analysis for each choice of pref vs null
    % sorting
    RN_RC(1:nRight(i),i) = sort(RN_RC(1:nRight(i),i));
    LN_RC(1:nRight(i),i) = sort(LN_RC(1:nRight(i),i));
    RN_LC(1:nLeft(i),i) = sort(RN_LC(1:nLeft(i),i));
    LN_LC(1:nLeft(i),i) = sort(LN_LC(1:nLeft(i),i));

    % right-preferring neuron
    for j=1:FRmax
        RN_RC_Index = find(j < RN_RC(:,i), 1, 'first');
        LN_RC_Index = find(j < LN_RC(:,i), 1, 'first');
        
        if ~isempty(RN_RC_Index) % preferred
            percentAbove_RN_RC(j,i) = 1-(RN_RC_Index-1)/nRight(i);
        end
        
        if ~isempty(LN_RC_Index)
            percentAbove_LN_RC(j,i) = 1-(LN_RC_Index-1)/nRight(i);
        end
    end
    
    % left-preferring neuron
    for j=1:FRmax
        RN_LC_Index = find(j < RN_LC(:,i), 1, 'first');
        LN_LC_Index = find(j < LN_LC(:,i), 1, 'first');
        
        if ~isempty(RN_LC_Index) % preferred
            percentAbove_RN_LC(j,i) = 1-(RN_LC_Index-1)/nLeft(i);
        end
        
        if ~isempty(LN_LC_Index)
            percentAbove_LN_LC(j,i) = 1-(LN_LC_Index-1)/nLeft(i);
        end
    end
end

%nLeft = ntrials.*ones(1,6)-nRight;
accuracyNeuro2 = zeros(6,2); % left choice is :,1 while right choice is :,2


% plot left choice ROC
figure(5)
subplot(2,1,1)
set(gca,'FontSize',14)
hold on
title('ROC Curve, Left choice')
for i=1:length(coherence)
    plot(squeeze(percentAbove_RN_LC(1:FRmax,i)), squeeze(percentAbove_LN_LC(1:FRmax,i)), 'LineWidth', 3);
    accuracyNeuro2(i,1) = (-1)*trapz(squeeze(percentAbove_RN_LC(1:FRmax,i)), squeeze(percentAbove_LN_LC(1:FRmax,i)));
end
ylabel('Left-preferring FR')
xlabel('Right-preferring FR')
legend('0','4','8','16','32','64')
hold off

% plot right choice ROC
subplot(2,1,2)
set(gca,'FontSize',14)
hold on
title('ROC Curve, Right choice')
for i=1:length(coherence)
    plot(squeeze(percentAbove_LN_RC(1:FRmax,i)), squeeze(percentAbove_RN_RC(1:FRmax,i)), 'LineWidth', 3);
    accuracyNeuro2(i,2) = (-1)*trapz(squeeze(percentAbove_LN_RC(1:FRmax,i)), squeeze(percentAbove_RN_RC(1:FRmax,i)));
end
legend('0','4','8','16','32','64')
axis([0,1,0,1])
xlabel('Left-preferring FR')
ylabel('Right-preferring FR')

hold off

% add title/axis labels and such
% NEED TO SHOW HISTOGRAM PROBABILITY DISTRIBUTIONS against coherence!!! 

