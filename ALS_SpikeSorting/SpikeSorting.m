load('waveforms.mat')
APsnippet = data.wf;
APtimestamp = data.stamps;

%% 1) plotting 100 waveforms
% snippet = 48 pts / 2ms -> sampling = 24kHz = 24 pts/ms

t = 0:1/24:2-1/24; % in ms
figure();
hold on
for i = 1:floor(80397/100+5):80397
    plot(t,APsnippet(i,:))
end
xlabel('Time (ms)')
ylabel('Amplitude')
set(gca,'FontSize',14)
hold off

%% 2) PCA
[coeff, score, eigenvals, tsquaredstat, percentExp] = pca(APsnippet);
pc1vals = APsnippet*coeff(:,1);
pc2vals = APsnippet*coeff(:,2);

% plotting
pcaHisto = hist3([pc1vals pc2vals], [50 50]);

figure()
imagesc(pcaHisto)
set(gca,'FontSize',14)
xlabel('PC 1')
ylabel('PC 2')
colorbar

%% 3 & 4) K-means clustering
k_value = 3;

clusters = kmeans([pc1vals pc2vals],k_value);
color = {'b.-','r.-','g.-'};

figure()
hold on
for i = 1:floor(80397/100+2):80397
    plot(t, APsnippet(i,:), color{clusters(i)})
end

xlabel('Time (ms)')
ylabel('Amplitude')
set(gca,'FontSize',14)

%% 5) Rotate PCA components
% all PC coeffs should be normalized to one, so they are the unit vecs.
% Rotation would just be plotting their values along the x/y axis used
% originally.

figure()
hold on
for i = 1:3
    plot(t, coeff(:,i),color{i},'LineWidth',3)
end
legend('PC 1', 'PC 2', 'PC 3')
xlabel('Time (ms)')
ylabel('Amplitude')
set(gca,'FontSize',14)

%% 6) Eigenspectrum evaluation
% visual search for the elbow
figure()
plot(percentExp,'LineWidth',3)
xlabel('Dimension #')
ylabel('Percent Variance Accounted For (%)')
set(gca,'FontSize',14)

% estimate: 5

cumPercentExp = cumsum(percentExp);  % cumulative summing
nintyPercentPoint = find(cumPercentExp > 90,1,'first');

figure()
plot(cumPercentExp,'LineWidth',3)
hold on
plot(1:48,90*ones(1,48),'r.--')
xlabel('# of Dimensions used')
ylabel('Percent Variance Accounted For (%)')
set(gca,'FontSize',14)

%% 7) Final sorting

% sampling was 24kHz
% max timestamp value is 4039.0
% you said the experiment was over an hour long
% therefore, timestamp must be in seconds
% 10s = 10000ms = 240000 points for snippets

tstart =  1900;
tend = 1910;
timestampStartIndex = find(APtimestamp >= tstart, 1, 'first');
timestampEndIndex = find(APtimestamp >= tend, 1, 'first');

% sampling at 24kHz
fs = 24000;
t = 1900 : 1/fs : 1950; % adding 47 pts in case AP is longer than t

% for plot, need to properly place in time and also color
% also need to only plot part of those exceeding our window
figure()
subplot(k_value+2,1,k_value+2)
plot(t(1:240000),zeros(1,240000),'k.-')
hold on
for i = timestampStartIndex : 1 : timestampEndIndex % for each timestamp
    startTime = find(t >= APtimestamp(i),1,'first');          % find its time-value
    timestampTime = t(startTime : startTime+47);    % set up 2ms
    plot(timestampTime, APsnippet(i,:), color{clusters(i)}) % plot
end
axis([1900,1910,-inf,inf])
xlabel('Time (s)')
ylabel('Combined APs')
set(gca,'FontSize',14)
hold off

% simplified version, combined rasters
clear i
subplot(k_value+2,1,k_value+1)
plot(t,zeros(1,length(t)),'k.-')
hold on
for i = timestampStartIndex : 1 : timestampEndIndex    % for each timestamp
    startTime = find(t >= APtimestamp(i),1,'first');   % find its time-value
    timestampTime = t(startTime+10 : startTime+14);    % set up short window
    plot(timestampTime, [0 1 1 1 0], color{clusters(i)}, 'LineWidth',1) % plot
end
axis([tstart,tend,0,1.25])
ylabel('Combined')
set(gca,'FontSize',14)
hold off

% separated rasters (lazy, unclean coding...sorry)
clear j
for j = 1:k_value
    subplot(k_value+2,1,j)
    plot(t,zeros(1,length(t)),'k.-')
    hold on
    for i = timestampStartIndex : 1 : timestampEndIndex    % for each timestamp
        if clusters(i) == j
            startTime = find(t >= APtimestamp(i),1,'first');   % find its time-value
            timestampTime = t(startTime+10 : startTime+14);    % set up short window
            plot(timestampTime, [0 1 1 1 0], color{clusters(i)}, 'LineWidth',1) % plot
        end
    end
    axis([tstart,tend,0,1.25])
    ylabel(['Neuron ' num2str(j)])
    set(gca,'FontSize',14)
    hold off
end

%% 8) ISI histograms (once again, ugly coding, very sorry)
neuron1times = APtimestamp(clusters == 1);
neuron2times = APtimestamp(clusters == 2);
neuron3times = APtimestamp(clusters == 3);

neuron1ISI = zeros(1,length(neuron1times)-1);
for i = 1:length(neuron1ISI)
    neuron1ISI(i) = neuron1times(i+1)-neuron1times(i);
end
neuron2ISI = zeros(1,length(neuron2times)-1);
for i = 1:length(neuron2ISI)
    neuron2ISI(i) = neuron2times(i+1)-neuron2times(i);
end
neuron3ISI = zeros(1,length(neuron3times)-1);
for i = 1:length(neuron3ISI)
    neuron3ISI(i) = neuron3times(i+1)-neuron3times(i);
end

% plotting (edit values for diff bin histograms)

figure()
subplot(3,1,1)
hold on
histogram(neuron1ISI,1000)
ylabel('Neuron 1')
axis([0,0.5,-inf,inf])
set(gca,'FontSize',14)
hold off
subplot(3,1,2)
hold on
histogram(neuron2ISI,1000)
axis([0,0.5,-inf,inf])
ylabel('Neuron 2')
xlabel('Interspike Interval (s)')
set(gca,'FontSize',14)
hold off

subplot(3,1,3)
hold on
histogram(neuron3ISI,10000)
axis([0,0.1,-inf,inf])
ylabel('Neuron 3')
xlabel('Interspike Interval (s)')
set(gca,'FontSize',14)
hold off







