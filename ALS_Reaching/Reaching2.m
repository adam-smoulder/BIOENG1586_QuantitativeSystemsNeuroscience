% RUN REACHING1 BEFORE THIS OR VARIABLES WILL NOT BE PRESENT

%% 1) Plotting rasters and hand movement activity

% selected trial
trialData = R(1011);
nCells = length(trialData.cells);

figure();
subplot(3,1,1)
hold on
raster = MultiNeuronRaster(trialData.trialLength, trialData.cells,1);
plot(squeeze(trialData.timeCueOnset.*ones(1,2*nCells+1)), 0:0.5:nCells, ...
    'k.--', 'LineWidth',1.5)
plot(squeeze(trialData.timeGoCuePHOTO.*ones(1,2*nCells+1)), 0:0.5:nCells, ...
    'k.--', 'LineWidth',1.5)
ylabel('Cell Rasters')
set(gca,'FontSize',14)
hold off

subplot(3,1,2)
hold on
plot(trialData.hhp,'LineWidth',2)
plot(trialData.vhp,'LineWidth',2)
plot(squeeze(trialData.timeCueOnset.*ones(1,10+1)), -80:10:20, ...
    'k.--', 'LineWidth',1.5)
plot(squeeze(trialData.timeGoCuePHOTO.*ones(1,10+1)), -80:10:20, ...
    'k.--', 'LineWidth',1.5)
axis([0,trialData.trialLength,-inf,inf])
set(gca,'FontSize',14)
ylabel('Hand Pos. (mm)');
hold off

subplot(3,1,3)
hold on
plot(trialData.hep,'LineWidth',2)
plot(trialData.vep,'LineWidth',2)
plot(squeeze(trialData.timeCueOnset.*ones(1,12+1)), -80:10:40, ...
    'k.--', 'LineWidth',1.5)
plot(squeeze(trialData.timeGoCuePHOTO.*ones(1,12+1)), -80:10:40, ...
    'k.--', 'LineWidth',1.5)
axis([0,trialData.trialLength,-inf,inf])
xlabel('Time (ms)')
set(gca,'FontSize',14)
ylabel('Eye Pos. (mm)');
legend('Horiz.','Vert.','Location','Best')
hold off

%% 2) PSTHs
% timeCueOnset-300, timeCueOnset+600

targ7Reaches = R(cue == 7); % only the [86,50] target
spikeData = zeros(length(targ7Reaches),length(targ7Reaches(1).cells),901); % dims: trial, neuron, time
windSize = 2;
%wind = gaussmf(0:windSize-1, [floor(windSize/4) floor(windSize/2)]);
wind = normpdf(0:1:900,0,1);
sepPSTHs = zeros(length(targ7Reaches),length(targ7Reaches(1).cells),901); % dims: trial, neuron, time

for i=1:size(spikeData,1)
    startTime = targ7Reaches(i).timeCueOnset-300;
    endTime = startTime+900;
    currentCells = targ7Reaches(i).cells;
    for j=1:length(currentCells)
        spikeTimes = currentCells(j).spikeTimes;
        for k=1:length(spikeTimes)
            if spikeTimes(k) > startTime && spikeTimes(k) < endTime
                spikeData(i,j,floor(spikeTimes(k))-startTime+1) = 1;
            end
        end
        convProd = conv(squeeze(spikeData(i,j,:)),squeeze(wind));
        sepPSTHs(i,j,:) = convProd(1:901);
    end
end

neuronPSTHs = squeeze(sum(sepPSTHs,1));

time = -300:1:600;
figure();
for i=1:length(targ7Reaches(1).cells)
    %scalingConst = max(wind)/length(leftReaches);
    subplot(length(targ7Reaches(1).cells),1,i)
    hold on
    plot(time, neuronPSTHs(i,:),'k.','Color',ColorSelection2(i),'LineWidth',2)
    ylabel(['Neuron ' num2str(i) ' Firing'])
    axis([-300,600,-inf,inf])
    set(gca,'FontSize',14)
    hold off
end

xlabel('Time (ms)')


%% 3) Tuning curves

% redefining spike times for all directions
allSpikeData = zeros(length(R),length(R(1).cells),501); % dims: trial, neuron, time
for i=1:size(allSpikeData,1)
    startTime = R(i).timeCueOnset-300;
    endTime = startTime+900;
    currentCells = R(i).cells;
    for j=1:length(currentCells)
        spikeTimes = currentCells(j).spikeTimes;
        for k=1:length(spikeTimes)
            if (spikeTimes(k) > startTime+100) && (spikeTimes(k) < startTime+600)
                allSpikeData(i,j,floor(spikeTimes(k))-startTime+100+1) = 1;
            end
        end
    end
end


dirCounters = zeros(1,8);
dirSums = zeros(4,8);

% summing spikes over trials and sorting to direction
for j=1:length(R) % for each trial
    currentCue = cue(j);
    trialSpikes = sum(squeeze(allSpikeData(j,:,:)),2);
    dirSums(:,currentCue) = dirSums(:,currentCue)+trialSpikes;
    dirCounters(currentCue) = dirCounters(currentCue)+1;
end

% averaging
tuningCurve = zeros(4,8);
for i=1:length(dirCounters)
    tuningCurve(:,i) = dirSums(:,i)./dirCounters(i);
end

% found using Excel
targAngles = [180+9.8411,180-30.1735,180+49.8991,180-69.9180,69.9180,360-49.8991,30.1735,360-9.8411];

% from CFtool:
fitCurves = zeros(4,360); % as opposed to overweight curves
a = [1.378, 5.728, 5.991, 0.832];
b = [1.235, 1, 1, 1];
c = [440, -32.73, -33.41, 9.377];
d = [15.21, 10.27, 9.297, 11.86];
x = 1:1:360;

for i=1:nCells
    fitCurves(i,:) = a(i).*cosd(b(i).*x+c(i).*ones(1,360))+d(i).*ones(1,360);
end

% plotting tuning curves
figure();
for i=1:nCells
    subplot(2,2,i)
    hold on
    plot(targAngles,squeeze(tuningCurve(i,:)),'ko','Color',ColorSelection2(i),'LineWidth',5)
    plot(x, fitCurves(i,:), 'k.-','LineWidth',2)
    ylabel(['Neuron' num2str(i) 'FR (sp/s)'])
    xlabel('Target Angle (deg)')
    set(gca,'FontSize',14)
    axis([0 360 0 20])
    hold off
end


