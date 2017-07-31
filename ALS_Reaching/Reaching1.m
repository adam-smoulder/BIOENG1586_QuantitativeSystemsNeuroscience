%% 1) load/organize data

% hhp = horizontal hand position
% vhp = vertical hand position
% zhp = z-axis hand position
% hep = horizontal eye position
% vep = vertical eye position

load('neuralData.mat')

%% 2) find target positions

% initial hand pos = (0,0), eye pos = (20,20)
TrialParams = [R.TrialParams];
allTarget1X = [TrialParams.target1X];
allTarget1Y = [TrialParams.target1Y];
allTargetXY = [allTarget1X' allTarget1Y'];

% find the 8 targets by inspection:
% disp(allTargetXY);

targets = [-98 -17; -86 50; -64 -76; -34 93; 34 93; 64 -76; 86 50; 98 -17];

%% 3) Generate plot for hand movement
%stuffWeCareAbout = struct();
figure(66);
hold on
for i = 1:length(R)
    % 1) Find the peak speed.
    horizHandPos = R(i).hhp;
    vertHandPos = R(i).vhp;
    R(i).handSpeed = sqrt(diff(horizHandPos).^2 + diff(vertHandPos).^2);
    R(i).peakSpeed = max(R(i).handSpeed(1000:end));
    % 2) Choose a fraction of the peak speed - usually, 20% is about right, but play with it and see what looks good
    R(i).threshSpeed = 0.2*R(i).peakSpeed;
    % 3) When the hand velocity first exceeds that threshold, call that the movement onset time. Limit your search to the time from timeGoCue until timeTargetAcquire.
    % adding 500 to ignore first bit of each trial
    R(i).timeMoveOnset = find(R(i).handSpeed(500:end) > R(i).threshSpeed, 1, 'first')+500;
    % 4) For the movement end time, it is the same calculation, applied looking forward from the peak velocity. (Movement end times are less important for this assignment.)
    R(i).timeMoveEnd = find(R(i).handSpeed > R(i).threshSpeed, 1, 'last');
    
    plot(horizHandPos(R(i).timeMoveOnset : R(i).timeMoveEnd), ...
        vertHandPos(R(i).timeMoveOnset : R(i).timeMoveEnd), ...
        'k.','Color',ColorSelection(allTargetXY(i,1)))
end

% plotting circles at target points
for i = 1:size(targets,1)
    plot(targets(i,1),targets(i,2),'k+-','LineWidth',5)
end

xlabel('Horiz. Pos (mm)')
ylabel('Vert. Pos (mm)')
set(gca,'FontSize',14)
hold off

% very odd data point comes from trial 1739

%% 5: Quantify reaction time

% Presentation of go cue 0> time movement begins, so
% timeMoveOnset-timeGoCue

rts = [R.timeMoveOnset] - [R.timeGoCue];

meanRT = mean(rts);
stdRT = std(rts);

% we got 289.7+/-22.45 ms

%% 6: Reaction time for each cue? Are they different?

dirRTs = zeros(length(R),8); % one column for each cue
dirCounter = ones(1,8);      % for counting trials to each target
cue = zeros(1,length(R));      % assigning cue for each trial

% assigning RTs
for i=1:length(rts)
    cue(i) = CueIdentifier(allTargetXY(i,1));  % determine direction cue
    dirRTs(dirCounter(cue(i)),cue(i)) = rts(i);   % assign RT to proper cue
    dirCounter(cue(i)) = dirCounter(cue(i))+1;    % increment direction cue count
end

% analyzing RTs
meanDirRTs = zeros(1, length(dirCounter));
stdDirRTs = zeros(1, length(dirCounter));
for i=1:length(dirCounter)
    meanDirRTs(i) = mean(dirRTs(1:dirCounter(i),i));
    stdDirRTs(i) = std(dirRTs(1:dirCounter(i),i));
end

% performing one-way ANOVA for significance in mean difference
[p, tab, stats] = anova1(rts, cue);

%% displaying the myriad of analyses we just did:

% table for means/stdevs
targNames = cell(8,1);
for i=1:length(targNames)
    targNames{i} = mat2str(targets(i,:));
end
Mean = meanDirRTs'; % for table
StDev = stdDirRTs';
statsTable = table(Mean, StDev,'RowNames',targNames');
disp(statsTable)

% quick ANOVA results;
disp(['p-value representing if all directions are equal: ' num2str(p)])
disp(tab);

% For checking individual significant differences:
multcompare(stats)


%% 7) Monitor refresh rate latency

% looking at difference between timeGoCue and timeGoCuePHOTO
monitorLatency = [R.timeGoCuePHOTO] - [R.timeGoCue];

meanLatency = mean(monitorLatency);
stdLatency = std(monitorLatency);
refreshRate = 1000/meanLatency;     %1000 is for ms -> s

% monitor latency is 18.9+/-4.89ms, refresh rate ~= 53Hz

%% 8) Overlaying hand plot with eye endpoint

figure(66);
hold on
for i = 1:length(R)
    % 1) Find avg. position
    horizEyePos = R(i).hep;
    vertEyePos = R(i).vep;
    movementHEP = horizEyePos(R(i).timeMoveOnset:R(i).timeMoveEnd);
    movementVEP = vertEyePos(R(i).timeMoveOnset:R(i).timeMoveEnd);
    avgHorizEyePos = mean(movementHEP);
    avgVertEyePos = mean(movementVEP);

    % plot average eye location for trial
    plot(avgHorizEyePos, avgVertEyePos, 'k+',...
        'Color',ColorSelection(allTargetXY(i,1))./1.4,'LineWidth',2)
end

xlabel('Horiz. Pos (mm)')
ylabel('Vert. Pos (mm)')
set(gca,'FontSize',14)