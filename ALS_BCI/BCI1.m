% BCI PART 1
%% Continuous Decoding
load('continuous1.mat')

% 1) Visualizing data
%comet(kin(:,1),kin(:,2))

% color as time, normal x and y
x = squeeze(kin(:,1))';
y = squeeze(kin(:,2))';
z = zeros(size(x));
time = (1:70:3103*70-1)/1000; % in s
surface([x;x],[y;y],[z;z],[time;time],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
xlabel('X-Position (mm)')
ylabel('Y-Position (mm)')
set(gca,'FontSize',14)
c = colorbar;
c.Label.String = 'Time (s)';
c.Label.FontSize = 14;

% that looked horrible. Let's just do 2D x and y against time
figure();
subplot(2,1,1)
hold on
plot(time,x,'LineWidth',2)
ylabel('Horz. Dist (mm)')
set(gca,'FontSize',14)
axis([0 max(time) -inf inf])
subplot(2,1,2)
plot(time,y,'LineWidth',2)
ylabel('Vert. Dist (mm)')
xlabel('Time (s)')
set(gca,'FontSize',14)
axis([0 max(time) -inf inf])

% still ugly...overlay and zoom?
figure();
hold on
plot(time,x,'LineWidth',2)
ylabel('Dist (mm)')
set(gca,'FontSize',14)
plot(time,y,'LineWidth',2)
xlabel('Time (s)')
set(gca,'FontSize',14)
axis([100 120 -inf inf])
legend('Horiz.','Vert')

% okay that gets the idea!

%% 2) Linear filter time
linFilt1 = mvregress(rate,kin);

shift = 2;
rate2 = [rate(1:end-shift,:) ones(length(rate)-shift,1)];
kin2 = kin(1+shift:end,:);
linFilt2 = mvregress(rate2, kin2);

%% 3) 
estKin2 = rate2*linFilt2;

kin2Dists = sqrt(kin2(:,1).^2+kin2(:,2).^2);
estKin2Dists = sqrt(estKin2(:,1).^2+estKin2(:,2).^2);

% inspired by Sam Waters, plot real vs. estimated
hold on
plot(kin2Dists,estKin2Dists,'k.')
plot(0:1:30,0:1:30,'r.--','LineWidth',3)
ylabel('Dist (mm)')
set(gca,'FontSize',14)
xlabel('Distance From Origin (mm^2)')
ylabel('Estimated Distance From Origin (mm^2)')

% MSE calc, corr coeff
corrCoeff2 = corr(kin2Dists,estKin2Dists,'rows','pairwise');
MSE2 = immse(kin2Dists,estKin2Dists);

corrCoeff2_horz = corr(kin2(:,1),estKin2(:,1),'rows','pairwise');
MSE2_horz = immse(kin2(:,1),estKin2(:,1));

corrCoeff2_vert = corr(kin2(:,2),estKin2(:,2),'rows','pairwise');
MSE2_vert = immse(kin2(:,2),estKin2(:,2));

MSEs = [MSE2; MSE2_horz; MSE2_vert];
corrCoeffs = [corrCoeff2; corrCoeff2_horz; corrCoeff2_vert];

RowNames = [{'Total'};{'Horz'};{'Vert'}];
errorTable = table(MSEs, corrCoeffs,'RowNames',RowNames');
disp(errorTable)

%% 4) Same filter, new data
data3 = load('continuous2.mat');
kin3 = data3.kin;
rate3 = data3.rate;

rate3 = [rate3(1:end-shift,:) ones(length(rate3)-shift,1)];
kin3 = kin3(1+shift:end,:);
estKin3 = rate3*linFilt2;

kin3Dists = sqrt(kin3(:,1).^2+kin3(:,2).^2);
estKin3Dists = sqrt(estKin3(:,1).^2+estKin3(:,2).^2);

% inspired by Sam Waters, plot real vs. estimated
hold on
plot(kin3Dists,estKin3Dists,'k.')
plot(0:1:30,0:1:30,'r.--','LineWidth',3)
ylabel('Dist (mm)')
set(gca,'FontSize',14)
xlabel('Distance From Origin (mm^2)')
ylabel('Estimated Distance From Origin (mm^2)')

% MSE calc, corr coeff
corrCoeff3 = corr(kin3Dists,estKin3Dists,'rows','pairwise');
MSE3 = immse(kin3Dists,estKin3Dists);

corrCoeff3_horz = corr(kin3(:,1),estKin3(:,1),'rows','pairwise');
MSE3_horz = immse(kin3(:,1),estKin3(:,1));

corrCoeff3_vert = corr(kin3(:,2),estKin3(:,2),'rows','pairwise');
MSE3_vert = immse(kin3(:,2),estKin3(:,2));

MSE_Train = MSEs;
r_Train = corrCoeffs;
MSE_Test = [MSE3; MSE3_horz; MSE3_vert];
r_Test = [corrCoeff3; corrCoeff3_horz; corrCoeff3_vert];

RowNames = [{'Total'};{'Horz'};{'Vert'}];
errorTable = table(MSE_Train, MSE_Test, r_Train, r_Test,'RowNames',RowNames');
disp(errorTable)

%% 5) Come back to this!!!
load('continuous1.mat')

shifts = 0:1:20;
kinShift = cell(1,length(shifts));
rateShift = cell(1,length(shifts));
linFiltShift = cell(1,length(shifts));
estKinShift = cell(1,length(shifts));
kinDistShift = cell(1,length(shifts));
estKinDistShift = cell(1,length(shifts));
corrCoeffShift = cell(1,length(shifts));
MSEShift = cell(1,length(shifts));
corrCoeffShift_horz = cell(1,length(shifts));
MSEShift_horz = cell(1,length(shifts));
corrCoeffShift_vert = cell(1,length(shifts));
MSEShift_vert = cell(1,length(shifts));


for i=1:length(shifts)
    kinShiftTemp = kin(1+shifts(i):end,:);
    rateShiftTemp = [rate(1:end-shifts(i),:) ones(length(rate)-shifts(i),1)];
    %linFiltShiftTemp = mldivide(kinShiftTemp,rateShiftTemp)';
    linFiltShiftTemp = mvregress(rateShiftTemp,kinShiftTemp);
    %linFiltShiftTemp = pinv(rateShiftTemp)*kinShiftTemp;
    estKinShiftTemp = rateShiftTemp*linFiltShiftTemp;
    kinDistShiftTemp = sqrt(kinShiftTemp(:,1).^2+kinShiftTemp(:,2).^2);
    estKinDistShiftTemp = sqrt(estKinShiftTemp(:,1).^2+estKinShiftTemp(:,2).^2);
    corrCoeffShiftTemp = corr(kinDistShiftTemp,estKinDistShiftTemp,'rows','pairwise');
    MSEShiftTemp = immse(kinDistShiftTemp,estKinDistShiftTemp);
    corrCoeffShiftTemp_horz = corr(kinShiftTemp(:,1),estKinShiftTemp(:,1),'rows','pairwise');
    MSEShiftTemp_horz = immse(kinShiftTemp(:,1),estKinShiftTemp(:,1));
    corrCoeffShiftTemp_vert = corr(kinShiftTemp(:,2),estKinShiftTemp(:,2),'rows','pairwise');
    MSEShiftTemp_vert = immse(kinShiftTemp(:,2),estKinShiftTemp(:,2));
    
    % saving
    kinShift{i} = kinShiftTemp;
    rateShift{i} = rateShiftTemp;
    linFiltShift{i} = linFiltShiftTemp;
    estKinShift{i} = estKinShiftTemp;
    kinDistShift{i} = kinDistShiftTemp;
    estKinDistShift{i} = estKinDistShiftTemp;
    corrCoeffShift{i} = corrCoeffShiftTemp;
    MSEShift{i} = MSEShiftTemp;
    corrCoeffShift_horz{i} = corrCoeffShiftTemp_horz;
    MSEShift_horz{i} = MSEShiftTemp_horz;
    corrCoeffShift_vert{i} = corrCoeffShiftTemp_vert;
    MSEShift_vert{i} = MSEShiftTemp_vert;
end

% RowNames = [{'Total'};{'Horz'};{'Vert'}];
% errorTable = table([cell2mat(MSEShift)' cell2mat(MSEShift_horz)' cell2mat(MSEShift_vert)']','RowNames',RowNames');
% disp(errorTable)

MSEdata = [shifts' (shifts.*0.070)' cell2mat(MSEShift)' cell2mat(MSEShift_horz)' cell2mat(MSEShift_vert)']';
corrCoeffData = [shifts' (shifts.*0.070)' cell2mat(corrCoeffShift)' cell2mat(corrCoeffShift_horz)' cell2mat(corrCoeffShift_vert)']';

figure();
subplot(2,1,1)
hold on
plot(shifts,MSEdata(3,:),'LineWidth',2)
plot(shifts,MSEdata(4,:),'LineWidth',2)
plot(shifts,MSEdata(5,:),'LineWidth',2)
legend('MSE','Horiz MSE', 'Vert MSE')
set(gca,'FontSize',14)
xlabel('Shift (1 index = 70ms)')
ylabel('MSE')

subplot(2,1,2)
hold on
plot(shifts,corrCoeffData(3,:),'LineWidth',2)
plot(shifts,corrCoeffData(4,:),'LineWidth',2)
plot(shifts,corrCoeffData(5,:),'LineWidth',2)
legend('r','Horiz r', 'Vert r')
set(gca,'FontSize',14)
xlabel('Shift (1 index = 70ms)')
ylabel('Correlation Coefficient')


%% Part 2: Bayesian Integration
% 1) Some "tuning curves"

load('spikeCounts.mat') % dims: trial, neuron, target

FRAvgs = squeeze(mean(SpikeCounts,1));

figure();
hold on
for i=1:36
    subplot(6,6,i)
    bar(FRAvgs(2*i+1,:))
    axis([0.5 5.5 -inf inf])
    if i == 13
        ylabel('Average FR (sp/s)')
    end
    if i == 33
        xlabel('Target #')
    end
    title(['Neuron ' num2str(2*i+1)])
    set(gca,'FontSize',14)
end
hold off

%% 2/3) Training and such

holdout = squeeze(SpikeCounts(1,:,1));

% vert cat each thing for each target/page
trainingData = vertcat(squeeze(SpikeCounts(2:end,:,1)),...
                   squeeze(SpikeCounts(:,:,2)),...
                   squeeze(SpikeCounts(:,:,3)),...
                   squeeze(SpikeCounts(:,:,4)),...
                   squeeze(SpikeCounts(:,:,5)));
nTrialsPerTarg = [15, 16, 16, 16, 16];
classNames = vertcat(ones(15,1),2*ones(16,1),3*ones(16,1),4*ones(16,1),...
                                                         5*ones(16,1));

trainingData = trainingData+(10^(-20))*randn(size(trainingData));
bayesModel = fitcnb(trainingData, classNames);
reconstructedStuff = predict(bayesModel,holdout);

newTestData = vertcat(squeeze(SpikeCounts(:,:,1)),...
                   squeeze(SpikeCounts(:,:,2)),...
                   squeeze(SpikeCounts(:,:,3)),...
                   squeeze(SpikeCounts(:,:,4)),...
                   squeeze(SpikeCounts(:,:,5)));

moreReconStuff = predict(bayesModel,newTestData);
actual = vertcat(ones(16,1),2*ones(16,1),3*ones(16,1),4*ones(16,1),...
                                                          5*ones(16,1));
confMat = confusionmat(moreReconStuff, actual)./16;

figure();
imagesc(confMat)
hold on
axis xy
ylabel('Predicted Target')
xlabel('Actual Target')
set(gca,'FontSize',14)
c = colorbar;
c.Label.String = 'Proportion Predicted';
c.Label.FontSize = 14;
