clear; % Clear variables
addpath('../data')
datasetNum = 4; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime, proj2Data] = init(datasetNum);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = 0.01*eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0;
vel = proj2Data.linearVel;
angVel2 = proj2Data.angVel;
%% Calculate Kalmann Filter
for i = 1:length(sampledTime)
    %% FILL IN THE FOR LOOP
        lv = vel(i,:);
        av = angVel2(i,:);
        z_t = [lv';av'];
        dt = sampledData(i).t-prevTime;      %difference in time 
        [covarEst,uEst] = pred_step(uPrev, covarPrev, sampledData(i).omg,sampledData(i).acc, dt);
        [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst);
        savedStates(:,i) = uCurr;   %store in savedStates
        prevTime = sampledData(i).t; 
        uPrev = uCurr;              %update uPrev
        covarPrev = covar_curr;     %update covarPrev
end

plotData(savedStates, sampledTime, sampledVicon, 2, datasetNum);