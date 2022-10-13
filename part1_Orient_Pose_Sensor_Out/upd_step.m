function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    R = eye(6)*.01; % the noise
    C =        [1,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0;
                0,1,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0;
                0,0,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0;
                0,0,0, 1,0,0, 0,0,0, 0,0,0, 0,0,0;
                0,0,0, 0,1,0, 0,0,0, 0,0,0, 0,0,0;
                0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0]; % The C matrix for the update
  % writing the update equations from the slide
    Kt_gain = covarEst*(C.')/(C*covarEst*(C.') + R);
    covar_curr = covarEst - Kt_gain*C*covarEst;
    uCurr = uEst + Kt_gain*(z_t-C*uEst);
end

