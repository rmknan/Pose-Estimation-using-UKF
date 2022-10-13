function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    

n = length(covarEst);
alpha = 0.01;
beta =2;
k=2;
updn = 2*n+1;
n1 = eye(3)*0.01;
n2 = ones(3,1)*0.01;
lambda = (alpha^2)*(n+k)-n;
chol_f = chol(covarEst,"lower");

upoint(:,1) = uEst;
XYZ = [-0.04, 0.0, -0.03]';
Yaw = pi/4;

%Similar to sigma points
for i=1:n
    upoint(:,2*i) = uEst + sqrt(n+lambda)*chol_f(:,i);
    upoint(:,2*i+1) = uEst - sqrt(n+lambda)*chol_f(:,i);
end

avel = z_t(4:6);
zit = ones(3,length(upoint));

rbc = eul2rotm([-Yaw,0,pi], 'ZYX');
rcb = rbc';

bw = rcb*avel;
ske_ = XYZ; 
skm =  [0 (-1*ske_(3)) ske_(2);          %skew symmetric matrix
        ske_(3) 0 (-1*ske_(1));
       (-1*ske_(2)) ske_(1) 0];
% iterating through the update points to find the non linear function
for j = 1:length(upoint)
    r = upoint(4,j);
    p = upoint(5,j);
    y = upoint(6,j);
    rwb = eul2rotm([y,p,r],'ZYX')';
    vwb = upoint(7:9,j);

    Va = rbc*(rwb*vwb);
    Vb = rbc*skm;
    Vc = Vb*bw;
    vel = Va-Vc;
    zit(:,j) = vel + n2;
end

wom = lambda/(n+lambda);
wim = 1/(2*(n+lambda));
woc = wom + (1-alpha^2+beta);
wic = wim;

zut = zeros(3,1);
% finding Zut
for h = 1:updn
        if h == 1
            zut = zut + wom.*zit(:,h);      
        else
            zut = zut + wim.*zit(:,h);     
        end
 end
 ct =zeros(15,3);        %create the matrix
 st = zeros(3,3);
 % finding Ct,St and A
 for h = 1:updn
        if h == 1
            A = zit(:,h)- zut;
            ct = ct + woc.*(upoint(:,h) - uEst)*A.';  
            st = st + woc.*(A)*A.' ;     

        else 
            A = zit(:,h)- zut;
            ct = ct + wic.*(upoint(:,h) - uEst)*A.';     
            st = st + wic.*(A)*A.' ; 
        end
  end
% The final steps as given in the slides
 st = st + n1;
 Kt = ct/st;
 uCurr = uEst + Kt*(z_t(1:3)-zut);
 covar_curr = covarEst - Kt*st*(Kt.');


end

