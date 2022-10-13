 function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time 
    k=1;
beta =2;
alpha = 0.001;
n = length(uPrev);

lambda = (alpha^2)*(n+k)-n;
chol_f = chol(covarPrev,"lower");
X = zeros(n,2*n+1);
X(:,1) = uPrev;
sigmap = 2*n+1;
g = [0;0;-9.80]; 
an = [0.0; 0.00; 0.00];
gn = [0.0; 0.00; 0.00];
Q = eye(15)*0.001;

for i=1:n
    X(:,2*i) = uPrev + sqrt(n+lambda)*chol_f(:,i);
    X(:,2*i+1) = uPrev - sqrt(n+lambda)*chol_f(:,i);
end

yi = ones(15,length(uPrev));

for j = 1:sigmap
    gb = [X(10,j);X(11,j);X(12,j)];
    ab = [X(13,j);X(14,j);X(15,j)];
    r = X(4,j);
    p = X(5,j);
    y = X(6,j);


gm = [-sin(p)          0      1;
            cos(p)*sin(r)  cos(r)  0;
            cos(p)*cos(r) -sin(r) 0];

R = eul2rotm([y,p,r],'ZYX');

p1 = [X(7,j);X(8,j);X(9,j)];
p2 = gm\(angVel-gb-gn);
p3 = g+R*(acc-ab-an);

pm = [p1(1);
      p1(2);
      p1(3);
      p2(1);
      p2(2);
      p2(3);
      p3(1);
      p3(2);
      p3(3)];
xd = vertcat(pm,gn,an);
yi(:,j) = X(:,j) + (dt*xd);
end

wom = lambda/(n+lambda);
wim = 1/(2*(n+lambda));
woc = wom + (1-alpha^2+beta);
wic = wim;

uEst = zeros(15,1);
covarEst = zeros(15,15);

for h=1:sigmap
    if h==1
        uEst = uEst+wom.*yi(:,h);
    else
        uEst = uEst+wim.*yi(:,h);
    end
end

for h=1:sigmap
    if h==1
        A = yi(:,h)- uEst;
        covarEst = covarEst+woc.*(A)*(A.');
    else
        A = yi(:,h)- uEst;
        covarEst = covarEst+wic.*(A)*(A.');

    end
end
covarEst = covarEst+Q;
end

