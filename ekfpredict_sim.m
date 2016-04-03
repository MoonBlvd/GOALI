function ekfpredict_sim(u)
% EKF-SLAM prediction for simulator process model

global Param;
global State;
theta = State.Ekf.mu(3);
n = size(State.Ekf.Sigma,1);
% --------------------------------------------
% Prediction step
Fx = zeros(3,3+State.Ekf.nL*2);
Fx(1:3,1:3) = eye(3);
tempG = [0 0 (-u(2))*sin(u(1)+theta);...
    0 0 (u(2))*cos(u(1)+theta);...
    0 0 0];
G = eye(3+State.Ekf.nL*2)+Fx'*tempG*Fx;
V = [(-u(2))*sin(u(1)+theta) cos(u(1)+theta) 0;...
    (u(2))*cos(u(1)+theta) sin(u(1)+theta) 0;...
    1 0 1];
alphas = Param.alphas;
M = [alphas(1)*u(1)^2+alphas(2)*u(2)^2 0 0;...
        0 alphas(3)*u(2)^2+alphas(4)*(u(1)^2+u(3)^2) 0;...
        0 0 alphas(1)*u(3)^2+alphas(2)*u(2)^2];
% --------------------------------------------
% EKF prediction of mean and covariance
State.Ekf.mu = prediction(State.Ekf.mu,u);
State.Ekf.mu(3) = minimizedAngle(State.Ekf.mu(3));
R = V*M*V';
State.Ekf.Sigma =  G*State.Ekf.Sigma*G'+Fx'*R*Fx;
end