function ekfpredict_vp(u,dt)
% EKF-SLAM prediction for Victoria Park process model

global Param;
global State;

x = State.Ekf.mu(1);
y = State.Ekf.mu(2);
theta = State.Ekf.mu(3);
a = Param.a; 
b = Param.b;
L = Param.L;
% Prediction step
%---------------------------------------------------------------
% EKF prediction of mean

tempMu = [u(1)*cos(theta)-(u(1)/L)*tan(u(2))*(a*sin(theta)+b*cos(theta));
             u(1)*sin(theta)+(u(1)/L)*tan(u(2))*(a*cos(theta)-b*sin(theta));
            (u(1)/L)*tan(u(2))].*dt;
State.Ekf.mu(State.Ekf.iR) = [x;y;theta] + tempMu; 
State.Ekf.mu(3) = minimizedAngle(State.Ekf.mu(3)); 
% --------------------------------------------
% EKF prediction of covariance
Fx = zeros(3,3+State.Ekf.nL*2);
Fx(1:3,1:3) = eye(3);
tempG = [0 0 -dt*(u(1)*sin(theta)+(u(1)/L)*tan(u(2))*(a*cos(theta) - b*sin(theta)));...
    0 0 dt*(u(1)*cos(theta)-(u(1)/L)*tan(u(2))*(a*sin(theta) + b*cos(theta)));...
    0 0 0];
G = eye(size(State.Ekf.Sigma))+Fx'*tempG*Fx;
V = [cos(theta)-(1/L)*tan(u(2))*(a*sin(theta)+b*cos(theta)) (u(1)/L)*sec(u(2)).^2*(a*sin(theta)+b*cos(theta));
     sin(theta)+(1/L)*tan(u(2))*(a*cos(theta)-b*sin(theta)) (u(1)/L)*sec(u(2)).^2*(a*cos(theta)-b*sin(theta));
     (1/L)*tan(u(2)) (u(1)/L)*sec(u(2)).^2];
State.Ekf.Sigma =  G*State.Ekf.Sigma*G' + G*Fx'*Param.Qf*Fx*G' + Fx'*V*Param.Qu*V'*Fx;
% sigma0 = State.Ekf.Sigma
end