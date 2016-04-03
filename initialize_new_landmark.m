function initialize_new_landmark(z,R)

global Param;
global State;
tempMu = State.Ekf.mu(1:2)+[z(1)*minimizedAngle(cos(z(2)+State.Ekf.mu(3)));z(1)*sin(minimizedAngle(z(2)+State.Ekf.mu(3)))];
State.Ekf.mu = [State.Ekf.mu;tempMu];
tempSigma = flintmax.*eye(2);% Initialized the diagnal of the Sigma of new landmark to be infinite(very large number)
State.Ekf.Sigma = blkdiag(State.Ekf.Sigma,tempSigma);% augmentation
State.Ekf.nL = State.Ekf.nL+1; % add one more landmark
State.Ekf.sL = [State.Ekf.sL z(3)]; % add the landmark ID
