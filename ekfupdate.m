function ekfupdate(z)
% EKF-SLAM update step for both simulator and Victoria Park data set

global Param;
global State;
Q = Param.R; % Sensing noise
% returns state vector indices pairing observations with landmarks
switch lower(Param.dataAssociation)
    case 'known'
        Li = da_known(z(3,:));
    case 'nn'
        Li = da_nn(z(1:2,:), Param.R);
    case 'jcbb'
        Li = da_jcbb(z(1:2,:), Param.R);
    otherwise
        error('unrecognized data association method: "%s"', Param.dataAssociation);
end
Fx = zeros(5,3+State.Ekf.nL*2);
Fx(1:3,1:3) = eye(3);
if strcmp(Param.updateMethod,'batch')
    for i = 1:length(Li)
        if Li(i) == 0
            initialize_new_landmark(z(:,i),Param.R);
            Li(i) = State.Ekf.nL;
        end
    end
end
if strcmp(Param.updateMethod,'sequential')
    for i = 1:length(Li)
        j = Li(i);
        if j == 0
            initialize_new_landmark(z(:,i),Param.R);
            j = State.Ekf.nL;
        elseif j == -1
            continue
        end
        delta = State.Ekf.mu(2+2*j:3+2*j) - State.Ekf.mu(1:2);
        q = delta'*delta;
        zhat = [sqrt(q);minimizedAngle(atan2(delta(2),delta(1))-State.Ekf.mu(3))];
        dz = z(1:2,i)-zhat;
        dz(2) = minimizedAngle(dz(2));
        H = (1/q)*[-sqrt(q)*delta(1) -sqrt(q)*delta(2) 0 sqrt(q)*delta(1) sqrt(q)*delta(2);...
                delta(2) -delta(1) -q -delta(2) delta(1)];
        Fx = zeros(5,3+State.Ekf.nL*2);
        Fx(1:3,1:3) = eye(3);
        Fx(4:5, 2*j+2:2*j+3) = eye(2);
        H = H * Fx;
        K = State.Ekf.Sigma*H'/(H*State.Ekf.Sigma*H'+Q);
        State.Ekf.mu = State.Ekf.mu+K*dz;
        State.Ekf.mu(3) = minimizedAngle(State.Ekf.mu(3));
        Joe = eye(size(State.Ekf.Sigma))-K*H;
        State.Ekf.Sigma = Joe*State.Ekf.Sigma*Joe' + K*Param.R*K';
    end
end
batchH = [];
batchDz = [];
batchQ = [];
if strcmp(Param.updateMethod,'batch')
    for i = 1:length(Li)
        j = Li(i);
        if j == -1
            continue
        end
        delta = State.Ekf.mu(2+2*j:3+2*j) - State.Ekf.mu(1:2);
        q = delta'*delta;
        zhat = [sqrt(q);minimizedAngle(atan2(delta(2),delta(1))-State.Ekf.mu(3))];
        dz = z(1:2,i)-zhat;
        dz(2) = minimizedAngle(dz(2));
        H = (1/q)*[-sqrt(q)*delta(1)   -sqrt(q)*delta(2)    0    sqrt(q)*delta(1)    sqrt(q)*delta(2);...
                   delta(2)            -delta(1)           -q    -delta(2)           delta(1)];
        Fx = zeros(3,3+State.Ekf.nL*2);
        Fx(1:3,1:3) = eye(3);
        Fx(4:5, 2*j+2:2*j+3) = eye(2);
        H = H * Fx;
        batchH = [batchH;H];
        batchDz = [batchDz;dz];
        batchQ = blkdiag(batchQ, Param.R);
    end
    S = batchH * State.Ekf.Sigma * batchH' + batchQ;
    K = State.Ekf.Sigma * batchH' / S;
    State.Ekf.mu = State.Ekf.mu + K * batchDz;
    State.Ekf.mu(3) = minimizedAngle(State.Ekf.mu(3));
    %Joeseph sigma
    Joe = eye(size(State.Ekf.Sigma))-K*batchH;
    State.Ekf.Sigma = Joe*State.Ekf.Sigma*Joe'+K*batchQ*K';

end
