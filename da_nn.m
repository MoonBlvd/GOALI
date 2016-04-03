function Li = da_nn(z, R)
% perform nearest-neighbor data association

global Param;
global State;
n = State.Ekf.nL;
Li = zeros(1,size(z,2));
X2 = chi2inv(0.99,2);
for i = 1:size(z,2)
    D2min = 1e15;%flintmax;
    nn = 0;
    for j = 1:n
        % calculate zhat and sigma for NN
        delta = State.Ekf.mu(2+2*j:3+2*j) - State.Ekf.mu(1:2);
        q = delta'*delta;
        zhat = [sqrt(q);minimizedAngle(atan2(delta(2),delta(1))-State.Ekf.mu(3))];
        dz = z(:,i) - zhat;
        dz(2) = minimizedAngle(dz(2));
        
        H = (1/q)*[-sqrt(q)*delta(1) -sqrt(q)*delta(2) 0 sqrt(q)*delta(1) sqrt(q)*delta(2);...
        delta(2) -delta(1) -q -delta(2) delta(1)];
        Fx = zeros(5,3+State.Ekf.nL*2);
        Fx(1:3,1:3) = eye(3);
        Fx(4:5, 2*j+2:2*j+3) = eye(2);
        %Fx = [[eye(3);zeros(2,3)],zeros(5,2*j-2),[zeros(3,2); eye(2)],zeros(5,2*(State.Ekf.nL-j))];
        H = H * Fx;
        P = H * State.Ekf.Sigma * H' + Param.R;
        % calculate mahalanobis distance
        Dij2 = (dz' / P) * dz;
        if Dij2 < D2min
            nn = j;
            D2min = Dij2;
        end
    end
    if D2min < X2
        Li(i) = nn;
    elseif D2min > X2+5
        Li(i) = 0;
    else
        Li(i) = -1;
    end
end


