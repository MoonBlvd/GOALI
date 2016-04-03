function Li = da_jcbb(z, R)
% perform joint-compatability branch and bound data association

global Param;
global State;
global m Best;
m = size(z,2);
H = [];
Best = zeros(1,m)
i = 1;
jcbb(H,i,z);
Li = Best;
end
function jcbb(H, i, z)
global m Best;
global State;
    if i > m
        if pairings(H) > pairings(Best)
            Best = H;
        end
    else
         for j = 1:State.Ekf.nL
             if individual_compatibility(z(:,i), j) && joint_compatibility(H, z, i, j)
                 jcbb([H j], i+1, z);
             end
         end
    end
    if pairings(H) + m - i >= pairings(Best)
        jcbb([H 0], i+1, z);
    end
end

function IC = individual_compatibility(z,j)
[~,~,D2] = mahalanobis(z,j);    
IC = D2 < chi2inv(0.99,2)
end

function JC = joint_compatibility(H, z, i, j)
    global Param;
    global State;
    matrixH = [];
    matrixDZ = [];
    Q = [];
    for k = 1 : i -1     
        if H(k) ~= 0
            [HH,dz,~] = mahalanobis(z(:,k),H(k));
            matrixH = [matrixH;HH];
            matrixDZ = [matrixDZ; dz];
            Q = blkdiag(Q, Param.R);
        end
    end
        [HH,dz,~] = mahalanobis(z(:,i),j);
        matrixH = [matrixH;HH];
        matrixDZ = [matrixDZ; dz];
        Q = blkdiag(Q, Param.R);
        CH = matrixH*State.Ekf.Sigma*matrixH'+Q;
        D2 = (matrixDZ' / CH)*matrixDZ;
        JC = D2 < chi2inv(0.99,2*i)
end

function p = pairings(H)
p = length(find(H>0));
end

function [H dz Dij2] = mahalanobis(z,j)
    global Param;
    global State;
    delta = State.Ekf.mu(2+2*j:3+2*j) - State.Ekf.mu(1:2);
    q = delta'*delta;
    zhat = [sqrt(q);minimizedAngle(atan2(delta(2),delta(1))-State.Ekf.mu(3))];
    dz = z - zhat;
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
end