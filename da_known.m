function Li = da_known(z)
% EKF-SLAM data association with known correspondences

global Param;
global State;
Li = zeros(1,size(z,2));
for i = 1:size(z,2)
    id = find(State.Ekf.sL == z(i));
    if id
        Li(i) = id;
    else
        Li(i) = 0;
    end
end
