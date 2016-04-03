function varargout = run(numSteps, choice, pauseLen, da,makeVideo)
% RUN PS3 EKF Feature-Based SLAM
%   RUN(ARG)
%   RUN(ARG, CHOICE, PAUSELEN)
%   RUN(ARG, CHOICE, PAUSELEN, DA)
%      ARG - is either the number of time steps, (e.g. 100 is a complete
%            circuit) or a data structure from a previous run.
%      CHOICE - is either 'sim' or 'vp' for simulator or Victoria Park
%               data set, respectively.
%      PAUSELEN - set to `inf`, to manually pause, o/w # of seconds to wait
%                 (e.g., 0.3 is the default)
%      DA - data assocation, is one of either:
%           'known' - only available in simulator
%           'nn'    - incremental maximum likelihood nearest neighbor
%           'nndg'  - nn double gate on landmark creation
%                     (throws away ambiguous observations)
%           'jcbb'  - joint compatability branch and bound
%
%   DATA = RUN(ARG, CHOISE, PAUSELEN, DA)
%      DATA - is an optional output and contains the data array generated
%             and/or used during the simulation.
%
%   Note: more parameters can be controlled in the run.m file itself via
%   fields of the Param structure.

%   (c) 2009-2015
%   Ryan M. Eustice
%   University of Michigan
%   eustice@umich.edu

addpath('./slamsim');
addpath('./vicpark');

if ~exist('pauseLen', 'var') || isempty(pauseLen)
    pauseLen = [];
end
if ~exist('makeVideo','var') || isempty(makeVideo)
    makeVideo = false;
end
clear global Param State Data;
global Param;
global State;
global Data;


% select which data association method to use in ekfupdate.m, choices are:
%   known - only available in simulator
%   nn    - incremental maximum likelhood nearest neighbor
%   nndg  - nn double gate on landmark creation (throws away ambiguous observations)
%   jcbb  - joint compatability branch and bound
if ~exist('da','var') || isempty(da)
    da = 'known';
end
Param.dataAssociation = da;

% select which update method to use in ekfupdate.m, choices are:
%   batch  - batch updates
%   seq    - sequential updates
Param.updateMethod = 'batch';
% Param.updateMethod = 'sequential';

% size of bounding box for VP data set plotting
Param.bbox = 0; % bbox = 20 [m] speeds up graphics

% Structure of global State variable
%===================================================
State.Ekf.t     = 0;          % time
State.Ekf.mu    = zeros(3,1); % robot initial pose
State.Ekf.Sigma = zeros(3,3); % robot initial covariance
State.Ekf.iR    = 1:3;        % 3 vector containing robot indices
State.Ekf.iM    = [];         % 2*nL vector containing map indices
State.Ekf.iL    = {};         % nL cell array containing indices of landmark i
State.Ekf.sL    = [];         % nL vector containing signatures of landmarks
State.Ekf.nL    = 0;          % scalar number of landmarks
%===================================================

switch lower(choice)
    case 'sim'
        Data = runsim(numSteps,pauseLen,makeVideo);
        if nargout > 0
            varargout{1} = Data;
        end
    case 'vp'
        runvp(numSteps,pauseLen,makeVideo);
    otherwise
        error('unrecognized selection "%s"', choice);
end
% Plot correlation coefficients
figure(2)
subplot(1,2,1);
symmetric = triu(State.Ekf.Sigma) + triu(State.Ekf.Sigma,1)';
corrM= corrcov(symmetric);
j = 1:State.Ekf.nL;
p = plot(1:State.Ekf.nL,corrM(1,2+2*j)); % correlation coefficients between X of Rob and X of Feature
hold on;
legendM = cell(1,State.Ekf.nL+1);
legendM{1} = ['X and features'];
for i = 1:State.Ekf.nL
	p = [p plot(1:State.Ekf.nL,corrM(2+2*i,2+2*j))];% correlation coefficients between X of all Features (including Auto parts)
    hold on;
    legendM{i+1} = ['Observation',num2str(i),' and all observations'];
end
legend(p,legendM);
title('Correlation coefficients ');
xlabel('Observation ID')
grid on;
subplot(1,2,2);
corrM = abs(eye(size(State.Ekf.Sigma, 1)) - corrM);
image = mat2gray(corrM);
imshow(image)
title('The correlation matrix in gray level');

% The determinant of the feature covariance matrix for all map features
figure(4)
detCov = zeros(1, State.Ekf.nL);
for i=1:State.Ekf.nL
    detCov(i) = det(State.Ekf.Sigma(2*i+2:2*i+3, 2*i+2:2*i+3));
end
plot(1:State.Ekf.nL, detCov.^(1/4),'--s','MarkerSize',20,'LineWidth',2);
h = ylabel('$\Sigma_{ii}^{\frac{1}{4}}$');
set(h, 'interpreter', 'latex')
xlabel('Observation id')
title('The determinant of the covariance matrix');