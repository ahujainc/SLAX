%-------------------------------------------------------------------------
% Function Name: mainfun_SLAX
% Date last modified: December 12, 2018
% Author: Vishal Ahuja
% PURPOSE: SLAX algorithm

% This function is the main function and calls the following functions:
% 1) fun_ap_sim_B,
% 2) fun_ap_imp_sim_B, and
% 3) fun_grid_gen (if needed for further grid refinement)

% It takes on the the parameter m (where m stands for the iteration number).
%-------------------------------------------------------------------------
function mainfun_SLAX(m)
% clear all;
% close all;
% clc;
% rng('shuffle')
% m = 1; % number of iterations to run;
global states num_states ab pat T I actions Del_T alpha srun1 srun2 ub_list H B d
% Initialize the computation times
elapsedtime1=0; % computation time for step 1
elapsedtime2=0; % computation time for step 2 (for all possible scenarios)

% -----------------Basic parameters--------------------------------
pat     = 20;       % total number of patients per period
T       = 10;       % number of time periods in the trial
N       = pat*T;    % number of patient observations
actions = pat+1;    % no. of actions (exceeds patient by 1)
Tx      = 2;        % number of treatments

% --------- Calibrated Parameters for SLAX---------------
% **********Can be Adjusted**********
alpha   = 0.4;      % weight on Std Dev in the objective function
srun1   = 5%3000;  % Simulation runs for outer loop (s'|s,a)
srun2   = 5%50;    % Simulation runs for inner loop (s'|s,a*)
H       = 3;        % length of restricted horizon for implementing SLAX
bsize   = 5;        % How many blocks of actions to consider
B       = pat/bsize;% Block Size


% ---Choice of starting priors for the treatments
list    = [4 1; 6 2; 1 0.5; 2 1; 0.5 0.5; 1 1; 2 2; 4 4; 6 6; 1 2; 0.5 1; 2 6; 1 4]; % List of all priors
sz      = size(list,1); % total number of priors

% The following refers to the pair number in the above "list" 
% For e.g. "pr1=6" refers to the sixth pair for tx A-> (alpha_A, beta_A) = (1,1)
% **********Can be Adjusted**********
pr1     = 6; % Prior for treatment A
pr2     = 6; % Prior for treatment B
ab = [list(pr1,:); list(pr2,:)]; % starting priors for two treatments

%---Create the set of simplices and grid points for base resolution (Z=0)--
d       = 3;        % dimensions of the state = I-1
I       = 4;        % Number of conditions: {As, Af, BS, Bf}

ends    = [0 1];    
points1 = zeros(2^d,d); % generate 8 corner points of the simplex

for i=0:d-1
    rep = repmat(ends,2^(d-i-1),2^i);
    points1(:,i+1) = rep(:);
end
states      = points1;          % set of eight grid points: H^tilde_Z=0
num_p       = size(points1,1);   

% Minimal Triangulation: 5 simplices (Mara, 1976); (2 3 5 8) is middle simplex
% See appendix A below if FURTHER GRID REFINEMENT IS REQUD (up to Z=5)
Sb      = [2 3 5 8; 1 2 3 5; 3 5 7 8; 2 5 6 8; 2 3 4 8];
nsx     = size(Sb,1);           % number of simplices
%-------------------------------------------------------------------------
% For further grid refinement,insert Appendix A here (see end of code)
%-------------------------------------------------------------------------
num_states= length(states);     % number of states = 8
Del_T   = Sb;                   % set of simplices X^tilde_Z=0

% % -- to save compuation time, we store the upper bound value functions for
% % each set of starting priors in a CSV file ahead of time 
% % Comment out the line below if you don't have this file
% ub_list=csvread('UBtable_i6_j6.csv'); % This file stores the upper bounds

%-------------------------------------------------------------------------
% STEP 1: Calculate & store approx vals-to-go in VF (uses simulation)
%-------------------------------------------------------------------------
tic;
VF_s         = fun_ap_sim_B; % calls the function
elapsedtime1 = elapsedtime1+toc;

%-------------------------------------------------------------------------
% STEP 2: Implement on ALL enumerated states (using simulation)
%-------------------------------------------------------------------------
tic;
approx_VF_sim = fun_ap_imp_sim_B;
elapsedtime2=elapsedtime2+toc;
disp('Iter. Exp proportion of succsesses Std Err')
ansr = [m approx_VF_sim./N]

% Display the computation times for each step at each iteration
disp ('Computation time for Step 1')
disp(elapsedtime1)
disp ('Computation time for Step 2')
disp(elapsedtime2)

% ***UNCOMMENT** the following line if you need to export the results to a CSV
% file (under the results folder)
dlmwrite('Results/VF_P20_T10_LJ.csv',ansr,'-append');
% end

% ----------END------------------------------

% %-----------APPENDIX A begins-------------------------------------------
% % ****UNCOMMENT BELOW if FURTHER GRID REFINEMENT IS REQUD (up to Z=5)***
% % -----the following calls the function "fun_grid_gen"
% statesz0         = points1; % H^tilde_Z=0 (Number of grid points at Z=0)
% Delz0            = Sb;      % X^tilde_Z=0 (Number of simplices at Z=0)
% Z=2;
% for z = 0:Z
%     % input to the function: fun_grid_gen
%     if z == 0
%         states_in = statesz0;
%         simplx_in = Delz0;
%     else
%         states_in = states_out;
%         simplx_in = simplx_out;
%     end
%     % output from the functionm: fun_grid_gen
%     [states_out simplx_out] = fun_grid_gen(states_in,simplx_in);
% end
% states = states_out;
% Sb = simplx_out;
% size(Sb)
% %-----------APPENDIX A ends--------------------------------------------
